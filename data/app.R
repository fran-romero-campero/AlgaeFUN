#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

## TODO
## Add button to upload files instead of inputting the genes in the text box
## Add buttom to clear text box with genes
## Add sanity check to test that the input genes follow the expected nomeclature
## Change width of table columns to fit mumber of genes (Enrichment too narrow)

## Add explanatory text before each table / figure
## Add download buttom for each table / figure

## Add nice logo and nice style to the web page

## Add GSEA analysis

## Add figures from KEGG enrichment

## To test the script:
# input <- list(microalgae = "otauri", pvalue = 0.05, analysis = "go", ontology = "BP", input_mode = "No")
# input <- list(microalgae = "otauri", pvalue = 0.05, analysis = "kegg", input_mode = "No")
# target.genes <- read.table(file="ostta/examples/no_iron_vs_iron_15H/activated/activated_genes.txt",as.is=T)[[1]]

# target.genes <- read.table(file="cre/examples/activated_genes.txt",as.is=T)[[1]]
# input <- list(microalgae = "creinhardtii", pvalue = 0.05, analysis = "go", ontology = "BP", input_mode = "No")

library(shinycssloaders)
library(shiny)
library(clusterProfiler)
library(pathview)

## Load microalgae annotation packages
library(org.Otauri.eg.db)
library(org.Creinhardtii.eg.db)
library(org.Dsalina.eg.db)

microalgae.names <- c("Ostreococcus tauri", "Chlamydomonas reinhardtii")
names(microalgae.names) <- c("otauri", "creinhardtii")

## Auxiliary functions
## Auxiliary function to compute enrichments
compute.enrichments <- function(gene.ratios, bg.ratios)
{
  gene.ratios.eval <- sapply(parse(text=gene.ratios),FUN = eval)
  bg.ratios.eval <- sapply(parse(text=bg.ratios),FUN = eval)
  enrichments <- round(x=gene.ratios.eval/bg.ratios.eval,digits = 2)
  enrichments.text <- paste(enrichments, " (", gene.ratios, "; ", bg.ratios, ")",sep="")

  return(enrichments.text)  
}

## Ostreococcus tauri gene link to ORCAE
## https://bioinformatics.psb.ugent.be/orcae/annotation/OsttaV2/current/ostta15g02520
ostta.gene.link <- function(gene.name)
{
  orcae.link <- paste0("https://bioinformatics.psb.ugent.be/orcae/annotation/OsttaV2/current/",gene.name)
  gene.link <- paste(c("<a href=\"",
                        orcae.link,
                        "\" target=\"_blank\">",
                        gene.name, "</a>"),
                      collapse="")
  return(gene.link)
}

## Gene link to Phytozome
## https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=0&searchText=Cre13.g569850&offset=0
phytozome.gene.link <- function(gene.name)
{
  phytozome.link <- paste(c("https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=0&searchText=",
                            gene.name,
                            "&offset=0"),collapse="")
  gene.link <- paste(c("<a href=\"",
                       phytozome.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## Gene Ontology term link
# http://amigo.geneontology.org/amigo/term/GO:0015979
go.link <- function(go.term)
{
  link <- paste0("http://amigo.geneontology.org/amigo/term/", go.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           go.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## KEGG pathway link
## https://www.genome.jp/kegg-bin/show_pathway?cre04136
kegg.pathway.link <- function(kegg.pathway)
{
  link <- paste0("https://www.genome.jp/kegg-bin/show_pathway?",kegg.pathway)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           kegg.pathway, "</a>"),
                         collapse = "")
  return(complete.link)
}

## KEGG module link
## https://www.genome.jp/kegg-bin/show_module?cre04136
kegg.module.link <- function(kegg.module)
{
  link <- paste0("https://www.genome.jp/kegg-bin/show_module?",kegg.module)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           kegg.module, "</a>"),
                         collapse = "")
  return(complete.link)
}




# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(#theme= "bootstrap.css",
  
  # Application title
  titlePanel("ALGAEFUN, microALGAE FUNctional annotation tool",windowTitle = "ALGAEFUN"),
  tags$br(),
  tags$br(),

  #Interface where the user can choose his/her preferencies, separated by columns
  fluidRow(
      column(width = 5,
             #Choose the target microalgae
        selectInput(inputId = "microalgae", label="Choose your favourite microalgae", 
                  choices=c("Ostreococcus tauri" = "otauri",
                            "Ostreococcus lucimarinus" = "olucimarinus",
                            "Coccomyxa subellipsoidea" = "csubellipsoidea",
                            "Dunaliella salina" = "dsalina",
                            "Chlamydomonas reinhardtii" = "creinhardtii",
                            "Bathycoccus prasinos" = "bathy")),
             #Choose a p-value
      numericInput(inputId = "pvalue", label= "Which will be your chosen p-value?", value= 0.05),
             #Choose the kind of analysis that you want us to execute 
      radioButtons(inputId = "analysis",
                   label="Choose your desirable analysis",
                   choices=c("GO terms enrichment" = "go",
                             "KEGG pathways enrichment analysis" = "kegg",
                             "Both" = "both"
                            )),
      
      conditionalPanel(condition= "input.analysis == 'go' || input.analysis == 'both'",
                       radioButtons(inputId = "ontology",
                                    label="Choose gene ontology:",
                                    choices = c("Biological process" = "BP",
                                                "Cellular Component" = "CC",
                                                "Mollecular Function" = "MF")))
      ),
      
      column(
        #The user can either insert your own background list or use ours. 
        radioButtons(inputId = "input_mode",
                   label = "Would you rather use your own background set?", 
                   choices = c("Yes", 
                               "No"),
                   selected = "No"),
        #This panel will only appear if the user chooses to use our background lists. 
      conditionalPanel(condition= "input.input_mode == 'No'",
      textAreaInput(inputId = "genes", label= "Insert a set of genes", width="250%", 
                    height = "400px",placeholder = "Insert set of genes",
                    value= "ostta11g02790
                    ostta10g02400
                    ostta09g02680
                    ostta10g03135
                    ostta17g00930"
      )
      ),
      #This panel wil only appear if the user wants to use his/her own background list. 
      conditionalPanel(condition= "input.input_mode == 'Yes'",
                       textAreaInput(inputId = "genes", label= "Insert a set of genes", width="250%", 
                                     height = "400px",placeholder = "Insert set of genes",
                                     value= "ostta11g02790
                    ostta10g02400
                    ostta09g02680
                    ostta10g03135
                    ostta17g00930"
                       ),
                       textAreaInput(inputId = "background", label= "Background set", width="250%", 
                                     height = "400px",placeholder = "Insert background list",
                                     value= "ostta11g02790
                    ostta10g02400
                    ostta09g02680
                    ostta10g03135
                    ostta17g00930"
                       )
                       #giving the chance to upload a file would also be interesting, 
                       #http://shiny.rstudio.com/gallery/upload-file.html
      ),
      
      actionButton(inputId = "go.button",label = "Have fun!", icon("send") ),
      width = 3)
      ),
    
    
    mainPanel(width = 13,
      # withSpinner(plotOutput(outputId = "plot", inline = TRUE, 
      #                        dblclick = "plot1_dblclick",
      #                        brush = brushOpts(
      #                                id = "plot1_brush",
      #                                resetOnNew = TRUE
      #                        ))),
      
      #Main panel containing the results organized in different tabs: GO map, Go terms data table, and 
      #KEGG pathway maps.
      tabsetPanel(type = "tabs",
                  tabPanel("GO map",
                           tags$br(),
                           tags$br(),
                           htmlOutput(outputId = "intro_go"),
                           tags$br(),
                           tags$br(),
                           dataTableOutput(outputId = "output_go_table"),
                           htmlOutput(outputId = "revigo"),
                           div(style= "overflow:scroll; height:500px;", 
                                       plotOutput(outputId = "go.plot", inline = TRUE)),
                           downloadButton(outputId= "downloadImage", "Get GO map"),
                           tags$br(),
                           tags$br(),
                           htmlOutput(outputId = "barplot_text"),
                           tags$br(),
                           tags$br(),
                           plotOutput(outputId = "bar.plot",inline=TRUE),
                           plotOutput(outputId = "dot.plot",inline=TRUE),
                           plotOutput(outputId = "emap.plot",inline=TRUE),
                           plotOutput(outputId = "cnet.plot",inline=TRUE)
                           ),
                           #align = "center"),
                  tabPanel("KEGG pathway", 
                           dataTableOutput(outputId = "output_pathway_table"),
                           plotOutput("keggpath"), 
                           dataTableOutput(outputId = "output_module_table"),
                           downloadButton(outputId= "downloadKEGGImage", "Get KEGG pathway image")),
                  tabPanel("Summary", dataTableOutput(outputId = "data"),
                           downloadButton(outputId= "downloadData", "Get GO terms of each gene"))#,
                           #htmlOutput(outputId = "revigo"))
                 )
              )
      ######tabpanels also can look like this: 
      #https://shiny.rstudio.com/gallery/navlistpanel-example.html
    ))

## Define server logic
server <- shinyServer(function(input, output) {
  observeEvent(input$go.button , {
    
    ## Select org.Db 
    if(input$microalgae == "otauri")
    {
      org.db <- org.Otauri.eg.db
    } else if (input$microalgae == "creinhardtii")
    {
      org.db <- org.Creinhardtii.eg.db
    }
    
    ## Extract genes
    target.genes <- as.vector(unlist(strsplit(input$genes, split="\n",
                                     fixed = TRUE)[1]))
    ## Select gene universe
    if(input$input_mode == "No")
    {
      gene.universe <- unique(select(org.db,columns = c("GO"),keys=keys(org.db,keytype = "GID"))[["GID"]])
      universe.text <- " default universe."
    } else 
    {
      gene.universe <- as.vector(unlist(strsplit(input$background, split="\n",
                                                 fixed = TRUE)[1]))
      universe.text <- paste(c(" custom universe (<i>", paste(gene.universe[1:3], collapse=" ")," </i> ...)."), collapse="")
    }

    ## GO term enirchment analysis
    if(input$analysis == "go" || input$analysis == "both")
    {
      
      ## Perform GO enrichment
      ## TODO q-value
      enrich.go <- enrichGO(gene          = target.genes,
                            universe      = gene.universe,
                            OrgDb         = org.db,
                            ont           = input$ontology,
                            pAdjustMethod = "BH",
                            pvalueCutoff  = input$pvalue,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE,
                            keyType = "GID")
      
      
      ## Generate ouput table
      enrich.go.result <- as.data.frame(enrich.go)
      #head(enrich.go.result)
      #enrich.go.result[1,]
      
      if(nrow(enrich.go.result) > 0)
      {
        ## Intro text for GO enrichment
        go.intro.text <- paste(c("This tab presents the results from the <b>GO enrichment analysis</b> 
                                      performed over the input target genes (<i>",
                                 paste(target.genes[1:3],collapse=" "),
                                 "</i> ...) from the microalgae <b> <i>", microalgae.names[input$microalgae],
                                 " </i> </b> with ", universe.text),collapse="") 
        output$intro_go <- renderText(expr = go.intro.text)
        
        ## GO term Description P-value Q-value Enrichment (SetRatio, BgRatio) Genes
        go.term.enrichments <- compute.enrichments(gene.ratios = enrich.go.result$GeneRatio,
                                                   bg.ratios = enrich.go.result$BgRatio)
        
        go.result.table <- data.frame(enrich.go.result$ID, enrich.go.result$Description,
                                      enrich.go.result$pvalue, enrich.go.result$qvalue,
                                      go.term.enrichments, 
                                      gsub(pattern = "/",replacement = " ",x = enrich.go.result$geneID),
                                      stringsAsFactors = FALSE)
        
        colnames(go.result.table) <- c("GO ID", "Description", "p-value", "q-value",
                                       "Enrichment (Target Ratio; BG Ration)","Genes")
        
        go.result.table.with.links <- go.result.table
        ## Add links to the genes
        genes.in.go.enrichment <- go.result.table$Genes
        
        if(input$microalgae == "otauri")
        {
          gene.link.function <- ostta.gene.link
        } else if(input$microalgae == "creinhardtii")
        {
          gene.link.function <- phytozome.gene.link
        }
        
        ## Add linkd to genes
        for(i in 1:length(genes.in.go.enrichment))
        {
          go.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(genes.in.go.enrichment[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }

        ## Add links to GO ids
        go.result.table.with.links[["GO ID"]] <- sapply(X = go.result.table.with.links[["GO ID"]], FUN = go.link)
        
        ## Output table with GO enrichment result
        output$output_go_table <- renderDataTable({
          go.result.table.with.links #go.result.table
        },escape=FALSE,options =list(pageLength = 5)) 
        
        ## Link to REVIGO 
        revigo.data <- paste(revigo.data <- apply(go.result.table[,c("GO ID", "q-value")], 1, paste, collapse = " "), collapse="\n")
        
        url1 <- a("here", href="#", onclick="document.revigoForm.submit();")
        url2 <- tags$form(
          name="revigoForm", action="http://revigo.irb.hr/", method="post", target="_blank", 
          tags$textarea(name="inputGoList", rows="1", cols="8", class="revigoText", 
                        style="visibility: hidden", revigo.data) 
        )
        
        output$revigo<- renderUI(
          tagList("Visualize output in REViGO", url1, url2)
        )
        
        ## GO plot
        output$go.plot <- renderPlot(
          width     = 940,
          height    = 900,
          res       = 120,
          expr = {
            goplot(enrich.go,showCategory = 10)
          })
        
        ## Barplot
        output$bar.plot <- renderPlot(
          width     = 870,
          height    = 600,
          res       = 120,
          expr = {
            barplot(enrich.go,drop=TRUE,showCategory = 10)
          })
        
        ## Dotplot
        output$dot.plot <- renderPlot(
          width     = 940,
          height    = 900,
          res       = 120,
          expr = {
            dotplot(enrich.go)
          })
        
        ##EMAP plot
        output$emap.plot <- renderPlot(
          width     = 940,
          height    = 900,
          res       = 120,
          expr = {
            emapplot(enrich.go)
          })
        
        ##CNET plot
        output$cnet.plot <- renderPlot(
          width     = 940,
          height    = 900,
          res       = 120,
          expr = {
            cnetplot(enrich.go)
          })
      } else
      {
        print("No GO term enrichment detected in the input set")
      }
      
    }
    
    ## KEGG pathways enrichment analysis
    if(input$analysis == "kegg"  || input$analysis == "both")
    {
      
      ## Update target genes and universe depending on the microalgae
      if(input$microalgae == "otauri")
      {
        target.genes <- paste0("OT_",target.genes)
        gene.universe <- paste0("OT_",gene.universe)
        organism.id <- "ota"
      } else if(input$microalgae == "creinhardtii")
      {
        cre.chlredraft.map <- select(org.Creinhardtii.eg.db,columns = c("CHLREDRAFT"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID"))
        head(cre.chlredraft.map)
        cre.ids <- cre.chlredraft.map$GID
        chlredraft.ids <- cre.chlredraft.map$CHLREDRAFT
        names(chlredraft.ids) <- cre.ids
        names(cre.ids) <- chlredraft.ids
        
        target.genes <- chlredraft.ids[target.genes]
        names(target.genes) <- NULL
        
        gene.universe <- chlredraft.ids[gene.universe]
        names(gene.universe) <- NULL
        
        organism.id <- "cre"
      }
      
      ## Compute KEGG pathway enrichment
      pathway.enrichment <- enrichKEGG(gene = target.genes, organism = organism.id, keyType = "kegg",
                                       universe = gene.universe,qvalueCutoff = input$pvalue)
      
      
      pathway.enrichment.result <- as.data.frame(pathway.enrichment)
      
      if(nrow(pathway.enrichment.result) > 0)
      {
        pathways.enrichment <- compute.enrichments(gene.ratios = pathway.enrichment.result$GeneRatio,
                                                   bg.ratios = pathway.enrichment.result$BgRatio)
        
        if (input$microalgae == "otauri")
        {
          kegg.enriched.genes <- gsub(pattern="OT_",replacement="",x=gsub(pattern = "/",replacement = " ",x = pathway.enrichment.result$geneID))
        } else if (input$microalgae == "creinhardtii")
        {
          kegg.enriched.genes <- pathway.enrichment.result$geneID
          for(i in 1:length(kegg.enriched.genes))
          {
            kegg.enriched.genes[i] <- paste(cre.ids[strsplit(kegg.enriched.genes[i],split="/")[[1]]],collapse=" ")
          }
        }

        pathways.result.table <- data.frame(pathway.enrichment.result$ID, pathway.enrichment.result$Description,
                                            pathway.enrichment.result$pvalue, pathway.enrichment.result$qvalue,
                                            pathways.enrichment, 
                                            kegg.enriched.genes,
                                            stringsAsFactors = FALSE)
        
        colnames(pathways.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                             "Enrichment (Target Ratio; BG Ration)","Genes")
        
        kegg.result.table.with.links <- pathways.result.table
        
        ## Add links to genes
        if(input$microalgae == "otauri")
        {
          gene.link.function <- ostta.gene.link
        } else if(input$microalgae == "creinhardtii")
        {
          gene.link.function <- phytozome.gene.link
        }
        
        for(i in 1:length(kegg.enriched.genes))
        {
          kegg.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(kegg.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }
        
        ## Add links to kegg pathways
        kegg.result.table.with.links[["KEGG ID"]] <- sapply(X=kegg.result.table.with.links[["KEGG ID"]],FUN = kegg.pathway.link)

        output$output_pathway_table <- renderDataTable({
          kegg.result.table.with.links
        },escape=FALSE,options =list(pageLength = 5)) 
      } else
      {
        print("No KEGG pathway enrichment detected in the input set")
      }
      
      modules.enrichment <- enrichMKEGG(gene = target.genes, organism = organism.id, keyType = "kegg")
      
      modules.enrichment.result <- as.data.frame(modules.enrichment)
      #modules.enrichment.result$
      
      if(nrow(modules.enrichment.result) > 0)
      {
        modules.enrichment <- compute.enrichments(gene.ratios = modules.enrichment.result$GeneRatio,
                                                   bg.ratios = modules.enrichment.result$BgRatio)
        
        if (input$microalgae == "otauri")
        {
          modules.enriched.genes <- gsub(pattern="OT_",replacement="",x=gsub(pattern = "/",replacement = " ",x = modules.enrichment.result$geneID))
        } else if (input$microalgae == "creinhardtii")
        {
          modules.enriched.genes <- modules.enrichment.result$geneID
          for(i in 1:length(modules.enriched.genes))
          {
            modules.enriched.genes[i] <- paste(cre.ids[strsplit(modules.enriched.genes[i],split="/")[[1]]],collapse=" ")
          }
        }
        
        
        
        modules.result.table <- data.frame(modules.enrichment.result$ID, modules.enrichment.result$Description,
                                           modules.enrichment.result$pvalue, modules.enrichment.result$qvalue,
                                           modules.enrichment, 
                                           modules.enriched.genes,
                                           stringsAsFactors = FALSE)
        
        colnames(modules.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                            "Enrichment (Target Ratio; BG Ration)","Genes")
        
        modules.result.table.with.links <- modules.result.table
        
        ## Add links to genes
        if(input$microalgae == "otauri")
        {
          gene.link.function <- ostta.gene.link
        } else if(input$microalgae == "creinhardtii")
        {
          gene.link.function <- phytozome.gene.link
        }
        
        for(i in 1:length(modules.enriched.genes))
        {
          modules.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(modules.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }
        
        ## Add links to kegg pathways
        modules.result.table.with.links[["KEGG ID"]] <- sapply(X=modules.result.table.with.links[["KEGG ID"]],FUN = kegg.module.link)
        
        output$output_module_table <- renderDataTable({
          modules.result.table.with.links
        },escape=FALSE,options =list(pageLength = 5)) 
        
        
        
      }
      
    }
    
    
    
    
    # p.valor <- input$pvalor
    # 
    # print(p.valor)
    # 
    # # p.valor<- 0.5
    # 
    # target.genes <- as.vector(unlist(
    #   strsplit(input$genes, split="\n",
    #            fixed = TRUE)[1]))
    # 
    # 
    # #Read annotation file
    # 
    # map.file <- paste(c("data/data_",input$mapfile,".map"), collapse = "")
    # 
    # geneID2GO <- readMappings(file=map.file)
    # 
    # #Set background
    # gene.names <- attributes(geneID2GO)[[1]]
    # gene.background <- rep(1, length(gene.names))
    # names(gene.background) <- gene.names
    # 
    # gene.background[target.genes] <- 0
    # 
    # #Every target gene in background = 0
    # gene.selec <- function(gene.list)
    # {
    #   return(gene.list == 0)
    # }
    # #Make the object topGOdata
    # sampleGOdata <- new("topGOdata",
    #                     description = "Microalgae session", ontology = "BP",
    #                     allGenes = gene.background, geneSel = gene.selec,
    #                     nodeSize = 10,
    #                     annot = annFUN.gene2GO, gene2GO = geneID2GO)
    # 
    # resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    # 
    # #Create the table where you'll write the most relevant 100 GO terms 
    # allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
    #                    ranksOf = "classicFisher", topNodes = 100, numChar=100)
    # 
    # #Save p-value column and create a vector where you'll save the corrected p-values.
    # all.pvalues <- as.vector(allRes[["classicFisher"]])
    # 
    # corrected.all.pvalues <- vector(length=length(all.pvalues), mode="numeric")
    # #Change value of the smallest p-value to 0, because they are so small that
    # #they can't be seen as numeric
    # for (i in 1:length(all.pvalues))
    # {
    #   if (is.na(as.numeric(all.pvalues[i])))
    #   {
    #     corrected.all.pvalues[i] <- 0
    #   }else
    #   {
    #     corrected.all.pvalues[i] <- as.numeric(all.pvalues[i])
    #   }
    # }
    # #Take just the p-values that are smaller than our p-value pre-selected.
    # filter.Res <- allRes[corrected.all.pvalues < p.valor,]
    # res.with.genes <- data.frame(filter.Res,genes=vector(length=nrow(filter.Res),mode="character"))
    # 
    # # #Loop that will save the GO terms of the target genes that have a significant
    # # #p-value.
    # #Save the GO terms related to each target gene
    # target.genes.annotation <- geneID2GO[target.genes]
    # 
    # #Save 100 GO terms selected
    # significative.GO.terms <- as.vector(filter.Res[["GO.ID"]])
    # 
    # #Create the list of all the possible GO terms 
    # offspring <- as.list(GOBPOFFSPRING)
    # 
    # #Create the empty vector where you'll save the genes in every GOterm selected.
    # genes.for.each.GO <- vector(length=length(significative.GO.terms),mode="character")
    # 
    # 
    # for (i in 1:length(significative.GO.terms))
    # {
    #   #for each GO terms of the 100 GO terms selected, get all the GO
    #   #terms related to it
    #   GO.genes <- vector(mode="character")
    #   k <- 1
    #   go.offspring <- offspring[[significative.GO.terms[i]]]
    #   for (j in 1:length(target.genes.annotation))
    #   {
    #     #for every GO term related to each target gene
    #     if(length(intersect(go.offspring,target.genes.annotation[[j]]))!=0)
    #     {
    #       #if there are common GO terms between the GO terms related to each
    #       #target gene and the GO terms related to each of the 100 GO terms selected,
    #       #save the gene in GO.genes vector that you created before.
    #       GO.genes[k] <- names(target.genes.annotation[j])
    #       k <- k+1
    #     }
    #   }
    #   GO.genes.string <- paste(unique(GO.genes), collapse=" ")
    #   
    #   genes.for.each.GO[i] <- GO.genes.string
    #   print(i)
    # }
    # 
    # res.with.genes[["genes"]] <- genes.for.each.GO
    # 
    # #Generate GOmap   
    # output$plot<- 
    #   renderPlot( 
    #     width     = 4700,
    #     height    = 4500,
    #     res       = 260,
    #     {
    #       
    #       ## Create the final image
    #       showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
    #       
    #     })
    # 
    # 
    # ## Download the final image
    # output$downloadImage<- downloadHandler(
    #   filename= "GO_map.png",
    #   content= function(file) {
    #     png(file, width     = 10,
    #         height    = 10,
    #         units     = "in",
    #         res       = 600)
    #     showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
    #     dev.off()
    #   })  
    # 
    # #Function to generate the links that will take you to a web about each gene.
    # createGeneLink<- function (gene.name)
    # {
    #   paste(c("<a href=\"https://www.ncbi.nlm.nih.gov/gene/?term=", 
    #           gene.name,
    #           "\" target=\"_blank\">",
    #           gene.name,
    #           "</a>"
    #   ), collapse= "")
    # }
    # #Function to generate the links that will take you to a web about each GO term.
    # 
    # createGOLink<- function (go.term)
    # {
    #   
    #   paste(c("<a href=\"http://amigo.geneontology.org/amigo/term/", 
    #           go.term,
    #           "\" target=\"_blank\">",
    #           go.term,
    #           "</a>"
    #   ), collapse= "")
    # }
    # #make a copy of res.with.genes before changing the data that it contains 
    # res.with.genes.1<- cbind (res.with.genes)
    # output$downloadData<- downloadHandler(
    #   filename= c('GO_terms.txt'),
    #   content= function(file) {
    #     write.table(res.with.genes.1, file, row.names = TRUE)
    #   })
    # 
    # #Use the previous functions before printing the final data table on screen      
    # res.with.genes$GO.ID <- sapply(res.with.genes$GO.ID, createGOLink)
    # 
    # for (i in 1:length(res.with.genes$genes))
    # {
    #   vector.genes <-strsplit(res.with.genes$genes[[i]], split=" ")[[1]]
    #   res.with.genes$genes[[i]] <- paste(sapply(vector.genes, createGeneLink),collapse= " ")
    # }
    # 
    # #Print the final data table on screen
    # output$data<- renderDataTable({
    #   res.with.genes
    # }, escape= FALSE)
    
    
    #     <a href="#" onclick="document.revigoForm.submit();"> Visualize output in REViGO </a><br>
    #       </p><form name="revigoForm" action="http://revigo.irb.hr/" method="post" target="_blank">
    #       <textarea name="inputGoList" rows="1" cols="8" class="revigoText" style="visibility: hidden">% List genereted using CircadiaNET
    #     % http://viridiplantae.ibvf.csic.es/circadiaNet/"),collapse = ""
    #     GO:0006811	0.0021
    #     GO:0006810	0.0035
    #     GO:0051234	0.0039
    #     GO:0051179	0.0046
    #     GO:0030001	0.0075
    #     GO:0006812	0.0081
    #     GO:0030163	0.0122
    #     GO:0006366	0.0159
    #     GO:0006508	0.0205
    #     GO:0009057	0.0213
    #     GO:1901575	0.0421
    #     </textarea>
    #     </form></p>
    # 
    #Link to revigo 
    # j<-1
    # revigo.data<-c()
    # #Generate the input data in the correct form: "GO p.value"
    # for (j in 1:length(res.with.genes.1$GO.ID))
    # {
    #   each.element<-paste(res.with.genes.1$GO.ID[j], res.with.genes.1$classicFisher[j], collapse = " ")
    #   revigo.data[j]<-paste(each.element,collapse = "")
    # }  
    # revigo.data<-paste(revigo.data, collapse="\n")
    # 
    # 
    # url1<- a("here", href="#", onclick="document.revigoForm.submit();")
    # url2<- tags$form(
    #   name="revigoForm", action="http://revigo.irb.hr/", method="post", target="_blank", 
    #   tags$textarea(name="inputGoList", rows="1", cols="8", class="revigoText", 
    #                 style="visibility: hidden", revigo.data) 
    # )
    # 
    # 
    # output$revigo<- renderUI(
    #   tagList("Visualize output in REViGO", url1, url2)
    # )
    
    
    #                    
    #                            
    
  })
  
  
})

# Run the application 
shinyApp(ui = ui, server = server)