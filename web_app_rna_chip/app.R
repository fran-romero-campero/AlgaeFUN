#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

## TODO
## Change width of table columns to fit mumber of genes (Enrichment too narrow)
## Add download buttom for each table / figure
## Add Warning messages when no enrichment is detected
## Add nice logo and nice style to the web page
## Add Functional annotation of genomic locations

## Issue with KEGG ids for Ngaditana
## expected: NGA_0076200,NGA_2122800,NGA_0378100,NGA_0380800,NGA_0239100,NGA_0126801

## To test the script:
# input <- list(microalgae = "otauri", pvalue = 0.05, analysis = "go", ontology = "BP", input_mode = "No")
# input <- list(microalgae = "otauri", pvalue = 0.05, analysis = "kegg", input_mode = "No")
# input <- list(microalgae = "ptricornutum", pvalue = 0.05, analysis = "kegg", input_mode = "No")
# input <- list(microalgae = "ngaditana", pvalue = 0.05, analysis = "kegg", input_mode = "No")
# input <- list(microalgae = "knitens", pvalue = 0.05, analysis = "kegg", input_mode = "No")
# input <- list(microalgae = "csubellipsoidea", pvalue = 0.05, analysis = "go", input_mode = "No")

# target.genes <- read.table(file="example_files/example_otauri.txt",as.is=T)[[1]]
# target.genes <- read.table(file="cre/examples/activated_genes.txt",as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_vcarteri.txt",as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_ptricornutum.txt", as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_ngaditana_1.txt", as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_knitens.txt", as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_csubellipsoidea.txt", as.is=T)[[1]]

# input <- list(microalgae = "creinhardtii", pvalue = 0.05, analysis = "go", ontology = "BP", input_mode = "No")

library(shinycssloaders)
library(shiny)
library(clusterProfiler)
library(pathview)
library(ChIPseeker)
library(ChIPpeakAnno)
library(rtracklayer)
library(seqinr)

## Load microalgae annotation packages
library(org.Otauri.eg.db)
library(org.Creinhardtii.eg.db)
library(org.Dsalina.eg.db)
library(org.Vcarteri.eg.db)
library(org.Ptricornutum.eg.db)
library(org.Ngaditana.eg.db)
library(org.Knitens.eg.db)
library(org.Csubellipsoidea.eg.db)

## Load microalgae genome annotation packages
library(TxDb.Otauri.JGI)
library(TxDb.Creinhardtii.Phytozome)
#
#
library(TxDb.Ptricornutum.Ensembl.Protists)

microalgae.names <- c("Ostreococcus tauri", 
                      "Chlamydomonas reinhardtii", 
                      "Dunaliella salina",
                      "Volvox carteri",
                      "Phaeodactylum tricornutum",
                      "Nannochloropsis gaditana",
                      "Klebsormidium nitens",
                      "Coccomyxa subellipsoidea",
                      "Bathycoccus prasinos")
names(microalgae.names) <- c("otauri", 
                             "creinhardtii", 
                             "dsalina", 
                             "vcarteri",
                             "ptricornutum", 
                             "ngaditana",
                             "knitens",
                             "csubellipsoidea",
                             "bprasinos")

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

#Bathy gene link to ORCAE
bathy.gene.link <- function(gene.name)
{
  orcae.link <- paste0("https://bioinformatics.psb.ugent.be/orcae/annotation/Bathy/current/",gene.name)
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

## Phaeodactylum tricornutum link to ENSEMBL PROTISTS
## https://protists.ensembl.org/Phaeodactylum_tricornutum/Gene/Summary?g=Phatr3_EG00535
phaeodactylum.gene.link <- function(gene.name)
{
  phatri.link <- paste(c("https://protists.ensembl.org/Phaeodactylum_tricornutum/Gene/Summary?g=",
                          gene.name),collapse="")
  gene.link <- paste(c("<a href=\"",
                       phatri.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## Nannochloropsis gaditana link to CRIBI Genomics
#http://www.nannochloropsis.org/gene/Naga_100001g4
ngaditana.gene.link <- function(gene.name)
{
  naga.link <- paste(c("http://www.nannochloropsis.org/gene/",
                         gene.name),collapse="")
  gene.link <- paste(c("<a href=\"",
                       naga.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## Klebsormidium nitens link to CGA Klebsormidium
#http://genome.annotation.jp/klebsormidium/nies-2285/genes/kfl01813_0010
knitens.gene.link <- function(gene.name)
{
  knitens.link <- paste(c("http://genome.annotation.jp/klebsormidium/nies-2285/genes/",
                       gene.name),collapse="")
  gene.link <- paste(c("<a href=\"",
                       knitens.link,
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

## Extract annotation from the result of chipseeker
extract.annotation <- function(annotation.str)
{
  return(strsplit(x = annotation.str,split=" \\(")[[1]][1])
}

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(#theme= "bootstrap.css",
  
  # Application title
  titlePanel("ALGAEFUN, microALGAE FUNctional annotation tool",windowTitle = "ALGAEFUN"),
  tags$br(),

  tags$p("Welcome to", tags$b("ALGAEFUN"),", the microalgae functional annotation tool. ALGAEFUN
         is a web-based tool and database for the functional annotation of gene
         sets and genomic locations from a wide collection of microalgae that
         include Chlamydomonas reinhardtii, Ostreococcus tauri, Phaeodactylum tricornutum
         and Nannochlorpsis gaditana."), 
  tags$br(), 
  
  radioButtons(inputId = "go_chip", width="100%",selected="",
               label="Choose between annotating genes sets geneared, 
                     for instance, from an RNA-seq analysis or genomics regions
                     obtained, for instance, from a Chip-seq analysis:",
               choices=c("Gene set annotation" = "gene_sets",
                         "Genomic regions annotations" = "genomic_regions"
               )),
  
  #Interface where the user can choose his/her preferencies, separated by columns
  fluidRow(
      column(width = 5,

        #Choose the target microalgae
        selectInput(inputId = "microalgae", label="Choose your favourite microalgae", 
                  choices=c("Ostreococcus tauri" = "otauri",
                            "Chlamydomonas reinhardtii" = "creinhardtii",
                            "Dunaliella salina" = "dsalina",
                            "Volvox Carteri" = "vcarteri",
                            "Phaeodactylum tricornutum" = "ptricornutum",
                            "Nannochloropsis gaditana" = "ngaditana",
                            "Ostreococcus lucimarinus" = "olucimarinus",
                            "Coccomyxa subellipsoidea" = "csubellipsoidea",
                            "Bathycoccus prasinos" = "bathy",
                            "Klebsormidium nitens" = "knitens")),


      #Choose a p-value
      conditionalPanel(condition = "input.go_chip == 'gene_sets'",
        numericInput(inputId = "pvalue", 
                     label= "Which will be your chosen p-value?", 
                     value= 0.05)),
 
      #Choose the kind of analysis that you want us to execute 
      conditionalPanel(condition = "input.go_chip == 'gene_sets'",
        radioButtons(inputId = "analysis",
                   label="Choose your desirable analysis",
                   choices=c("GO terms enrichment" = "go",
                             "KEGG pathways enrichment analysis" = "kegg",
                             "Both" = "both"
                            ))),
      
      conditionalPanel(condition= "(input.analysis == 'go' || input.analysis == 'both') &&
                                   input.go_chip == 'gene_sets'",
                       radioButtons(inputId = "ontology",
                                    label="Choose gene ontology:",
                                    choices = c("Biological process" = "BP",
                                                "Cellular Component" = "CC",
                                                "Mollecular Function" = "MF"))),
       
      #Choose a promoter length
      conditionalPanel(condition = "input.go_chip == 'genomic_regions'",
                       sliderInput(inputId = "promoter_length", 
                                    label= "Choose the distance in base pairs around the 
                                    Transcriptional Start Site to determine gene promoters", 
                                    min=100, max=2000,value=1000,step=100)),
      
      #Choose a promoter length
      conditionalPanel(condition = "input.go_chip == 'genomic_regions'",
                       checkboxGroupInput(inputId = "selected_genomic_features", 
                                   label= "A gene will be considered as a target of an input genomic locus
                                   when it overlaps one of the following gene parts:",
                                   choices = c("Promoter","5' UTR","Exon","Intron","3' UTR")))
      
      ),
      
      column( width = 7,
        #Choose the kind of analysis that you want us to execute 
        # tags$b("Choose between annotating genes sets geneared, 
        #              for instance, from an RNA-seq analysis or genomics regions
        #              obtained, for instance, from a Chip-seq analysis:"),
        
        conditionalPanel(condition = "input.go_chip == 'gene_sets'",
          #The user can either insert his/her own background list or use ours. 
          radioButtons(inputId = "input_mode",
                       label = "Would you rather use your own background set?", 
                       choices = c("Yes", 
                                   "No"),
                        selected = "No"),

          #This panel will only appear if the user chooses to use our background lists. 
          actionButton(inputId = "example_genes",label = "Example"),
          textAreaInput(inputId = "genes", label= "Insert a set of genes", width="200%", 
                        height = "200px",placeholder = "Insert set of genes",
                        value= "
                        ostta11g02790
                        ostta10g02400
                        ostta09g02680
                        ostta10g03135
                        ostta17g00930"
          ),
      
        actionButton(inputId = "clear_gene_set",label = "Clear"),
        fileInput(inputId = "gene_set_file",label = "Choose File with Gene Set to Upload")),

        #This panel will only appear if the user wants to use his/her own background list. 
        conditionalPanel(condition = "input.input_mode == 'Yes' && input.go_chip == 'gene_sets'",
                         textAreaInput(inputId = "background", label= "Background set", width="200%", 
                                       height = "100px",placeholder = "Insert background list",
                                       value= "ostta11g02790
                                               ostta10g02400
                                               ostta09g02680
                                               ostta10g03135
                                               ostta17g00930"
                       ),
                       
                       actionButton(inputId = "clear_universe_set",label = "Clear"),
                       fileInput(inputId = "gene_universe_file",label = "Choose File with Custom Gene Universe to Upload")
                       
      ),
      
      conditionalPanel(condition = "input.go_chip == 'genomic_regions'",
        #This panel will only appear if the user chooses to use our background lists. 
        actionButton(inputId = "example_genomic_regions",label = "Example"),
        textAreaInput(inputId = "genomic_regions", 
                      label= "Insert a set of genomic regions", 
                      width="200%", height = "200px", 
                      placeholder = "Insert set of genomic regions",
                      value= "
                              chromosome_1	23846	24380
                              chromosome_1	30768	31215
                              chromosome_1	41249	41658
                              chromosome_1	41846	42231
                              chromosome_1	53798	54602
                              chromosome_1	60199	60649
                              chromosome_1	64652	65507
                              chromosome_1	69907	70578
                              chromosome_1	88474	89012
                              chromosome_1	89532	90424
                              chromosome_1	110158	110801
                              chromosome_1	113590	113972
                              chromosome_1	115923	117189
                              chromosome_1	121463	121887
                              chromosome_1	132220	132619
                              chromosome_1	137155	138578
                              chromosome_1	151990	152123
                              chromosome_1	152429	153198
                              chromosome_1	166710	167389
                              chromosome_1	178641	180196"
                       ),
                       
        actionButton(inputId = "clear_genomic_regions",label = "Clear"),
        fileInput(inputId = "genomic_regions_file",label = "Choose File with the Genomic Regions to Upload",width = "100%"),
        actionButton(inputId = "genomic_button",label = "Have fun!", icon("send") )                       
                       ),
      conditionalPanel(condition = "input.go_chip == 'gene_sets'",
                       actionButton(inputId = "go.button",label = "Have fun!", icon("send") )                       
      )
    )
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
      conditionalPanel(condition = "input.go_chip == 'gene_sets'",
            tabsetPanel(type = "tabs",
                  tabPanel("GO map",
                           tags$br(), tags$br(),
                           htmlOutput(outputId = "gene_sanity_go"),
                           htmlOutput(outputId = "wrong_genes_go"),
                           htmlOutput(outputId = "intro_go"),
                           htmlOutput(outputId = "textGOTable"),
                           tags$br(), tags$br(),
                           dataTableOutput(outputId = "output_go_table"),
                           htmlOutput(outputId = "revigo"),
                           downloadButton(outputId= "downloadData", "Get GO terms of each gene"),
                           tags$br(),
                           htmlOutput(outputId = "go_graph"),
                           tags$br(), tags$br(),
                           div(style= "overflow:scroll; height:500px;", 
                                       plotOutput(outputId = "go.plot", inline = TRUE)),
                           downloadButton(outputId= "downloadImage", "Get this GO plot"),
                           tags$br(), tags$br(),
                           htmlOutput(outputId = "barplot_text"),
                           tags$br(),
                           plotOutput(outputId = "bar.plot",inline=TRUE),
                           downloadButton(outputId= "downloadbarplot", "Get this plot"),
                           tags$br(), tags$br(),
                           htmlOutput(outputId = "dotplot_text"),
                           tags$br(),
                           plotOutput(outputId = "dot.plot",inline=TRUE),
                           downloadButton(outputId= "downloadotplot", "Get this plot"),
                           tags$br(), tags$br(),
                           htmlOutput(outputId = "emapplot_text"),
                           tags$br(),
                           plotOutput(outputId = "emap.plot",inline=TRUE),
                           downloadButton(outputId= "downloademapplot", "Get this plot"),
                           tags$br(), tags$br(),
                           htmlOutput(outputId = "cnetplot_text"),
                           tags$br(),
                           plotOutput(outputId = "cnet.plot",inline=TRUE),
                           downloadButton(outputId= "downloadcnetplot", "Get this plot")
                           
                           ),
                           #align = "center"),
                  tabPanel("KEGG pathway", 
                           tags$br(), tags$br(),
                           htmlOutput(outputId = "gene_sanity_kegg"),
                           htmlOutput(outputId = "wrong_genes_kegg"),
                           htmlOutput(outputId = "intro_kegg"),
                           htmlOutput(outputId = "textKEGGTable"),
                           tags$br(), tags$br(),
                           htmlOutput(outputId = "kegg_pathway_table_text"),
                           tags$br(),
                           dataTableOutput(outputId = "output_pathway_table"),
                           br(), br(), br(),
                           htmlOutput(outputId = "textKEGGImage"),
                           br(), br(),
                           uiOutput(outputId = "kegg_selectize"),
                           imageOutput("kegg_image"),
                           br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                           br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                           br(), br(), br(), br(), br(),
                           htmlOutput(outputId = "text_module_kegg"),
                           br(), br(),
                           dataTableOutput(outputId = "output_module_table"),
                           uiOutput(outputId = "kegg_module_selectize"),
                           imageOutput("kegg_module_image"),
                           downloadButton(outputId= "downloadKEGGImage", "Get KEGG pathway image"))#,
                  # tabPanel("Summary", dataTableOutput(outputId = "data"),
                  #          downloadButton(outputId= "downloadData", "Get GO terms of each gene"))#,
                           #htmlOutput(outputId = "revigo"))
                 )
              ),
      ######tabpanels also can look like this: 
      #https://shiny.rstudio.com/gallery/navlistpanel-example.html
      
      conditionalPanel(condition = "input.go_chip == 'genomic_regions'",
                       plotOutput(outputId = "annotation.pie.chart",inline=TRUE),
                       plotOutput(outputId = "distance.to.tss",inline=TRUE),
                       uiOutput(outputId = "annotated_genes")
      )
      
      
      )
  

))

## Define server logic
server <- shinyServer(function(input, output, session) {

  ## Clear content of gene set text area
  observeEvent(input$clear_gene_set, {
    updateTextAreaInput(session=session, inputId = "genes",value = "")
  })
  
  ## Clear content of genomic regions set text area
  observeEvent(input$clear_genomic_regions, {
    updateTextAreaInput(session=session, inputId = "genomic_regions",value = "")
  })

  ## Clear content of universe set text area
  observeEvent(input$clear_universe_set, {
    updateTextAreaInput(session=session, inputId = "background",value = "")
  })
  
  ## Add an example of gene set to text area
  observeEvent(input$example_genes, {
    example.file <- paste(c("example_files/example_",input$microalgae,".txt"),collapse="")
    example.genes <- read.table(file = example.file,header = F,as.is = T)[[1]]
    updateTextAreaInput(session=session, inputId = "genes",value = paste(example.genes,collapse="\n"))
  })
  

  ## Add an example of gene set to text area
  observeEvent(input$example_genomic_regions, {
    example.file <- paste(c("example_files/example_genomic_regions_",input$microalgae,".txt"),collapse="")
    example.genomic.regions <- read.table(file = example.file,header = F,as.is = T)
    example.text <- NULL
    for(i in 1:nrow(example.genomic.regions))
    {
      example.text <- paste(example.text,paste(example.genomic.regions[i,],collapse = "\t"),sep="\n")
    }
    print(example.text)
    updateTextAreaInput(session=session, inputId = "genomic_regions",value = example.text)
  })
  
  ## Actions to perform after click the go button
  observeEvent(input$go.button , {
    
    ## Select org.Db 
    if(input$microalgae == "otauri")
    {
      org.db <- org.Otauri.eg.db
      microalgae.genes <- read.table(file = "universe/otauri_universe.txt",as.is = T)[[1]]
    } else if (input$microalgae == "creinhardtii")
    {
      org.db <- org.Creinhardtii.eg.db
      microalgae.genes <- read.table(file = "universe/cre_universe.txt",as.is = T)[[1]]
    } else if (input$microalgae == "dsalina")
    {
      org.db <- org.Dsalina.eg.db
      microalgae.genes <- read.table(file = "universe/dusal_universe.txt",as.is = T)[[1]]
    } else if (input$microalgae == "vcarteri")
    {
      org.db <- org.Vcarteri.eg.db
      microalgae.genes <- read.table(file = "universe/vocar_universe.txt",as.is = T)[[1]]
    } else if (input$microalgae == "ptricornutum")
    {
      org.db <- org.Ptricornutum.eg.db
      microalgae.genes <- read.table(file = "universe/phatri_universe.txt",as.is = T)[[1]]
    } else if (input$microalgae == "ngaditana")
    {
      org.db <- org.Ngaditana.eg.db
      microalgae.genes <- read.table(file = "universe/naga_universe.txt",as.is = T)[[1]]
    } else if (input$microalgae == "knitens")
    {
      org.db <- org.Knitens.eg.db
      microalgae.genes <- read.table(file = "universe/klebsor_universe.txt",as.is = T)[[1]]
    }else if (input$microalgae == "bathy")
    {
      org.db <- org.Bprasinos.eg.db
      microalgae.genes <- read.table(file = "universe/bathy_universe.txt",as.is = T)[[1]]
    }
    
    ## Extract genes from text box or uploaded file
    if(is.null(input$gene_set_file))
    {
      target.genes <- as.vector(unlist(strsplit(input$genes, split="\n",
                                                fixed = TRUE)[1]))
    } else
    {
      target.genes <- read.table(file=input$gene_set_file$datapath, header = F,as.is = TRUE)[[1]]
    }
    
    ## Select gene universe
    if(input$input_mode == "No")
    {
      gene.universe <- unique(select(org.db,columns = c("GO"),keys=keys(org.db,keytype = "GID"))[["GID"]])
      universe.text <- " default universe."
    } else 
    {
      if(is.null(input$gene_universe_file))
      {
        gene.universe <- as.vector(unlist(strsplit(input$background, split="\n",
                                                   fixed = TRUE)[1]))
      } else
      {
        gene.universe <- read.table(file=input$gene_universe_file$datapath, header = F,as.is = TRUE)[[1]]
      }
      
      universe.text <- paste(c(" custom universe (<i>", paste(gene.universe[1:3], collapse=" ")," </i> ...)."), collapse="")
    }

    ## Sanity test
    wrong.genes <- setdiff(target.genes,microalgae.genes)
    percent.wrong.genes <- length(wrong.genes) / length(target.genes)

    if(length(wrong.genes) > 0)
    {
      gene.sanity.text <- paste(c("<b> ERROR: The following genes have been detected in your target gene set that do not match any gene
        in our annotation for ", input$microalgae, ". You may consider removing them. Alternatively, check
        that the right microalgae has been selected and that your gene names follow the right nomenclature. Click on the Example button
        on top of the text area to input the target gene set to obtain an example of the gene nomenclature expected. If your target gene set
        has been generated from an RNA-seq experiment please consider using our tool MARACAS to obtain it or use our reference genome sequence
        and annotation files for your own analysis."),collapse="")

      output$gene_sanity_go <-renderText(expr = gene.sanity.text)
      output$wrong_genes_go <- renderText(expr = paste(wrong.genes,collapse="\n"))
      output$gene_sanity_kegg <-renderText(expr = gene.sanity.text)
      output$wrong_genes_kegg <- renderText(expr = paste(wrong.genes,collapse="\n"))
    }

    ## GO term enirchment analysis
    if((input$analysis == "go" || input$analysis == "both") && (length(wrong.genes) == 0))
    {
      ## Intro text for GO enrichment
      go.intro.text <- paste(c("This tab presents the results from the <b>GO enrichment analysis</b> 
                                      performed over the input target genes (<i>",
                               paste(target.genes[1:3],collapse=" "),
                               "</i> ...) from the microalgae <b> <i>", microalgae.names[input$microalgae],
                               " </i> </b> with ", universe.text),collapse="") 
      
      output$intro_go <- renderText(expr = go.intro.text)

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

      if(nrow(enrich.go.result) > 0)
      {
        
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
        } else if(input$microalgae == "creinhardtii" | input$microalgae == "dsalina" | input$microalgae == "vcarteri")
        {
          gene.link.function <- phytozome.gene.link
        } else if(input$microalgae == "ptricornutum")
        {
          gene.link.function <- phaeodactylum.gene.link
        } else if(input$microalgae == "ngaditana")
        {
          gene.link.function <- ngaditana.gene.link
        } else if (input$microalgae == "knitens")
        {
          gene.link.function <- knitens.gene.link
        }
        
        ## Add link to genes
        for(i in 1:length(genes.in.go.enrichment))
        {
          go.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(genes.in.go.enrichment[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }

        ## Add links to GO ids
        go.result.table.with.links[["GO ID"]] <- sapply(X = go.result.table.with.links[["GO ID"]], FUN = go.link)
        
        
        ## Introductory text for GO enrichment table
        go.table.text <- "The table below summarizes the result of the GO term
enrichment analysis. Each row represents a GO term significantly enriched in the target
gene set with respect to the selected gene universe. The first column represents the GO term
identifier. The second column contains a human readable description. For more details on the 
corresponding GO term, click on the identifier in the first column. The third and fourth 
column presents the p-value and q-value (adjusted p-value or FDR) capturing the level
of significance. The fifth column displays the corresponding enrichment value E (m/n; M/N) where
n is the number of genes with annotation from the target set, N is the number of genes with
annotation from the gene universe, m is the number of genes from the target set annotated with the
corresponding GO term and M is the number of genes from the gene universe annotated with
the GO term associated with the corresponding row. The enrichment is then computed as
E = (m/n) / (M/N). Finally, the last column, contains the genes from the target set
annotated with the GO term represented in the corresponding row." 
        
        output$textGOTable <- renderText(expr = go.table.text)

        ## Output table with GO enrichment result
        output$output_go_table <- renderDataTable({
          go.result.table.with.links #go.result.table
        },escape=FALSE,options =list(pageLength = 5)) 
        ## Download result
        output$downloadData<- downloadHandler(
          filename= function() {
              paste("godata-",microalgae.names[input$microalgae] , ".csv", sep="")
            },
          content= function(file) {
            write.csv(go.result.table,
                      file,
                      row.names=TRUE
                      )
          })
        
        ## Link to REVIGO 
        # revigo.data <- paste(revigo.data <- apply(go.result.table[,c("GO ID", "q-value")], 1, paste, collapse = " "), collapse="\n")
        # 
        # url1 <- a("here", href="#", onclick="document.revigoForm.submit();")
        # url2 <- tags$form(
        #   name="revigoForm", action="http://revigo.irb.hr/", method="post", target="_blank",
        #   tags$textarea(name="inputGoList", rows="1", cols="8", class="revigoText",
        #                 style="visibility: hidden", revigo.data)
        # )
        # 
        # output$revigo<- renderUI(
        #   tagList("The enriched GO terms above may be redundant. Visualize these results in REViGO in order to remove redundancy. Click", url1, url2)
        # )
        # 
        
        go.graph.text <- "The following acyclic graph represents the GO term enrichment
        in the target gene set. Each node stands for a GO term. The color of each node
        indicates the level of significance from grey, non-significant, to intense red,
        highly significant. An arrow is drawn from GO term A to GO term B when A is a more
        general GO term than B or B is more specific than A."
        
        output$go_graph <- renderText(expr = go.graph.text)
        
        ## GO plot
        output$go.plot <- renderPlot(
          width     = 940,
          height    = 900,
          res       = 120,
          expr = {
            goplot(enrich.go,showCategory = 10)
          })
        
        output$barplot_text <- renderText("In the following barplot each bar represents a significantly enriched 
        GO term. The length of the bar corresponds to the number of genes in the
        target set annotated with the given GO term. The bar color captures the level
        of significance from blue, less significant, to red, more significant.")
        
        ## Barplot
        output$bar.plot <- renderPlot(
          width     = 870,
          height    = 600,
          res       = 120,
          expr = {
            barplot(enrich.go,drop=TRUE,showCategory = 10)
          })
        
        output$dotplot_text <- renderText("In the following dotplot each dot represents a significantly enriched 
        GO term. The x-position of the dot corresponds to the ratio between the number of genes annotated with the
corresponding GO term and the total number of annotated genes in the target set. The dot color captures the level
        of significance from blue, less significant, to red, more significant.")
        
        ## Dotplot
        output$dot.plot <- renderPlot(
          width     = 870,
          height    = 600,
          res       = 120,
          expr = {
            dotplot(enrich.go)
          })
        
        output$emapplot_text <- renderText("The following figure consists of an enrichment map where nodes represent enriched GO terms. The
        size of a node is proportional to the number of genes annotated with the corresponding GO term in the target set.
The node colors represent the level of significance from less signficant in blue to more significant in red. Edges are drawn
between two nodes when the corresponding GO terms are semantically related.")
        
        ##EMAP plot
        output$emap.plot <- renderPlot(
          width     = 870,
          height    = 600,
          res       = 120,
          expr = {
            emapplot(enrich.go)
          })
        
        output$cnetplot_text <- renderText("The following figure corresponds to a gene-concept network. The beige
nodes represents GO terms and the grey nodes genes. An edge is drawn from a gene to a GO term when the gene is annotated
with the corresponding gene. The size of nodes representing GO terms is proportional to the number of genes annotated
with the corresponding GO term.")
        
        
        ##CNET plot
        output$cnet.plot <- renderPlot(
          width     = 870,
          height    = 600,
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
    if( (input$analysis == "kegg"  || input$analysis == "both") && (length(wrong.genes) == 0))
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
        cre.ids <- cre.chlredraft.map$GID
        chlredraft.ids <- cre.chlredraft.map$CHLREDRAFT
        names(chlredraft.ids) <- cre.ids
        names(cre.ids) <- chlredraft.ids
        
        target.genes <- chlredraft.ids[target.genes]
        names(target.genes) <- NULL
        
        gene.universe <- chlredraft.ids[gene.universe]
        names(gene.universe) <- NULL
        
        organism.id <- "cre"
      } else if(input$microalgae == "vcarteri")
      {
        vocar.volcadraft.map <- select(org.Vcarteri.eg.db,columns = c("VOLCADRAFT"),keys=keys(org.Vcarteri.eg.db,keytype = "GID"))
        vocar.ids <- vocar.volcadraft.map$GID
        volcadraft.ids <- vocar.volcadraft.map$VOLCADRAFT
        names(volcadraft.ids) <- vocar.ids
        names(vocar.ids) <- volcadraft.ids
        
        target.genes <- volcadraft.ids[target.genes]
        names(target.genes) <- NULL
        
        gene.universe <- volcadraft.ids[gene.universe]
        names(gene.universe) <- NULL
        
        organism.id <- "vcn"
      } else if(input$microalgae == "ptricornutum")
      {
        phatri.draft.map <- select(org.Ptricornutum.eg.db,columns = c("PHATRIDRAFT"),keys=keys(org.Ptricornutum.eg.db,keytype = "GID"))
        phatri.ids <- phatri.draft.map$GID
        phatridraft.ids <- phatri.draft.map$PHATRIDRAFT
        names(phatridraft.ids) <- phatri.ids
        names(phatri.ids) <- phatridraft.ids
        
        target.genes <- phatridraft.ids[target.genes]
        names(target.genes) <- NULL
        
        gene.universe <- phatridraft.ids[gene.universe]
        names(gene.universe) <- NULL
        
        organism.id <- "pti"
      } else if(input$microalgae == "ngaditana")
      {
        naga.draft.map <- select(org.Ngaditana.eg.db,columns = c("NAGADRAFT"),keys=keys(org.Ngaditana.eg.db,keytype = "GID"))
        naga.ids <- naga.draft.map$GID
        nagadraft.ids <- naga.draft.map$NAGADRAFT
        names(nagadraft.ids) <- naga.ids
        names(naga.ids) <- nagadraft.ids
        
        target.genes <- nagadraft.ids[target.genes]
        names(target.genes) <- NULL
        
        gene.universe <- nagadraft.ids[gene.universe]
        names(gene.universe) <- NULL
        
        organism.id <- "ngd"
      } else if(input$microalgae == "knitens")
      {
        knitens.ko <- select(org.Knitens.eg.db,columns = c("KO"),keys=keys(org.Knitens.eg.db,keytype = "GID"))
        ko.universe <- knitens.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(knitens.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = input$pvalue))
        
        for(i in 1:nrow(pathway.enrichment))
        {
          current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
          
          current.genes <- c()
          for(j in 1:length(current.Ks))
          {
            current.genes <- c(current.genes,subset(knitens.ko, KO == current.Ks[j])$GID)
          }
          
          pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
        }
      }
      
      ## Compute KEGG pathway enrichment
      if (input$microalgae != "knitens")
      {
        pathway.enrichment <- enrichKEGG(gene = target.genes, organism = organism.id, keyType = "kegg",
                                         universe = gene.universe,qvalueCutoff = input$pvalue)
      }
      
      pathway.enrichment.result <- as.data.frame(pathway.enrichment)

      if(nrow(pathway.enrichment.result) > 0)
      {
        kegg.intro.text <- paste(c("This tab presents the results from the <b>KEGG pathways/modules enrichment analysis</b> 
                                      performed over the input target genes (<i>",
                                 paste(gsub(pattern="OT_",replacement="",x=target.genes[1:3]),collapse=" "),
                                 "</i> ...) from the microalgae <b> <i>", microalgae.names[input$microalgae],
                                 " </i> </b> with ", universe.text),collapse="") 
        output$intro_kegg <- renderText(expr = kegg.intro.text)
        
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
        } else if (input$microalgae == "vcarteri")
        {
          kegg.enriched.genes <- pathway.enrichment.result$geneID
          for(i in 1:length(kegg.enriched.genes))
          {
            kegg.enriched.genes[i] <- paste(vocar.ids[strsplit(kegg.enriched.genes[i],split="/")[[1]]],collapse=" ")
          }
        } else if (input$microalgae == "ptricornutum")
        {
          kegg.enriched.genes <- pathway.enrichment.result$geneID
          for(i in 1:length(kegg.enriched.genes))
          {
            kegg.enriched.genes[i] <- paste(phatri.ids[strsplit(kegg.enriched.genes[i],split="/")[[1]]],collapse=" ")
          }
        } else if (input$microalgae == "ngaditana")
        {
          kegg.enriched.genes <- pathway.enrichment.result$geneID
          for(i in 1:length(kegg.enriched.genes))
          {
            kegg.enriched.genes[i] <- paste(naga.ids[strsplit(kegg.enriched.genes[i],split="/")[[1]]],collapse=" ")
          }
        } else if (input$microalgae == "knitens")
        {
          kegg.enriched.genes <- pathway.enrichment.result$geneID
          for(i in 1:length(kegg.enriched.genes))
          {
            kegg.enriched.genes[i] <- paste(strsplit(kegg.enriched.genes[i],split="/")[[1]],collapse=" ")
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
        } else if(input$microalgae == "creinhardtii" || input$microalgae == "vcarteri")
        {
          gene.link.function <- phytozome.gene.link
        } else if(input$microalgae == "ptricornutum")
        {
          gene.link.function <- phaeodactylum.gene.link
        } else if(input$microalgae == "ngaditana")
        {
          gene.link.function <- ngaditana.gene.link
        } else if(input$microalgae == "knitens")
        {
          gene.link.function <- knitens.gene.link
        }
        
        for(i in 1:length(kegg.enriched.genes))
        {
          kegg.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(kegg.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }
        
        ## Add links to kegg pathways
        kegg.result.table.with.links[["KEGG ID"]] <- sapply(X=kegg.result.table.with.links[["KEGG ID"]],FUN = kegg.pathway.link)

        
        ## Introductory text for GO enrichment table
        kegg.table.text <- "The table below summarizes the result of the KEGG pathway
enrichment analysis. Each row represents a pathway significantly enriched in the target
gene set with respect to the selected gene universe. The first column represents the KEGG pathway 
identifier. The second column contains a human readable description. For more details on the 
corresponding pathway, click on the identifier in the first column. The third and fourth 
column presents the p-value and q-value (adjusted p-value or FDR) capturing the level
of significance. The fifth column displays the corresponding enrichment value E (m/n; M/N) where
n and N are the number of genes associated to any pathway in the target set and in the gene universe
respectively; m and M are the number of genes associated to the pathway represented in the corresponding
row in the target gene set and in the gene universe respectively. The enrichment is then computed as
E = (m/n) / (M/N). Finally, the last column, contains the genes from the target gene set
assocated to the enriched pathway represented in the corresponding row." 
        
        output$textKEGGTable <- renderText(expr = kegg.table.text)

        output$output_pathway_table <- renderDataTable({
          kegg.result.table.with.links
        },escape=FALSE,options =list(pageLength = 5)) 
      } else
      {
        print("No KEGG pathway enrichment detected in the input set")
      }
      
      ## Figures for KEGG pathway enrichment analysis
      
      ## Prepare gene set for representation
      if(input$microalgae == "knitens")
      {
        genes.pathway <- rep(0, length(ko.universe))
        names(genes.pathway) <- ko.universe
        
        genes.pathway[target.ko] <- 1
      } else 
      {
        genes.pathway <- rep(0, length(gene.universe))
        names(genes.pathway) <- gene.universe
        
        genes.pathway[target.genes] <- 1
      }

      pathways.for.select <- paste(pathways.result.table[["KEGG ID"]], pathways.result.table[["Description"]], sep=" - ")
      
      kegg.image.text <- "<b> The enriched pathways detected above can be visualized using the dropdown menu below. 
      Genes in the target set associated to the corresponding pathway will be highlighted as red rectangles: </b>"

      output$textKEGGImage <- renderText(expr = kegg.image.text)
      
      output$kegg_selectize <- renderUI({
          selectInput(inputId = "kegg_pathway", label="Choose Pathway for Representation",multiple = FALSE,selected = pathways.for.select[1],
                      choices=pathways.for.select)
      })
    
      if( input$microalgae == "knitens")
      {
        modules.enrichment <- enrichMKEGG(gene = target.ko, universe = ko.universe, organism = "ko", keyType = "kegg",minGSSize = 4)
      } else
      {
        modules.enrichment <- enrichMKEGG(gene = target.genes, universe = gene.universe, organism = organism.id, keyType = "kegg",minGSSize = 4)
      }

      modules.enrichment.result <- as.data.frame(modules.enrichment)

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
        } else if(input$microalgae == "vcarteri")
        {
          modules.enriched.genes <- modules.enrichment.result$geneID
          for(i in 1:length(modules.enriched.genes))
          {
            modules.enriched.genes[i] <- paste(vocar.ids[strsplit(modules.enriched.genes[i],split="/")[[1]]],collapse=" ")
          }
        } else if(input$microalgae == "ptricornutum")
        {
          modules.enriched.genes <- modules.enrichment.result$geneID
          for(i in 1:length(modules.enriched.genes))
          {
            modules.enriched.genes[i] <- paste(phatri.ids[strsplit(modules.enriched.genes[i],split="/")[[1]]],collapse=" ")
          }
        } else if(input$microalgae == "ngaditana")
        {
          modules.enriched.genes <- modules.enrichment.result$geneID
          for(i in 1:length(modules.enriched.genes))
          {
            modules.enriched.genes[i] <- paste(naga.ids[strsplit(modules.enriched.genes[i],split="/")[[1]]],collapse=" ")
          }
        } else if(input$microalgae == "knitens")
        {
          modules.enriched.genes <- modules.enrichment.result$geneID
          for(i in 1:length(modules.enriched.genes))
          {
            modules.enriched.genes[i] <- paste(strsplit(modules.enriched.genes[i],split="/")[[1]],collapse=" ")
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
        } else if(input$microalgae == "creinhardtii" || input$microalgae == "vcarteri")
        {
          gene.link.function <- phytozome.gene.link
        } else if(input$microalgae == "ptricornutum")
        {
          gene.link.function <- phaeodactylum.gene.link
        } else if(input$microalgae == "ngaditana")
        {
          gene.link.function <- ngaditana.gene.link
        } else if(input$microalgae == "knitens")
        {
          gene.link.function <- knitens.gene.link
        }
        
        for(i in 1:length(modules.enriched.genes))
        {
          modules.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(modules.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }
        
        ## Add links to kegg pathways
        modules.result.table.with.links[["KEGG ID"]] <- sapply(X=modules.result.table.with.links[["KEGG ID"]],FUN = kegg.module.link)
        
        
        ## Introductory text for GO enrichment table
        module.table.text <- "The table below summarizes the result of the KEGG module
        enrichment analysis. Modules are sometimes easier to interpret than pathways. Each row represents a module significantly enriched in the target
        gene set with respect to the selected gene universe. The first column represents the KEGG module 
        identifier. The second column contains a human readable description. For more details on the 
        corresponding module, click on the identifier in the first column. The third and fourth 
        column presents the p-value and q-value (adjusted p-value or FDR) capturing the level
        of significance. The fifth column displays the corresponding enrichment value E (m/n; M/N) where
        n and N are the number of genes associated to any module in the target set and in the gene universe
        respectively; m and M are the number of genes associated to the module represented in the corresponding
        row in the target gene set and in the gene universe respectively. The enrichment is then computed as
        E = (m/n) / (M/N). Finally, the last column, contains the genes from the target gene set
        associated to the enriched module represented in the corresponding row." 
        
        output$text_module_kegg <- renderText(expr = module.table.text)

        output$output_module_table <- renderDataTable({
          modules.result.table.with.links
        },escape=FALSE,options =list(pageLength = 5)) 
      }
      
    }

    enriched.pathway.id <- reactive({ 
        if(is.null(input$kegg_pathway))
        {
          return()
        } else
        {
          return(strsplit(input$kegg_pathway,split=" - ")[[1]][1] )
        }
      })
      
    observeEvent(enriched.pathway.id(), {
      
      if(input$microalgae == "creinhardtii")
      {
        organism.id <- "cre"
      } else if(input$microalgae == "otauri")
      {
        organism.id <- "ota"
      } else if(input$microalgae == "vcarteri")
      {
        organism.id <- "vcn"
      } else if(input$microalgae == "ptricornutum")
      {
        organism.id <- "pti"
      } else if(input$microalgae == "ngaditana")
      {
        organism.id <- "ngd"
      } else if(input$microalgae == "knitens")
      {
        organism.id <- "ko"
      }
      
        output$kegg_image <- renderImage({
          pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
                   pathway.id = enriched.pathway.id(),
                   species = organism.id,
                   limit = list(gene=max(abs(genes.pathway)), cpd=1),
                   gene.idtype ="kegg")
          
          list(src = paste(c(enriched.pathway.id(),"pathview","png"), collapse="."),
               contentType="image/png",width=1200,height=900)
        },deleteFile = T)
      })

  })
  
  ## Actions to perform after click on the genomic button
  observeEvent(input$genomic_button,{

    ## Select txdb 
    if(input$microalgae == "otauri") 
    {
      txdb <- TxDb.Otauri.JGI
    } else if (input$microalgae == "creinhardtii")
    {
      txdb <- TxDb.Creinhardtii.Phytozome
    } else if (input$microalgae == "dsalina")
    {
      ## TODO
    } else if (input$microalgae == "vcarteri")
    {
      ## TODO
    } else if (input$microalgae == "ptricornutum")
    {
      txdb <- TxDb.Ptricornutum.Ensembl.Protists
    } else if (input$microalgae == "ngaditana")
    {
      ## TODO
    } else if (input$microalgae == "knitens")
    {
      ## TODO
    } else if (input$microalgae == "bathy")
    {
      ## TODO
    }
    
    ## Extract genomic regions from text box or uploaded file
    if(is.null(input$genomic_regions_file))
    {
      genomic.regions <- as.vector(unlist(strsplit(input$genomic_regions, split="\n",
                                                fixed = TRUE)[1]))
      
      chrs <- vector(mode = "character", length=length(genomic.regions))
      start.points <- vector(mode = "character", length=length(genomic.regions))
      end.points <- vector(mode = "character", length=length(genomic.regions))
      
      for(i in 1:length(genomic.regions))
      {
        current.splitted.row <- strsplit(genomic.regions[i],split="\\s+")[[1]]
        current.splitted.row <- current.splitted.row[current.splitted.row != ""]
        chrs[i] <- current.splitted.row[1]
        start.points[i] <- current.splitted.row[2]
        end.points[i] <- current.splitted.row[3]
      }
      
      genomic.df <- data.frame(chr=chrs,start=start.points,end=end.points)
      genomic.df <- genomic.df[complete.cases(genomic.df),]
      print(head(genomic.df))
      genomic.regions <- makeGRangesFromDataFrame(df = genomic.df, 
                                                  seqnames.field = "chr",
                                                  start.field = "start",
                                                  end.field = "end")
    } else
    {
      genomic.regions <- readPeakFile(peakfile = input$genomic_regions_file$datapath,header=FALSE)
    }
    
    ## Define promoter region around TSS
    promoter <- getPromoters(TxDb=txdb, 
                             upstream=input$promoter_length, 
                             downstream=input$promoter_length)
    
    ## Annotate genomic loci
    peakAnno <- annotatePeak(peak = genomic.regions, tssRegion=c(-input$promoter_length, input$promoter_length),
                             TxDb=txdb)
    
    ## Plot pie chart with annotation
    output$annotation.pie.chart <- renderPlot(width = 940, height = 900, res = 120, {
        plotAnnoPie(peakAnno)
      })
    
    ## Plot distance to tss
    output$distance.to.tss <- renderPlot(width = 940, height = 450, res = 120, {
      plotDistToTSS(peakAnno,
                    title="Distribution of genomic loci relative to TSS",
                    ylab = "Genomic Loci (%) (5' -> 3')")
      
    })
    
    ## Extract genes 
    peak.annotation <- as.data.frame(peakAnno)
    simple.annotation <- sapply(X = as.vector(peak.annotation$annotation),FUN = extract.annotation)
    names(simple.annotation) <- NULL
    
    genes.promoter <- peak.annotation$geneId[simple.annotation == "Promoter"]
    genes.5utr <- peak.annotation$geneId[simple.annotation == "5' UTR"]
    genes.3utr <- peak.annotation$geneId[simple.annotation == "3' UTR"]
    genes.exon <- peak.annotation$geneId[simple.annotation == "Exon"]
    genes.intron <- peak.annotation$geneId[simple.annotation == "Intron"]
    
    ## Select final gene set
    genes <- c()
    if( "Promoter" %in% input$selected_genomic_features )
    {
      genes <- c(genes,genes.promoter)
    }
    if( "5' UTR" %in% input$selected_genomic_features )
    {
      genes <- c(genes,genes.5utr)
    }
    if( "3' UTR" %in% input$selected_genomic_features )
    {
      genes <- c(genes,genes.3utr)
    }
    if( "Exon" %in% input$selected_genomic_features )
    {
      genes <- c(genes,genes.exon)
    }
    if( "Intron" %in% input$selected_genomic_features )
    {
      genes <- c(genes,genes.intron)
    }

    genes <- unique(genes)
    
    print(length(genes))
    print(genes[1:6])
    
    output$annotated_genes <- renderUI({
      selectInput(inputId = "selected_annotated_gene", 
                  label="Choose a gene to inspect its mark profile",
                  multiple = FALSE,selected = genes[1],
                  choices=genes)
    })
    
    selected.annotated.gene.id <- reactive({ 
      if(is.null(input$selected_annotated_gene))
      {
        return()
      } else
      {
        return(input$selected_annotated_gene)
      }
    })
    
    observeEvent(selected.annotated.gene.id(), {
      print("VAAAAAMOOOOOOS!!!!")
      
      
      # if(input$microalgae == "creinhardtii")
      # {
      #   organism.id <- "cre"
      # } else if(input$microalgae == "otauri")
      # {
      #   organism.id <- "ota"
      # } else if(input$microalgae == "vcarteri")
      # {
      #   organism.id <- "vcn"
      # } else if(input$microalgae == "ptricornutum")
      # {
      #   organism.id <- "pti"
      # } else if(input$microalgae == "ngaditana")
      # {
      #   organism.id <- "ngd"
      # } else if(input$microalgae == "knitens")
      # {
      #   organism.id <- "ko"
      # }
      # 
      # output$kegg_image <- renderImage({
      #   pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
      #            pathway.id = enriched.pathway.id(),
      #            species = organism.id,
      #            limit = list(gene=max(abs(genes.pathway)), cpd=1),
      #            gene.idtype ="kegg")
      #   
      #   list(src = paste(c(enriched.pathway.id(),"pathview","png"), collapse="."),
      #        contentType="image/png",width=1200,height=900)
      # },deleteFile = T)
    })
    

  })
})

# Run the application 
shinyApp(ui = ui, server = server)