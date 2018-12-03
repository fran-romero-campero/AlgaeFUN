#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shinycssloaders)
library(shiny)
library(topGO)

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
        selectInput(inputId = "mapfile", label="Choose your favourite microalgae", 
                  choices=c("Ostreococcus tauri" = "otauri",
                            "Ostreococcus lucimarinus" = "olucimarinus",
                            "Coccomyxa subellipsoidea" = "csubellipsoidea",
                            "Dunaliella salina" = "dsalina",
                            "Chlamydomonas reinhardtii" = "creinhardtii",
                            "Bathycoccus prasinos" = "bathy")),
             #Choose a p-value
      numericInput(inputId = "pvalor", label= "Which will be your chosen p-value?", value= 0.05),
             #Choose the kind of analysis that you want us to execute 
      radioButtons(inputId = "analysis_asked_for",
                   label="Choose your desirable analysis",
                   choices=c("KEGG analysis",
                             "GO terms enrichment",
                             "Both of them"
                            ))      
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
                  tabPanel("GO map", div(style= "overflow:scroll; height:500px;", 
                                       plotOutput(outputId = "plot", inline = TRUE)),
                           downloadButton(outputId= "downloadImage", "Get GO map")),
                  tabPanel("KEGG pathway", plotOutput("keggpath"), 
                           downloadButton(outputId= "downloadKEGGImage", "Get KEGG pathway image")),
                  tabPanel("Summary", dataTableOutput(outputId = "data"),
                           downloadButton(outputId= "downloadData", "Get GO terms of each gene"),
                           htmlOutput(outputId = "revigo"))
                 )
              )
      ######tabpanels also can look like this: 
      #https://shiny.rstudio.com/gallery/navlistpanel-example.html
    ))
##################################Server not modified yet###################################

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  observeEvent(input$go.button , {
    p.valor <- input$pvalor
    
    print(p.valor)
    
    # p.valor<- 0.5
    
    target.genes <- as.vector(unlist(
      strsplit(input$genes, split="\n",
               fixed = TRUE)[1]))
    
    
    #Read annotation file
    
    map.file <- paste(c("data/data_",input$mapfile,".map"), collapse = "")
    
    geneID2GO <- readMappings(file=map.file)
    
    #Set background
    gene.names <- attributes(geneID2GO)[[1]]
    gene.background <- rep(1, length(gene.names))
    names(gene.background) <- gene.names
    
    gene.background[target.genes] <- 0
    
    #Every target gene in background = 0
    gene.selec <- function(gene.list)
    {
      return(gene.list == 0)
    }
    #Make the object topGOdata
    sampleGOdata <- new("topGOdata",
                        description = "Microalgae session", ontology = "BP",
                        allGenes = gene.background, geneSel = gene.selec,
                        nodeSize = 10,
                        annot = annFUN.gene2GO, gene2GO = geneID2GO)
    
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    
    #Create the table where you'll write the most relevant 100 GO terms 
    allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                       ranksOf = "classicFisher", topNodes = 100, numChar=100)
    
    #Save p-value column and create a vector where you'll save the corrected p-values.
    all.pvalues <- as.vector(allRes[["classicFisher"]])
    
    corrected.all.pvalues <- vector(length=length(all.pvalues), mode="numeric")
    #Change value of the smallest p-value to 0, because they are so small that
    #they can't be seen as numeric
    for (i in 1:length(all.pvalues))
    {
      if (is.na(as.numeric(all.pvalues[i])))
      {
        corrected.all.pvalues[i] <- 0
      }else
      {
        corrected.all.pvalues[i] <- as.numeric(all.pvalues[i])
      }
    }
    #Take just the p-values that are smaller than our p-value pre-selected.
    filter.Res <- allRes[corrected.all.pvalues < p.valor,]
    res.with.genes <- data.frame(filter.Res,genes=vector(length=nrow(filter.Res),mode="character"))
    
    # #Loop that will save the GO terms of the target genes that have a significant
    # #p-value.
    #Save the GO terms related to each target gene
    target.genes.annotation <- geneID2GO[target.genes]
    
    #Save 100 GO terms selected
    significative.GO.terms <- as.vector(filter.Res[["GO.ID"]])
    
    #Create the list of all the possible GO terms 
    offspring <- as.list(GOBPOFFSPRING)
    
    #Create the empty vector where you'll save the genes in every GOterm selected.
    genes.for.each.GO <- vector(length=length(significative.GO.terms),mode="character")
    
    
    for (i in 1:length(significative.GO.terms))
    {
      #for each GO terms of the 100 GO terms selected, get all the GO
      #terms related to it
      GO.genes <- vector(mode="character")
      k <- 1
      go.offspring <- offspring[[significative.GO.terms[i]]]
      for (j in 1:length(target.genes.annotation))
      {
        #for every GO term related to each target gene
        if(length(intersect(go.offspring,target.genes.annotation[[j]]))!=0)
        {
          #if there are common GO terms between the GO terms related to each
          #target gene and the GO terms related to each of the 100 GO terms selected,
          #save the gene in GO.genes vector that you created before.
          GO.genes[k] <- names(target.genes.annotation[j])
          k <- k+1
        }
      }
      GO.genes.string <- paste(unique(GO.genes), collapse=" ")
      
      genes.for.each.GO[i] <- GO.genes.string
      print(i)
    }
    
    res.with.genes[["genes"]] <- genes.for.each.GO
    
    #Generate GOmap   
    output$plot<- 
      renderPlot( 
        width     = 4700,
        height    = 4500,
        res       = 260,
        {
          
          ## Create the final image
          showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
          
        })
    
    
    ## Download the final image
    output$downloadImage<- downloadHandler(
      filename= "GO_map.png",
      content= function(file) {
        png(file, width     = 10,
            height    = 10,
            units     = "in",
            res       = 600)
        showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
        dev.off()
      })  
    
    #Function to generate the links that will take you to a web about each gene.
    createGeneLink<- function (gene.name)
    {
      paste(c("<a href=\"https://www.ncbi.nlm.nih.gov/gene/?term=", 
              gene.name,
              "\" target=\"_blank\">",
              gene.name,
              "</a>"
      ), collapse= "")
    }
    #Function to generate the links that will take you to a web about each GO term.
    
    createGOLink<- function (go.term)
    {
      
      paste(c("<a href=\"http://amigo.geneontology.org/amigo/term/", 
              go.term,
              "\" target=\"_blank\">",
              go.term,
              "</a>"
      ), collapse= "")
    }
    #make a copy of res.with.genes before changing the data that it contains 
    res.with.genes.1<- cbind (res.with.genes)
    output$downloadData<- downloadHandler(
      filename= c('GO_terms.txt'),
      content= function(file) {
        write.table(res.with.genes.1, file, row.names = TRUE)
      })
    
    #Use the previous functions before printing the final data table on screen      
    res.with.genes$GO.ID <- sapply(res.with.genes$GO.ID, createGOLink)
    
    for (i in 1:length(res.with.genes$genes))
    {
      vector.genes <-strsplit(res.with.genes$genes[[i]], split=" ")[[1]]
      res.with.genes$genes[[i]] <- paste(sapply(vector.genes, createGeneLink),collapse= " ")
    }
    
    #Print the final data table on screen
    output$data<- renderDataTable({
      res.with.genes
    }, escape= FALSE)
    
    
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
    j<-1
    revigo.data<-c()
    #Generate the input data in the correct form: "GO p.value"
    for (j in 1:length(res.with.genes.1$GO.ID))
    {
      each.element<-paste(res.with.genes.1$GO.ID[j], res.with.genes.1$classicFisher[j], collapse = " ")
      revigo.data[j]<-paste(each.element,collapse = "")
    }  
    revigo.data<-paste(revigo.data, collapse="\n")
    
    
    url1<- a("here", href="#", onclick="document.revigoForm.submit();")
    url2<- tags$form(
      name="revigoForm", action="http://revigo.irb.hr/", method="post", target="_blank", 
      tags$textarea(name="inputGoList", rows="1", cols="8", class="revigoText", 
                    style="visibility: hidden", revigo.data) 
    )
    
    
    output$revigo<- renderUI(
      tagList("Visualize output in REViGO", url1, url2)
    )
    
    
    #                    
    #                            
    
  })
  
  
})

# Run the application 
shinyApp(ui = ui, server = server)