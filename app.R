## Authors: Ana B. Romero-Losada
##          Pedro de los Reyes 
##          Francisco J. Romero-Campero <fran@us.es>

## Contact & Maintainer: Francisco J. Romero-Campero <fran@us.es>

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
# input <- list(microalgae = "ptricornutum", promoter_length = 1000, genomic_regions_file = "example_files/example_genomic_regions_ptricornutum.txt")
# input <- list(microalgae = "ptricornutum", genomic_regions_file = "example_files/example_genomic_regions_ptricornutum.txt", bw_file= "example_files/example_ptricornutum.bw" ,promoter_length = 1000, selected_genomic_features = "Promoter")
# input <- list(microalgae = "creinhardtii", genomic_regions_file = "example_files/example_genomic_regions_creinhardtii_2.txt", bw_file= "example_files/example_creinhardtii.bw" ,promoter_length = 1000, selected_genomic_features = "Promoter")

# target.genes <- read.table(file="example_files/example_otauri.txt",as.is=T)[[1]]
# target.genes <- read.table(file="cre/examples/activated_genes.txt",as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_vcarteri.txt",as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_ptricornutum.txt", as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_ngaditana_1.txt", as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_knitens.txt", as.is=T)[[1]]
# target.genes <- read.table(file="example_files/example_csubellipsoidea.txt", as.is=T)[[1]]
# target.genes <- read.table(file="/home/fran/tmp/ld_sd/results/activated_genes.txt",as.is=T)[[1]]
# target.genes <- read.table(file="/home/fran/tmp/ld_sd/results/repressed_genes.txt",as.is=T)[[1]]
# input <- list(microalgae = "creinhardtii", pvalue = 0.05, analysis = "go", ontology = "BP", input_mode = "No")

## Increase max file size allowed to upload to 100MB
options(shiny.maxRequestSize=100*1024^2)

## Load necessary packages
library(shinycssloaders)
library(shiny)
library(clusterProfiler)
library(pathview)
library(ChIPseeker)
library(ChIPpeakAnno)
library(rtracklayer)
library(seqinr)
library(shinythemes)
library(shinyjs)
## Load microalgae annotation packages
library(org.Otauri.eg.db) ##install.packages(pkgs = "./packages/annotation_packages/org.Otauri.eg.db/",repos = NULL,type="source")
library(org.Creinhardtii.eg.db)
library(org.Dsalina.eg.db)
library(org.Vcarteri.eg.db)
library(org.Ptricornutum.eg.db)
library(org.Ngaditana.eg.db)
library(org.Knitens.eg.db)
library(org.Csubellipsoidea.eg.db)
library(org.Hlacustris.eg.db)
library(org.Czofingiensis.eg.db)
library(org.Bprasinos.eg.db)

## TODO lucimarinus

## Load microalgae genome annotation packages
library(TxDb.Otauri.JGI)
library(TxDb.Creinhardtii.Phytozome)
library(TxDb.Ptricornutum.Ensembl.Protists)
library(TxDb.Hlacustris.NCBI)
library(TxDb.Czofingiensis.Phytozome)
library(TxDb.Bprasinos.Orcae)
library(TxDb.Csubellipsoidea.Phytozome)
library(TxDb.Dsalina.Phytozome)
library(TxDb.Vcarteri.Phytozome)

microalgae.names <- c("Ostreococcus tauri", 
                      "Chlamydomonas reinhardtii", 
                      "Dunaliella salina",
                      "Volvox carteri",
                      "Phaeodactylum tricornutum",
                      "Nannochloropsis gaditana",
                      "Klebsormidium nitens",
                      "Coccomyxa subellipsoidea",
                      "Bathycoccus prasinos",
                      "Haematococcus lacustris",
                      "Chromochloris zofingiensis",
                      "Ostreococcus lucimarinus")
names(microalgae.names) <- c("otauri", 
                             "creinhardtii", 
                             "dsalina", 
                             "vcarteri",
                             "ptricornutum", 
                             "ngaditana",
                             "knitens",
                             "csubellipsoidea",
                             "bprasinos",
                             "hlacustris",
                             "zofi",
                             "olucimarinus")

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

## Auxiliary function to split a string using commas
split.commas <- function(annotation.str)
{
  return(strsplit(annotation.str,split=",")[[1]])
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

#NCBI link
ncbi.gene.link <- function(gene.name)
{
  ncbi.link <- paste0("https://www.ncbi.nlm.nih.gov/protein/",gene.name)
  gene.link <- paste(c("<a href=\"",
                       ncbi.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

#C zofingiensis gene link
zofi.gene.link <- function(gene.name)
{
  zofi.link <- paste(c("https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=5211&searchText=",
                            gene.name,
                            "&offset=0"),collapse="")
  gene.link <- paste(c("<a href=\"",
                       zofi.link,
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

## KO link
#https://www.genome.jp/dbget-bin/www_bget?ko:K00276
ko.link <- function(ko.term)
{
  link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ko:", ko.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           ko.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## KOG link
## https://www.ncbi.nlm.nih.gov/Structure/cdd/KOG3720
kog.link <- function(kog.term)
{
  link <- paste0("https://www.ncbi.nlm.nih.gov/Structure/cdd/KOG3720", kog.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           kog.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## ENZYME link
## https://www.brenda-enzymes.org/enzyme.php?ecno=1.4.3.21
enzyme.link <- function(ec.term)
{
  link <- paste0("https://www.brenda-enzymes.org/enzyme.php?ecno=", ec.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           ec.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## PANTHER link
## http://www.pantherdb.org/panther/family.do?clsAccession=PTHR10638
panther.link <- function(panther.term)
{
  link <- paste0("http://www.pantherdb.org/panther/family.do?clsAccession=", panther.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           panther.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## PFAM link
## http://pfam.xfam.org/family/PF02728
pfam.link <- function(pfam.term)
{
  link <- paste0("http://pfam.xfam.org/family/", pfam.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           pfam.term, "</a>"),
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

## Function to compute the reverse complement of a DNA sequence
reverse.complement <- function(dna.sequence)
{
  return(c2s(comp(rev(s2c(dna.sequence)),forceToLower = FALSE)))
}

## Load Position Weight Matrices needed in DNA binding motifs search
## Open file connection
con <- file("jaspar_motifs/pfm_plants_20180911.txt",open = "r")

## Empty list for storing PWM
motifs.pwm <- vector(mode="list",length = 453)
motif.ids <- vector(mode="character",length=453)
motif.names <- vector(mode="character",length=453)

## Load 64 PWM
for(j in 1:453)
{
  ## First line contains motif id and name
  first.line <- readLines(con,1)
  
  motif.ids[j] <- strsplit(first.line,split=" ")[[1]][1]
  motif.names[j] <- strsplit(first.line,split=" ")[[1]][2]
  
  ## Next four line contians probabilites for each nucleotide
  a.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  c.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  g.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  t.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  
  ## Construct PWM
  motif.pwm <- matrix(nrow = 4,ncol=length(a.row))
  
  motif.pwm[1,] <- a.row
  motif.pwm[2,] <- c.row 
  motif.pwm[3,] <- g.row
  motif.pwm[4,] <- t.row
  
  rownames(motif.pwm) <- c("A","C","G","T")
  
  motifs.pwm[[j]] <- prop.table(motif.pwm,2)
}

## Close file connection
close(con)

## Naming list with PWM
names(motifs.pwm) <- motif.names
names(motif.ids) <- motif.names

# Define UI
ui <- shinyUI(fluidPage(#theme= "bootstrap.css",
  theme = shinytheme("flatly"),
  
  fluidRow(
    column(
      width = 2,
      img(src='logo_1.png', align = "center", width=200),
      tags$br(),
      radioButtons(inputId = "navigation_bar", width="100%",selected="home",
                   label="",
                   choices=c(
                     "Home" = "home",
                     "MARACAS, MicroAlgae RnA-seq and Chip-seq AnalysiS" = "maracas",
                     "Gene Set Functional Analysis" = "genes",
                     "Genomic Loci Functional Analysis" = "chip",
                     "Tutorials" = "tutorials",
                     "GitHub repository" = "github",
                     "Citation and Contact" = "citation"
                   ))),
    column(
      width = 8,
      tags$div(align = "center", 
               tags$h1(tags$b("ALGAEFUN,"), "microALGAE FUNctional annotation tool, with ",tags$b("MARACAS,"), "MicroAlgae RnA-seq and Chip-seq AnalysiS")),
      tags$br(),tags$br(),
      conditionalPanel(condition = "input.navigation_bar == 'home'",
                       tags$div(align = "justify", "Welcome to", tags$b("ALGAEFUN")," with ", tags$b("MARACAS"), "a microalgae web based tool for the analysis of ", 
                                tags$b("RNA-seq"), "and ", tags$b("ChIP-seq"), "data and the", tags$b("functional annotation"), "of the resulting gene sets and genomic loci. ",
                                tags$b("ALGAEFUN"), "with ", tags$b("MARACAS"), "supports the analysis for a wide collection 
               of microalgae that includes", tags$i("Chlamydomonas reinhardtii, Ostreococcus tauri, Phaeodactylum tricornutum"), "and ", 
                                tags$i("Nannochlorpsis gaditana."), "Please select from the navigation bar on the left the type of analysis you want to perform. You can also 
               see our", tags$b("video tutorial"), "on how to analyse RNA-seq and
               ChIP-seq data as well as on how to functionally annotate gene sets and genomic loci. Our
               code is freely available at", tags$b("Github."), "Please cite our work if you find it useful in your research.")
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'genes'",
                        tags$div(align="justify", tags$b("AlgaeFUN"), "allows researchers to perform", tags$b("functional annotation"), 
                       "over gene sets.", tags$b("Gene Ontology (GO) enrichment"), "analysis as well as", tags$b("KEGG (Kyoto Encyclopedia
                       of Genes and Genomes) pathway enrichment"), "analysis are supported. The gene set of interest can be obtained, for example,
                       as the result of a differential expression analysis carried out using", tags$b("MARACAS."), " See our", tags$b("video tutorial"),
                       "for details or follow the next steps to perform your analysis:",
                                tags$ol(
                                  tags$li("In the left panel choose your ", tags$b("microalgae")," of interest, the type of enrichment analysis 
                                          to perform and the", tags$b("p-value threshold.")),
                                  tags$li("Insert your ", tags$b("gene set"), " in the text box or load it from a file using the",
                                          tags$b("Browse …"), " button. An example can be loaded by clicking on the ",  tags$b("Example"), " button. 
                                          Click on", tags$b("Clear"), " button to remove the loaded gene set."),
                                  tags$li("Users can choose between the default", tags$b("background"), " gene provided by AlgaeFUN of a custom one 
                                          that can be specified."),
                                  tags$li("Click on the ", tags$b("Have Fun"), " button to perform the specified functional enrichment analysis. The
                                          results will be shown in the different tabs below.")
                                 )
      )),
    ),
    column(
      width = 2,
      img(src='logo_ibvf.jpg', align = "center", width=100),
      img(src='logo_us.png', align = "center", width=100),
      tags$br(),tags$br(),tags$br(),
      img(src='logo_csic.jpg', align = "center", width=100)
    )
  ),

  
  tags$br(),tags$br(),
  
  #Interface where the user can choose his/her preferencies, separated by columns
  fluidRow(
      column(width = 4,

      
      conditionalPanel(condition = "input.navigation_bar == 'genes'",
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
                            "Klebsormidium nitens" = "knitens",
                            "Haematococcus lacustris" = "hlacustris",
                            "Chlomochloris zofingiensis" = "zofi")),

        #Choose the kind of analysis that you want us to execute 
        radioButtons(inputId = "analysis",
                   label="Choose your desirable analysis",
                   choices=c("GO terms enrichment" = "go",
                             "KEGG pathways enrichment analysis" = "kegg",
                             "Both" = "both"
                            ))),
      
      conditionalPanel(condition= "(input.analysis == 'go' || input.analysis == 'both') &&
                                   input.navigation_bar == 'genes'",
                       radioButtons(inputId = "ontology",
                                    label="Choose gene ontology:",
                                    choices = c("Biological process" = "BP",
                                                "Cellular Component" = "CC",
                                                "Mollecular Function" = "MF"))),
      #Choose a p-value
      conditionalPanel(condition = "input.navigation_bar == 'genes' || input.navigation_bar == 'chip'",
                       
        numericInput(inputId = "pvalue", 
                     label= "Which will be your chosen p-value?", 
                     value= 0.05)
      ),
      
       
      #Choose a promoter length
      conditionalPanel(condition = "input.navigation_bar == 'chip'",
                       sliderInput(inputId = "promoter_length", 
                                    label= "Choose the distance in base pairs around the 
                                    Transcriptional Start Site to determine gene promoters", 
                                    min=100, max=2000,value=1000,step=100)),
      
      #Choose a promoter length
      conditionalPanel(condition = "input.navigation_bar == 'chip'",
                       checkboxGroupInput(inputId = "selected_genomic_features", 
                                   label= "A gene will be considered as a target of an input genomic locus
                                   when it overlaps one of the following gene parts:",
                                   choices = c("Promoter","5' UTR","Exon","Intron","3' UTR")))
      
      ),
      
      column( width = 8,
        #Choose the kind of analysis that you want us to execute 
        conditionalPanel(condition = "input.navigation_bar == 'genes'",

          #This panel will only appear if the user chooses to use our background lists. 
          textAreaInput(inputId = "genes", label= "Insert a set of genes", width="200%", 
                        height = "200px",placeholder = "Insert set of genes",
                        value= ""
          ),
          actionButton(inputId = "example_genes",label = "Example"),
          actionButton(inputId = "clear_gene_set",label = "Clear"),
          fileInput(inputId = "gene_set_file",label = "Choose File with Gene Set to Upload"),
        
          #The user can either insert his/her own background list or use ours. 
          radioButtons(inputId = "input_mode", width = "100%",
                       label = "Would you rather use your own background set?", 
                       choices = c("Yes", 
                                   "No"),
                       selected = "No"),
          #This panel will only appear if the user wants to use his/her own background list. 
          conditionalPanel(condition = "input.input_mode == 'Yes'",
                         textAreaInput(inputId = "background", label= "Background set", width="200%", 
                                       height = "100px",placeholder = "Insert background list",
                                       value= ""
                         ),
                         
                         actionButton(inputId = "clear_universe_set",label = "Clear"),
                         fileInput(inputId = "gene_universe_file",
                                   label = "Choose File with Custom Gene Universe to Upload",
                                   width = "100%")
                         
        )),

      conditionalPanel(condition = "input.navigation_bar == 'chip'",
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
        fileInput(inputId = "bw_file",label = "Choose BigWig File to Upload for Profile Representations:", width= "100%"),
        actionButton(inputId = "genomic_button",label = "Have fun!", icon("send") )                       
                       ),
      conditionalPanel(condition = "input.navigation_bar == 'genes'",
                       actionButton(inputId = "go.button",label = "Have fun!", icon("send") )                       
      )
    )
  ),
    
    tags$br(), tags$br(),
    mainPanel(width = 13,

      #Main panel containing the results organized in different tabs: GO map, Go terms data table, and 
      #KEGG pathway maps.

      conditionalPanel(condition = "(input.navigation_bar == 'genes')",
         tabsetPanel(type ="tabs",
            tabPanel(tags$b("GO ENRICHMENT"),
                  shinyjs::useShinyjs(),
                  hidden(div(id='loading.enrichment.go',h3('Please be patient, computing GO enrichment ...'))), 
                  hidden(div(id='ready.enrichment.go',h3('Your GO enrichment is ready!'))), 
                  htmlOutput(outputId = "gene_sanity_go"),
                  htmlOutput(outputId = "wrong_genes_go"),
                  htmlOutput(outputId = "intro_go"),
                  tags$br(),
                  tags$br(),
                  tabsetPanel(type = "tabs",
                     tabPanel(tags$b("GO Enrichment Table"),
                              tags$br(),
                              htmlOutput(outputId = "textGOTable"),
                              tags$br(), tags$br(),
                              dataTableOutput(outputId = "output_go_table"),
                              uiOutput(outputId = "download_ui_for_go_table"),
                              htmlOutput(outputId = "revigo")
                     ),
                     tabPanel(tags$b("GO Map"),
                              tags$br(),
                              htmlOutput(outputId = "go_graph"),
                              tags$br(), tags$br(),
                              div(style= "overflow:scroll; height:500px; text-align: center;", 
                                  plotOutput(outputId = "go.plot", inline = T)),
                              tags$br(),
                              tags$br(), tags$br()),
                     tabPanel(tags$b("GO Barplot"),
                              tags$br(),
                             htmlOutput(outputId = "barplot_text"),
                              tags$br(),
                             div(style= "text-align: center;",
                                 plotOutput(outputId = "bar.plot",inline=TRUE)),
                             tags$br(),
                             tags$br(), tags$br()),
                     tabPanel(tags$b("GO Dotplot"),
                            tags$br(),
                              htmlOutput(outputId = "dotplot_text"),
                              tags$br(),
                              div(style= "text-align: center;",
                                  plotOutput(outputId = "dot.plot",inline=TRUE)),
                              tags$br(),
                              tags$br(), tags$br()),
                    tabPanel(tags$b("GO Emap"),
                             tags$br(),
                             htmlOutput(outputId = "emapplot_text"),
                             tags$br(),
                             div(style= "text-align: center;",
                                 plotOutput(outputId = "emap.plot",inline=TRUE)),
                             tags$br(),
                             tags$br(), tags$br()),
                    tabPanel(tags$b("GO Concept Map"),
                             tags$br(),
                             htmlOutput(outputId = "cnetplot_text"),
                             tags$br(),
                             div(style= "text-align: center;",
                                 plotOutput(outputId = "cnet.plot",inline=TRUE)),
                             tags$br(),
                             tags$br(),tags$br())
                    ) # close tabsetPanel for go result
                  ), # close tabPanel for GO ENRICHMENT
            tabPanel(tags$b("KEGG PATHWAY ENRICHMENT"),
                     shinyjs::useShinyjs(),
                     hidden(div(id='loading.enrichment.kegg',h3('Please be patient, computing KEGG pathway enrichment ...'))), 
                     hidden(div(id='ready.enrichment.kegg',h3('Your KEGG enrichment is ready!'))), 
                     htmlOutput(outputId = "gene_sanity_kegg"),
                     htmlOutput(outputId = "wrong_genes_kegg"),
                     htmlOutput(outputId = "intro_kegg"),
                     tags$br(), tags$br(),
                     tabsetPanel(type = "tabs",
                                 tabPanel(tags$b("KEGG Pathway Enrichment Table"), 
                                          tags$br(),
                                          htmlOutput(outputId = "textKEGGTable"),
                                          htmlOutput(outputId = "kegg_pathway_table_text"),
                                          tags$br(),
                                          dataTableOutput(outputId = "output_pathway_table"),
                                          br(), br(), br()),
                                 tabPanel(tags$b("KEGG Pathway Visualization"),
                                          tags$br(), tags$br(),
                                          htmlOutput(outputId = "textKEGGImage"),
                                          br(), br(),
                                          uiOutput(outputId = "kegg_selectize"),
                                          div(style= "text-align: center;",
                                              imageOutput("kegg_image", inline = T)),
                                          tags$br(), tags$br(),
                                          br(), br()),
                                 tabPanel(tags$b("KEGG Module Enrichment Table"),
                                          tags$br(), tags$br(),
                                          htmlOutput(outputId = "text_module_kegg"),
                                          br(), br(),
                                          dataTableOutput(outputId = "output_module_table"),
                                          tags$br(), tags$br(),
                                          tags$br(), tags$br(),
                                          uiOutput(outputId = "kegg_module_selectize"),
                                          imageOutput("kegg_module_image"),
                                          tags$br(), tags$br())
                     )
                     
            )
            )
      ),

      conditionalPanel(condition = "input.navigation_bar == 'chip'",
                       tabsetPanel(type = "tabs",
                                   tabPanel(tags$b("Annotated Genes Table"),
                       tags$br(), tags$br(),
                       hidden(div(id='loading.chip',h3('Please be patient, computing genomic loci annotation ...'))), 
                       hidden(div(id='ready.chip',h3('Your genomic loci annotation is ready!'))),
                       htmlOutput(outputId = "textTableAnnotatedGenes"),
                       tags$br(), tags$br(),
                       dataTableOutput(outputId = "output_gene_chip_table"),
                       tags$br(), tags$br()),
                                   tabPanel("Pie chart",
                       tags$br(), tags$br(),
                       htmlOutput(outputId = "piechart.text"),
                       tags$br(), tags$br(),
                       plotOutput(outputId = "annotation.pie.chart",inline=TRUE),
                       tags$br(), tags$br()), 
                                  tabPanel("Distance to TSS visualization",
                       tags$br(), tags$br(),
                       plotOutput(outputId = "distance.to.tss",inline=TRUE),
                       tags$br(), tags$br()),
                                  tabPanel("TSS signal visualization",
                       tags$br(), tags$br(),
                       plotOutput(outputId = "tss_signal"),
                       tags$br(), tags$br()),
                                  tabPanel("Mark inspection of annotated genes",
                       tags$br(), tags$br(),
                       uiOutput(outputId = "annotated_genes"),
                       tags$br(), tags$br()))#,
                      # plotOutput(outputId = "individual_gene_profile")
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
      gene.link.function <- ostta.gene.link
    } else if (input$microalgae == "creinhardtii")
    {
      org.db <- org.Creinhardtii.eg.db
      microalgae.genes <- read.table(file = "universe/cre_universe.txt",as.is = T)[[1]]
      gene.link.function <- phytozome.gene.link
    } else if (input$microalgae == "dsalina")
    {
      org.db <- org.Dsalina.eg.db
      microalgae.genes <- read.table(file = "universe/dusal_universe.txt",as.is = T)[[1]]
      gene.link.function <- phytozome.gene.link
    } else if (input$microalgae == "vcarteri")
    {
      org.db <- org.Vcarteri.eg.db
      microalgae.genes <- read.table(file = "universe/vocar_universe.txt",as.is = T)[[1]]
      gene.link.function <- phytozome.gene.link
    } else if (input$microalgae == "ptricornutum")
    {
      org.db <- org.Ptricornutum.eg.db
      microalgae.genes <- read.table(file = "universe/phatri_universe.txt",as.is = T)[[1]]
      gene.link.function <- phaeodactylum.gene.link
    } else if (input$microalgae == "ngaditana")
    {
      org.db <- org.Ngaditana.eg.db
      microalgae.genes <- read.table(file = "universe/naga_universe.txt",as.is = T)[[1]]
      gene.link.function <- ngaditana.gene.link
    } else if (input$microalgae == "knitens")
    {
      org.db <- org.Knitens.eg.db
      microalgae.genes <- read.table(file = "universe/klebsor_universe.txt",as.is = T)[[1]]
      gene.link.function <- knitens.gene.link
    }else if (input$microalgae == "bathy")
    {
      org.db <- org.Bprasinos.eg.db
      microalgae.genes <- read.table(file = "universe/bathy_universe.txt",as.is = T)[[1]]
    }else if (input$microalgae == "csubellipsoidea")
    {
      org.db <- org.Csubellipsoidea.eg.db
      microalgae.genes <- read.table(file = "universe/cocsu_universe.txt",as.is = T)[[1]]
      gene.link.function <- phytozome.gene.link
    }else if (input$microalgae == "hlacustris")
    {
      org.db <- org.Hlacustris.eg.db
      microalgae.genes <- read.table(file = "universe/hlacustris_universe.txt",as.is = T)[[1]]
      gene.link.function <- ncbi.gene.link
    }else if (input$microalgae == "zofi")
    {
      org.db <- org.Czofingiensis.eg.db
      microalgae.genes <- read.table(file = "universe/zofi_universe.txt",as.is = T)[[1]]
      gene.link.function <- zofi.gene.link
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
      go.intro.text <- paste(c("The tabs below present the results from the <b>GO enrichment analysis</b> 
                                      performed over the input target genes (<i>",
                               paste(target.genes[1:3],collapse=" "),
                               "</i> ...) from the microalgae <b> <i>", microalgae.names[input$microalgae],
                               " </i> </b> with ", universe.text),collapse="") 
      
      output$intro_go <- renderText(expr = go.intro.text)

      ## Perform GO enrichment
      ## TODO q-value
      shinyjs::showElement(id = 'loading.enrichment.go')
      shinyjs::hideElement(id = 'ready.enrichment.go')
      
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

        ## Add link to genes
        for(i in 1:length(genes.in.go.enrichment))
        {
          go.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(genes.in.go.enrichment[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }

        ## Add links to GO ids
        go.result.table.with.links[["GO ID"]] <- sapply(X = go.result.table.with.links[["GO ID"]], FUN = go.link)
        
        shinyjs::hideElement(id = 'loading.enrichment.go')
        shinyjs::showElement(id = 'ready.enrichment.go')
        
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
        
        ## Generate UI to download go enrichment table and creating a downlodable table
        output$download_ui_for_go_table<- renderUI(
          tagList(downloadButton(outputId= "downloadGOTable", "Download GO Enrichment Table"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
        )
        
        ## Download result
        output$downloadGOTable<- downloadHandler(
          filename= function() {
              paste("go_enrichment_table_",microalgae.names[input$microalgae] , ".tsv", sep="")
            },
          content= function(file) {
            write.table(x = go.result.table,quote = F,sep = "\t",
                        file=file,row.names=FALSE,col.names=TRUE)
          })
        
        ## Link to REVIGO
        revigo.data <- paste(revigo.data <- apply(go.result.table[,c("GO ID", "q-value")], 1, paste, collapse = " "), collapse="\n")

        url1 <- tags$a("here", href="#", onclick="document.revigoForm.submit();")
        url2 <- tags$form(
          name="revigoForm", action="http://revigo.irb.hr/", method="post", target="_blank",
          tags$textarea(name="inputGoList", rows="1", cols="8", class="revigoText",
                        style="visibility: hidden", revigo.data)
        )

        output$revigo<- renderUI(
          tagList("The enriched GO terms above may be redundant. Visualize these results in REViGO in order to remove redundancy. Click", url1, url2)
        )

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
      shinyjs::showElement(id = 'loading.enrichment.kegg')
      shinyjs::hideElement(id = 'ready.enrichment.kegg')
      
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
      }else if(input$microalgae == "hlacustris")
      {
        hlacustris.ko <- select(org.Hlacustris.eg.db,columns = c("KO"),keys=keys(org.Hlacustris.eg.db,keytype = "GID"))
        ko.universe <- hlacustris.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(hlacustris.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = input$pvalue))
        
        for(i in 1:nrow(pathway.enrichment))
        {
          current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
          
          current.genes <- c()
          for(j in 1:length(current.Ks))
          {
            current.genes <- c(current.genes,subset(hlacustris.ko, KO == current.Ks[j])$GID)
          }
          
          pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
        }}else if(input$microalgae == "zofi")
        {
          zofi.ko <- select(org.Czofingiensis.eg.db,columns = c("KO"),keys=keys(org.Czofingiensis.eg.db,keytype = "GID"))
          ko.universe <- zofi.ko$KO
          ko.universe <- ko.universe[!is.na(ko.universe)]
          
          target.ko <- subset(zofi.ko,GID %in% target.genes)$KO
          target.ko <- target.ko[!is.na(target.ko)]
          
          pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = input$pvalue))
          
          for(i in 1:nrow(pathway.enrichment))
          {
            current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
            
            current.genes <- c()
            for(j in 1:length(current.Ks))
            {
              current.genes <- c(current.genes,subset(zofi.ko, KO == current.Ks[j])$GID)
            }
            
            pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
          }}
        
      
      ## Compute KEGG pathway enrichment
      if (input$microalgae != "hlacustris" | input$microalgae != "knitens" | input$microalgae != "zofi" )
      {
        pathway.enrichment <- enrichKEGG(gene = target.genes, organism = organism.id, keyType = "kegg",
                                         universe = gene.universe,qvalueCutoff = input$pvalue)
      }
      shinyjs::showElement(id = 'ready.enrichment.kegg')
      shinyjs::hideElement(id = 'loading.enrichment.kegg')
      
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
        }else if (input$microalgae == "hlacustris")
        {
          kegg.enriched.genes <- pathway.enrichment.result$geneID
          for(i in 1:length(kegg.enriched.genes))
          {
            kegg.enriched.genes[i] <- paste(strsplit(kegg.enriched.genes[i],split="/")[[1]],collapse=" ")
          }
        }else if (input$microalgae == "zofi")
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
        } else if(input$microalgae == "creinhardtii" || input$microalgae == "vcarteri" || input$microalgae == "cocsu" || input$microalgae == "dsalina")
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
        }else if(input$microalgae == "hlacustris")
        {
          gene.link.function <- ncbi.gene.link
        }else if(input$microalgae == "zofi")
        {
          gene.link.function <- zofi.gene.link
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
      if(input$microalgae == "knitens" | input$microalgae == "hlacustris" | input$microalgae == "zofi")
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
    
      if( input$microalgae == "knitens" | input$microalgae == "hlacustris" | input$microalgae == "zofi")
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
        } else if(input$microalgae == "knitens" | input$microalgae == "hlacustris" | input$microalgae == "zofi")
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
        } else if(input$microalgae == "creinhardtii" || input$microalgae == "vcarteri"  || input$microalgae == "cocsu" || input$microalgae == "dsalina")
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
        }else if(input$microalgae == "hlacustris")
        {
          gene.link.function <- ncbi.gene.link
        }else if(input$microalgae == "zofi")
        {
          gene.link.function <- zofi.gene.link
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
      } else if(input$microalgae == "knitens" | input$microalgae == "hlacustris")
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
      gene.link.function <- ostta.gene.link
      txdb <- TxDb.Otauri.JGI
      org.db <- org.Otauri.eg.db
    } else if (input$microalgae == "creinhardtii")
    {
      gene.link.function <- phytozome.gene.link
      txdb <- TxDb.Creinhardtii.Phytozome
      org.db <- org.Creinhardtii.eg.db
    } else if (input$microalgae == "dsalina")
    {
      gene.link.function <- phytozome.gene.link
      txdb <- TxDb.Dsalina.eg.db
      org.db <- org.Dsalina.eg.db
    } else if (input$microalgae == "vcarteri")
    {
      gene.link.function <- phytozome.gene.link
      txdb <- TxDb.Vcarteri.eg.db
      org.db <- org.Vcarteri.eg.db
    } else if (input$microalgae == "ptricornutum")
    {
      gene.link.function <- phaeodactylum.gene.link
      txdb <- TxDb.Ptricornutum.Ensembl.Protists
      org.db <- org.Ptricornutum.eg.db
    } else if (input$microalgae == "ngaditana")
    {
      gene.link.function <- ngaditana.gene.link
      org.db <- org.Ngaditana.eg.db
      ## TODO
    } else if (input$microalgae == "knitens")
    {
      gene.link.function <- knitens.gene.link
      org.db <- org.Knitens.eg.db
      ## TODO
    } else if (input$microalgae == "bathy")
    {
      gene.link.function <- bathy.gene.link
      org.db <- org.Bprasinos.eg.db
      txdb <- TxDb.Bprasinos.Orcae
    }else if (input$microalgae == "cocsu")
    {
      gene.link.function <- phytozome.gene.link
      org.db <- org.Csubellipsoidea.eg.db
      txdb <- TxDb.Csubellipsoidea.Phytozome
    }else if (input$microalgae == "zofi")
    {
      gene.link.function <- zofi.gene.link
      org.db <- org.Czofingiensis.eg.db
      txdb <- TxDb.Czofingiensis.Phytozome
    }else if (input$microalgae == "hlacustris")
    {
      gene.link.function <- ncbi.gene.link
      org.db <- org.Hlacustris.eg.db
      txdb <- TxDb.Hlacustris.NCBI
    }else if (input$microalgae == "olucimarinus")
    {
      #TODO
    }
    
    shinyjs::showElement(id = 'loading.chip')
    shinyjs::hideElement(id = 'ready.chip')
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
      #genomic.regions <- readPeakFile(peakfile = input$genomic_regions_file$datapath,header=FALSE)
      genomic.regions <- readPeakFile(peakfile = input$genomic_regions_file$datapath,header=FALSE)
        #readPeakFile(peakfile = input$genomic_regions_file$datapath,header=FALSE)
    }
    
    ## Define promoter region around TSS
    promoter <- getPromoters(TxDb=txdb, 
                             upstream=input$promoter_length, 
                             downstream=input$promoter_length)
    
    ## Annotate genomic loci
    peakAnno <- annotatePeak(peak = genomic.regions, tssRegion=c(-input$promoter_length, input$promoter_length),
                             TxDb=txdb)
    
    
    ## Introductory text for GO enrichment table
    output$piechart.text <- renderText(expr = "<b>The following piechart summarrepresents the distribution of the 
    gene parts overlaped by the input genomic loci.</b>")
    
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

    ## Output table with gene annotation
    annotations <- intersect(c("GO", "KO", "KOG", "ENZYME", "PANTHER","PFAM"),columns(org.db))
    microalgae.annotation <- select(org.db,columns = annotations,keys=keys(org.db,keytype = "GID"))
    genes.annotation <- subset(microalgae.annotation,GID %in% genes)
    
    genes.annotation.download <- data.frame(matrix(nrow=length(genes),ncol=(length(annotations)+1)))
    genes.annotation.links <- genes.annotation.download 
      
    colnames(genes.annotation.download)[1] <- "Gene ID"
    colnames(genes.annotation.links)[1] <- "Gene ID"
    
    for(i in 1:length(annotations))
    {
      current.annotation <- annotations[i]
      
      colnames(genes.annotation.download)[i+1] <- current.annotation
      colnames(genes.annotation.links)[i+1] <- current.annotation
      
      if(current.annotation == "GO")
      {
        annotation.link <- go.link
      } else if (current.annotation == "KO")
      {
        annotation.link <- ko.link
      } else if (current.annotation == "KOG")
      {
        annotation.link <- kog.link  
      } else if (current.annotation == "ENZYME")
      {
        annotation.link <- enzyme.link
      } else if (current.annotation == "PANTHER")
      {
        annotation.link <- panther.link
      } else if (current.annotation == "PFAM")
      {
        annotation.link <- pfam.link
      }
      
      j<-1
      for(j in 1:length(genes))
      {
        current.gene <- genes[j]
        genes.annotation.download[j,1] <- current.gene
        genes.annotation.links[j,1] <- gene.link.function(current.gene)
        current.gene.annotation <- sapply(unique(subset(genes.annotation, GID == current.gene)[[current.annotation]]),split.commas)
        if(is.na(current.gene.annotation[1]))
        {
          genes.annotation.download[j,(i+1)] <- ""
          genes.annotation.links[j,(i+1)] <- ""
        } else
        {
          genes.annotation.download[j,(i+1)] <- paste(current.gene.annotation,collapse=" ")
          genes.annotation.links[j,(i+1)] <- paste(sapply(current.gene.annotation,annotation.link),collapse=" ")
        }
      }
    }

    shinyjs::showElement(id = 'ready.chip')
    shinyjs::hideElement(id = 'loading.chip')
    ## Introductory text for GO enrichment table
    annotated.genes.table.text <- "<b>The table below enumerates the potential gene targets associated with the input
    genomic loci. A gene is associated as a target of a genomic locus when it overlaps at least one of the selected
    gene parts (i.e. promoter, exon, etc.). Each row represents a gene with its available annotation. Click on the 
    gene id to access information from the corresponding data base. Click on the annotation ids to access a detailed 
    description. You can select the number of genes displayed in each page of the table with the <i>Show</i> 
    dropdown menu below. You can use the <i>Search</i> box to identify a specific gene or annotation. The search
    boxes at the end of the table can be used to search in a specific column. Finally, this table can be downloaded
    by clicking on the button below.</b>"

    output$textTableAnnotatedGenes <- renderText(expr = annotated.genes.table.text)
    
    ## Output table with annotated genes
    output$output_gene_chip_table <- renderDataTable({
      genes.annotation.links
    },escape=FALSE,options =list(pageLength = 10)) 

    ## Plot signal around tss
    ## Extraction of the genomic features of the specified genes.
    genes.data <- subset(genes(txdb), gene_id %in% genes)
    
    ## Extraction of the TSS 
    genes.tss <- resize(genes.data, width=1, fix='start')
    
    ## Centering around TSS with promoter length
    around.genes.tss <- genes.tss
    start(around.genes.tss) <- start(genes.tss) - input$promoter_length
    end(around.genes.tss) <- end(genes.tss) + input$promoter_length
    
    # print("around TSS")
    # print(around.genes.tss)
    # print(input$bw_file$data)
    # print(input$bw_file)
    ## Importing bigWig file
    cvglists <- sapply(input$bw_file$data, import,
                       format="BigWig",
                       which=around.genes.tss,
                       as="RleList")
      # sapply(input$bw_file$data, import,
      #                  format="BigWig",
      #                  which=around.genes.tss,
      #                  as="RleList")

    ## Extracting the signal around TSS with promoter length
    number.tiles <- 2*input$promoter_length/20
    tss.sig <- featureAlignedSignal(cvglists, around.genes.tss,
                                    upstream=input$promoter_length,
                                    downstream=input$promoter_length,
                                    n.tile=number.tiles)
    
    ## Extraction of the TES
    genes.tes <- resize(genes.data, width=1, fix='end')
    
    ## Centering around TES 
    around.genes.tes <- genes.tes
    start(around.genes.tes) <- start(genes.tes) - input$promoter_length
    end(around.genes.tes) <- end(genes.tes) + input$promoter_length
    
    print("around TES")
    print(around.genes.tes)
    
    
    ## Extracting the signal 2kb from the center with 50 tile (This may be slow)
    number.tiles.tes <- 2 * input$promoter_length /20
    tes.sig <- featureAlignedSignal(cvglists, around.genes.tes, 
                                    upstream=input$promoter_length, 
                                    downstream=input$promoter_length,
                                    n.tile=number.tiles.tes) 
    
    
    print("TSS signal:")
    print(tss.sig)
    output$tss_signal <- renderPlot({
      
    par(mfrow=c(1,2))
    ## Plotting the results around TSS (you may have to change the names of the conditions)
    profile.around.tss <- colMeans(tss.sig[[1]],na.rm = TRUE)
    profile.around.tes <- colMeans(tes.sig[[1]],na.rm = TRUE)
    max.y <- max(c(profile.around.tss, profile.around.tes))
    plot(profile.around.tss,type="l",col="blue",lwd=3,ylab="",
         cex.lab=2,axes=FALSE,xlab="",main="Signal around TSS",cex.main=1.5,ylim=c(0,max.y))
    polygon(c(1,1:length(profile.around.tss),length(profile.around.tss)),
            c(0,profile.around.tss,0),col="lightblue")
  
    axis(side = 1,
         labels = c(-input$promoter_length,-input$promoter_length/2,
                    "TSS",
                    input$promoter_length/2,input$promoter_length),
         at = c(1,number.tiles/4,number.tiles/2,3*number.tiles/4,number.tiles),lwd=2,cex=1.5,las=2,cex=2)
      
    plot(profile.around.tes,type="l",col="red4",lwd=3,ylab="",
         cex.lab=2,axes=FALSE,xlab="",main="Signal around TES",cex.main=1.5,ylim=c(0,max.y))#,ylim=c(0,830))
    polygon(c(1,1:length(profile.around.tes),length(profile.around.tes)),
            c(0,profile.around.tes,0),col="bisque")
      
    axis(side = 1,
         labels = c(-input$promoter_length,-input$promoter_length/2,
                    "TES",
                    input$promoter_length/2,input$promoter_length),
         at = c(1,number.tiles.tes/4,number.tiles.tes/2,3*number.tiles.tes/4,number.tiles.tes),lwd=2,cex=1.5,las=2,cex=2)

    })
    
    ## Extract info from genes for profile representations
    genes.data.df <- as.data.frame(genes.data)
    exons.data <- as.data.frame(exons(txdb))
    cds.data <- as.data.frame(cds(txdb))

    output$annotated_genes <- renderUI({
      fluidRow(
        column (6, 
                selectInput(inputId = "selected_annotated_gene", 
                            label="Choose a gene to inspect its mark profile",
                            multiple = FALSE,selected = genes[1],
                            choices=genes), 
                ## Selectize to choose target gene to represent
                selectizeInput(inputId = "selected.motifs",
                               label = "Select Motifs",
                               choices = motif.names,
                               multiple = TRUE),
                
                ## Checkbox to select all available motifs
                checkboxInput(inputId = "all.motifs",
                              label = "Select All Motifs:",
                              value = FALSE),
                
                ## Numeric input for PWM score
                numericInput(inputId = "min.score.pwm", 
                             label = "Minimum Score for Motif Identification:",
                             value = 100, 
                             min = 80,
                             max = 100,
                             step = 5),
                actionButton(inputId = "individual_gene_mark",label = "Go")
        ),
        column(6,
               plotOutput(outputId = "individual_gene_profile")
        )
      )
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
    
    observeEvent(input$individual_gene_mark, { #selected.annotated.gene.id(), {
#      print("VAAAAAMOOOOOOS!!!!")
      gene.name <- selected.annotated.gene.id()
      
      target.gene.body <- genes.data.df[gene.name,]
      target.gene.chr <- as.character(target.gene.body$seqnames)
      target.gene.start <- target.gene.body$start
      target.gene.end <- target.gene.body$end
      
      target.gene.strand <- as.character(target.gene.body$strand)
      
      ## Extract cds annotation
      cds.data.target.gene <- subset(cds.data, 
                                     seqnames == target.gene.chr & 
                                       (start >= target.gene.start & 
                                          end <= target.gene.end))
      
      ## Extract exons annotation
      exons.data.target.gene <- subset(exons.data, 
                                       seqnames == target.gene.chr & 
                                         (start >= target.gene.start & 
                                            end <= target.gene.end))
      
      ## Determine the genome range to plot including promoter, gene body and 5' UTR
      ## This depends on whether the gene is on the forward or reverse strand
      range.to.plot <- target.gene.body
      
      if(target.gene.strand == "+")
      {
        range.to.plot$start <- range.to.plot$start - input$promoter_length
        range.to.plot$end <- range.to.plot$end + input$promoter_length
      } else if (target.gene.strand == "-")
      {
        range.to.plot$end <- range.to.plot$end + input$promoter_length
        range.to.plot$start <- range.to.plot$start - input$promoter_length
      }
      
      ## Compute the length of the genome range to represent
      current.length <- range.to.plot$end - range.to.plot$start
      
      ## Compute profile in gene
      selected.bigwig.files <- input$bw_file$data
      selected.bed.files <- input$genomic_regions_file$data
      
      ## Since ChIPpeakAnno needs more than one region to plot our region
      ## is duplicated 
      regions.plot <- GRanges(rbind(range.to.plot,range.to.plot))
      
      ## Import signal from the bigwig files
      cvglists <- sapply(selected.bigwig.files, import, 
                         format="BigWig", 
                         which=regions.plot, 
                         as="RleList")
      
      ## Compute signal in the region to plot
      chip.signal <- featureAlignedSignal(cvglists, regions.plot, 
                                          upstream=ceiling(current.length/2), 
                                          downstream=ceiling(current.length/2),
                                          n.tile=current.length) 
      
      ## Compute mean signal 
      if(target.gene.strand == "+")
      {
        chip.signal.mean <- colMeans(chip.signal[[1]],na.rm = TRUE)
      } else if (target.gene.strand == "-")
      {
        chip.signal.mean <- rev(colMeans(chip.signal[[1]],na.rm = TRUE))
      }
      
      ## Normalization
      chip.signal.mean <- 20 * chip.signal.mean / max(chip.signal.mean)
      
      ## Determine upper limit of the graph
      upper.lim <- 21
      
      ## Height to draw DNA strand
      gene.height <- -25
      cord.x <- 1:current.length
      
      ## Exon width to plot
      exon.width <- 2
      
      ## Cds width to plot
      cds.width <- 3
      
      ## Width of the rectangule representing the peak reagion
      peak.width <- 1
      
      ## Extract exons for target gene
      exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      ## Transform exon coordinates to current range
      min.pos <- min(exons.data.target.gene$start)
      
      if(target.gene.strand == "+")
      {
        exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$promoter_length
        exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$promoter_length
      } else if(target.gene.strand == "-")
      {
        exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$tes_length
        exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$tes_length
      }
      
      ## Extract cds for target gene
      cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      ## Transform cds coordinates to current range
      if(target.gene.strand == "+")
      {
        cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$promoter_length
        cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$promoter_length
      } else if (target.gene.strand == "-")
      {
        cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$tes_length
        cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$tes_length
      }
      
      output$individual_gene_profile <- renderPlot(width = 940, height = 700, res = 120, {
        plot(cord.x, rep(gene.height,length(cord.x)),type="l",col="black",lwd=3,ylab="",
             cex.lab=2,axes=FALSE,xlab="",main="",cex.main=2,
             ylim=c(-35,upper.lim),xlim=c(-3000,max(cord.x)))
        
        ## Represent exons
        for(i in 1:nrow(exons.data.target.gene))
        {
          # Determine start/end for each exon
          current.exon.start <- exons.data.target.gene$start[i]
          current.exon.end <- exons.data.target.gene$end[i]
          
          ## Determine coordinates for each exon polygon and represent it
          exon.x <- c(current.exon.start,current.exon.end,current.exon.end,current.exon.start)
          exon.y <- c(gene.height + exon.width, gene.height + exon.width, gene.height - exon.width, gene.height - exon.width)
          
          polygon(x = exon.x, y = exon.y, col = "blue",border = "blue")
        }
        
        for(i in 1:nrow(cds.data.target.gene))
        {
          # Determine current cds start/end
          current.cds.start <- cds.data.target.gene$start[i]
          current.cds.end <- cds.data.target.gene$end[i]
          
          # Determine curret cds coordinates for the polygon and represent it
          cds.x <- c(current.cds.start,current.cds.end,current.cds.end,current.cds.start)
          cds.y <- c(gene.height + cds.width, gene.height + cds.width, gene.height - cds.width, gene.height - cds.width)
          
          polygon(x = cds.x, y = cds.y, col = "blue",border = "blue")
        }
        
        ## Draw arrow to represent transcription direction 
        if(target.gene.strand == "+")
        {
          lines(c(input$promoter_length,input$promoter_length,input$promoter_length+100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
          lines(c(input$promoter_length+50,input$promoter_length+100),y=c(gene.height+6,gene.height+5),lwd=3)
          lines(c(input$promoter_length+50,input$promoter_length+100),y=c(gene.height+4,gene.height+5),lwd=3)
        } else if (target.gene.strand == "-")
        {
          lines(c(current.length - input$promoter_length, current.length - input$promoter_length, current.length - input$promoter_length-100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
          lines(c(current.length - input$promoter_length-50, current.length - input$promoter_length - 100),y=c(gene.height + 6, gene.height + 5),lwd=3)
          lines(c(current.length - input$promoter_length-50, current.length - input$promoter_length - 100),y=c(gene.height + 4, gene.height + 5),lwd=3)
        }
        
        ## Draw promoter range
        if(target.gene.strand == "+")
        {
          axis(side = 1,labels = c(- input$promoter_length, - input$promoter_length / 2,"TSS"),at = c(1,input$promoter_length/2,input$promoter_length),lwd=2,cex=1.5,las=2,cex=2)
        } else if(target.gene.strand == "-")
        {
          axis(side = 1,labels = c("TSS",- input$promoter_length / 2,- input$promoter_length),at = c(current.length-input$promoter_length,current.length-input$promoter_length/2, current.length),lwd=2,cex=1.5,las=2,cex=2)
        }
        
        ## Draw gene name
        text(x = current.length / 2, y = -33 , 
             labels = bquote(italic(.(gene.name))),cex = 1.7,font = 3)
        
        ## Extract bed file name 1 and read it
        current.peaks <- read.table(file=selected.bed.files,header = F, as.is = T)
        peak.coordinates <- subset(current.peaks, V1 == range.to.plot$seqnames & V2 >= range.to.plot$start & V3 <= range.to.plot$end) 
        current.peaks.to.plot <- peak.coordinates[,2:3]
        
        ## Transform coordinates 
        current.peaks.to.plot <- current.peaks.to.plot - range.to.plot$start
        
        ## Check if there are peaks for the target gene
        if(nrow(current.peaks.to.plot) > 0)
        {
          #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
          for(j in 1:nrow(current.peaks.to.plot))
          {
            ## Extract start and end point of each peak region
            current.peak.start <- current.peaks.to.plot[j,1]
            current.peak.end <- current.peaks.to.plot[j,2]
            
            ## Computer coordinates for polygon and draw it
            peak.x <- c(current.peak.start,current.peak.end,
                        current.peak.end,current.peak.start)
            peak.y <- c(peak.width - 12,   peak.width - 12, 
                        - peak.width - 12, - peak.width - 12)  
            
            polygon(x = peak.x, y = peak.y, col = "bisque", border = "darkred",lwd=2)
          }
        }
        
        line.colors <- "blue"
        area.colors <- "lightblue"
        
        ## Draw profiles
        ## Compute base line for current TF
        current.base.line <- - 10
        
        ## Represent signal from the current TF
        lines(chip.signal.mean+current.base.line,type="l",col=line.colors,lwd=3)
        
        ## Determine polygon coordinates and represent it
        cord.y <- c(current.base.line,chip.signal.mean+current.base.line,current.base.line)
        cord.x <- 1:length(cord.y)
        polygon(cord.x,cord.y,col=area.colors)
        
        ## Load reference genome for the chosen microalgae
        microalgae.genome.data <- read.fasta(file = paste(c("genomes/",input$microalgae,".fa"),collapse=""),
                                             seqtype = "DNA")
        microalgae.genome <- getSequence(microalgae.genome.data)
        names(microalgae.genome) <- getName(microalgae.genome.data)

        ## Determine TFBS motifs to search for
        if(input$all.motifs)
        {
          selected.motifs.pwm <- motifs.pwm
        } else
        {
          selected.motifs.pwm <- motifs.pwm[input$selected.motifs]
        }
        
        selected.motif.names <- names(selected.motifs.pwm)
        selected.motif.ids <- motif.ids[selected.motif.names]
        
        ## Initialize data frame containing TF binding sequences in the peak regions
        df.hits <- data.frame(0,"","","")
        colnames(df.hits) <- c("position","id","name","seq")
        
        ## Identify TF binding DNA motifs 
        if(nrow(current.peaks.to.plot) > 0)
        {
          #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
          for(j in 1:nrow(current.peaks.to.plot))
          {
            ## Genomic coordinates of the current peak
            peak.chr <- peak.coordinates[j, 1]
            peak.start <- peak.coordinates[j, 2]
            peak.end <- peak.coordinates[j, 3]
            
            ## Extract start and end point of each peak region in our plot
            current.peak.start <- current.peaks.to.plot[j,1]
            current.peak.end <- current.peaks.to.plot[j,2]
            
            
            ## Extract peak sequence
            peak.sequence <- c2s(microalgae.genome[[peak.chr]][peak.start:peak.end])
            peak.rev.comp.sequence <- reverse.complement(peak.sequence)
            
            for(k in 1:length(selected.motifs.pwm))
            {
              print(k)
              motif.pwm <- selected.motifs.pwm[[k]]
              
              hits.fw <- matchPWM(motif.pwm, peak.sequence, 
                                  min.score = paste0(input$min_score_pwm,"%"))
              hits.fw
              hits.fw.seqs <- as.data.frame(hits.fw)[[1]]
              hits.fw <- as(hits.fw, "IRanges")
              hits.fw.start <- start(hits.fw)
              hits.fw.end <- end(hits.fw)
              
              if(length(hits.fw.start) > 0)
              {
                df.hits.fw <- data.frame(((hits.fw.start+hits.fw.end)/2) + current.peak.start,
                                         rep(selected.motif.ids[k],length(hits.fw.start)),
                                         rep(selected.motif.names[k],length(hits.fw.start)),
                                         hits.fw.seqs)
                colnames(df.hits.fw)  <- c("position","id","name","seq")
                df.hits <- rbind(df.hits,df.hits.fw)
              }
              
              hits.rev <- matchPWM(motif.pwm, peak.rev.comp.sequence, 
                                   min.score = paste0(input$min.score.pwm,"%"))
              hits.rev.seqs <- as.data.frame(hits.rev)[[1]]
              hits.rev.seqs <- sapply(hits.rev.seqs,reverse.complement)
              names(hits.rev.seqs) <- NULL
              
              hits.rev <- as(hits.rev, "IRanges")
              hits.rev.start <- nchar(peak.sequence) - end(hits.rev) + 1
              hits.rev.end <- nchar(peak.sequence) - start(hits.rev) + 1
              
              if(length(hits.rev.start) > 0)
              {
                df.hits.rev <- data.frame(((hits.rev.start+hits.rev.end)/2) + current.peak.start,
                                          rep(selected.motif.ids[k],length(hits.rev.start)),
                                          rep(selected.motif.names[k],length(hits.rev.start)),
                                          hits.rev.seqs)
                colnames(df.hits.rev)  <- c("position","id","name","seq")
                df.hits <- rbind(df.hits,df.hits.rev)
              }
            }
          }
        }
        
        ## Remove first line of the data frame added just for technical reason
        df.hits <- df.hits[-1,]
        nrow(df.hits)
        
        ## Draw TF binding sites
        detected.tfbs <- unique(as.vector(df.hits$name))
        
        ## TF binding sites colors and symbol shapes
        symbol.shapes <- c(17, 18, 19, 15)
        symbol.color <- c("blue", "red", "darkgreen", "magenta")
        
        number.of.shapes <- ceiling(length(detected.tfbs) / length(symbol.color))
        necessary.shapes <- rep(symbol.shapes[1:number.of.shapes],each = length(detected.tfbs)/number.of.shapes)
        necessary.colors <- rep(symbol.color,number.of.shapes)
        
        if(length(detected.tfbs) > 0)
        {
          for(i in 1:length(detected.tfbs))
          {
            current.tfbs <- detected.tfbs[i]
            current.shape <- necessary.shapes[i]
            current.color <- necessary.colors[i]
            
            positions <- subset(df.hits, name == current.tfbs)
            
            for(j in 1:nrow(positions))
            {
              pos.to.draw <- positions$position[j]
              
              points(x = pos.to.draw, y = -15,
                     pch = current.shape, col = current.color, cex = 1)
            }
          }
          
          ## Add legend for TFBS
          legend.step <- 5
          for(i in 1:length(detected.tfbs))
          {
            points(x = -3000, y = upper.lim - (i-1)*legend.step, 
                   pch=necessary.shapes[i], col = necessary.colors[i],cex = 1)
            
            
            current.seq <- as.character(subset(df.hits,name == detected.tfbs[i])[["seq"]][[1]])
            current.label <- paste(c(detected.tfbs[i], "  -  ", current.seq ),collapse="")
            
            text(x = -2900, y = upper.lim - (i-1)*legend.step, labels = current.label,
                 adj = 0,cex = 0.7)
          }
        }
      })
    })
  })
})

# Run the application 
shinyApp(ui = ui, server = server)