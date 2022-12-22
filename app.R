## Authors: Ana B. Romero-Losada
##          Francisco J. Romero-Campero <fran@us.es>
##          Christina Arvanitidou
##          Marcos Ramos-Gonzalez
##          Pedro de los Reyes 


## Contact & Maintainer: Francisco J. Romero-Campero <fran@us.es>

## Increase max file size allowed to upload to 100MB
options(shiny.maxRequestSize=100*1024^2)

## Load necessary packages
library(shinycssloaders)
library(shiny)
#library(clusterProfiler)
#library(pathview)
#library(ChIPseeker)
#library(ChIPpeakAnno)
#library(rtracklayer)
#library(seqinr)
library(shinythemes)
library(shinyjs)
## Load microalgae annotation packages
# library(org.Otauri.eg.db) 
##install.packages(pkgs = "./packages/annotation_packages/org.Vcarteri.eg.db/",repos = NULL,type="source")
##install.packages(pkgs = c("org.Bprasinos.eg.db", "org.Csubellipsoidea.eg.db", "org.Dsalina.eg.db", "org.Knitens.eg.db", "org.MpusillaCCMP1545.eg.db", "org.Otauri.eg.db", "org.Smuscicola.eg.db", "org.Creinhardtii.eg.db", "org.Czofingiensis.eg.db", "org.Hlacustris.eg.db", "org.Mendlicherianum.eg.db", "org.Ngaditana.eg.db", "org.Ptricornutum.eg.db", "org.Vcarteri.eg.db"),repos = NULL,type="source")
##install.packages(pkgs = c("TxDb.Bprasinos.ORCAE", "TxDb.Czofingiensis.Phytozome", "TxDb.Knitens.Phycocosm", "TxDb.Ngaditana.JGI", "TxDb.Smuscicola.pub", "TxDb.Creinhardtii.Phytozome", "TxDb.Dsalina.Phytozome", "TxDb.Mendlicherianum.pub", "TxDb.Otauri.JGI", "TxDb.Vcarteri.Phytozome", "TxDb.Csubellipsoidea.Phytozome", "TxDb.Hlacustris.NCBI", "TxDb.MpusillaCCMP1545.Phytozome", "TxDb.Ptricornutum.Ensembl.Protists"),repos = NULL,type="source")
# library(org.MpusillaCCMP1545.eg.db)
# library(org.Bprasinos.eg.db)
# library(org.Csubellipsoidea.eg.db)
# library(org.Creinhardtii.eg.db)
# library(org.Vcarteri.eg.db)
# library(org.Dsalina.eg.db)
# library(org.Hlacustris.eg.db)
# library(org.Czofingiensis.eg.db)
# library(org.Knitens.eg.db)
# library(org.Mendlicherianum.eg.db)
# library(org.Smuscicola.eg.db)
# library(org.Ptricornutum.eg.db)
# library(org.Ngaditana.eg.db)

## Load microalgae genome annotation packages
# library(TxDb.Otauri.JGI)
# library(TxDb.MpusillaCCMP1545.Phytozome)
# library(TxDb.Bprasinos.Orcae)
# library(TxDb.Csubellipsoidea.Phytozome)
# library(TxDb.Creinhardtii.Phytozome)
# library(TxDb.Vcarteri.Phytozome)
# library(TxDb.Dsalina.Phytozome)
# library(TxDb.Hlacustris.NCBI)
# library(TxDb.Czofingiensis.Phytozome)
# library(TxDb.Knitens.Phycocosm)
# library(TxDb.Mendlicherianum.pub)
# library(TxDb.Smuscicola.pub)
# library(TxDb.Ptricornutum.Ensembl.Protists)
# library(TxDb.Ngaditana.JGI)
library(ape)
library(ggtree)
library(stringr)
library(seqinr)
library(MetBrewer)
library(glue)
library(ggplot2)
library(shiny)
library(shinythemes)
library(patchwork)

microalgae.names <- c("Ostreococcus tauri", 
                      "Micromonas pusilla CCMP1545",
                      "Bathycoccus prasinos",
                      "Coccomyxa subellipsoidea",
                      "Chlamydomonas reinhardtii", 
                      "Volvox carteri",
                      "Dunaliella salina",
                      "Haematococcus lacustris",
                      "Chromochloris zofingiensis",
                      "Klebsormidium nitens",
                      "Mesotaenium endlicherianum",
                      "Spirogloea muscicola",
                      "Phaeodactylum tricornutum",
                      "Nannochloropsis gaditana")
names(microalgae.names) <- c("otauri", 
                             "mpusilla",
                             "bprasinos",
                             "csubellipsoidea",
                             "creinhardtii", 
                             "vcarteri",
                             "dsalina",
                             "hlacustris",
                             "czofingiensis",
                             "knitens",
                             "mendlicherianum",
                             "smuscicola",
                             "ptricornutum",
                             "ngaditana")

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

## Micromonas pusilla CCMP1545 gene link to Phytozome
mpusilla.gene.link <- function(gene.name)
{
  phytozome.link <- paste0("https://phytozome-next.jgi.doe.gov/report/gene/MpusillaCCMP1545_v3_0/",gene.name)
  gene.link <- paste(c("<a href=\"",
                       phytozome.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## Coccomixa subellipsoidea gene link to Phytozome
csubellipsoidea.gene.link <- function(gene.name)
{
  phytozome.link <- paste0("https://phytozome-next.jgi.doe.gov/report/gene/CsubellipsoideaC_169_v2_0/",gene.name)
  gene.link <- paste(c("<a href=\"",
                       phytozome.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## Chlamydomonas reinhardtii gene link to Phytozome
chlamy.gene.link <- function(gene.name)
{
  phytozome.link <- paste0("https://phytozome-next.jgi.doe.gov/report/gene/Creinhardtii_v5_6/",gene.name)
  gene.link <- paste(c("<a href=\"",
                       phytozome.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## Dunaliella salina gene link to Phytozome
dsalina.gene.link <- function(gene.name)
{
  phytozome.link <- paste0("https://phytozome-next.jgi.doe.gov/report/gene/Dsalina_v1_0/",gene.name)
  gene.link <- paste(c("<a href=\"",
                       phytozome.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## Volvox carteri gene link to Phytozome
vcarteri.gene.link <- function(gene.name)
{
  phytozome.link <- paste0("https://phytozome-next.jgi.doe.gov/report/gene/Vcarteri_v2_1/",gene.name)
  gene.link <- paste(c("<a href=\"",
                       phytozome.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## Chromochloris zofingiensis gene link to Phytozome
czofingiensis.gene.link <- function(gene.name)
{
  phytozome.link <- paste0("https://phytozome-next.jgi.doe.gov/report/gene/Czofingiensis_v5_2_3_2/",gene.name)
  gene.link <- paste(c("<a href=\"",
                       phytozome.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

# Ngaditana link
ngaditana.gene.link <- function(gene.name)
{
  return("https://nannochloropsis.4ngs.com/page/listgenes")
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

## no gene link
no.gene.link <- function(gene.name)
{
  return(gene.name)
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
                     "Gene Set Functional Analysis" = "genes",
                     "Genomic Loci Functional Analysis" = "chip",
                     "MARACAS, MicroAlgae RnA-seq and Chip-seq AnalysiS" = "maracas",
                     "Funtree, Phylogenomic Analysis of Genes in Viridiplantae" = "funtree",
                     "Tutorials" = "tutorials",
                     "GitHub repository" = "github",
                     "Citation and Contact" = "citation"
                   ))),
    column(
      width = 8,
      tags$div(align = "center", 
               tags$h1(tags$b("ALGAEFUN with MARACAS"), tags$br()),
               tags$h2("microALGAE FUNctional enrichment tool for MicroAlgae RnA-seq and Chip-seq AnalysiS")),
      tags$br(),tags$br(),
      conditionalPanel(condition = "input.navigation_bar == 'home'",
                       tags$div(align = "justify", "Welcome to", tags$b("ALGAEFUN")," with ", tags$b("MARACAS"), "a microalgae web based tool for the analysis of ", 
                                tags$b("RNA-seq"), "and ", tags$b("ChIP-seq"), "data and the", tags$b("functional annotation"), "of the resulting gene sets and genomic loci. ",
                                tags$b("ALGAEFUN"), "with ", tags$b("MARACAS"), "supports the analysis for a wide collection 
               of microalgae that includes", tags$i("Chlamydomonas reinhardtii, Ostreococcus tauri, Phaeodactylum tricornutum"), "and ", 
                                tags$i("Nannochlorpsis gaditana."), "Please select from the navigation bar on the left the type of analysis you want to perform. You can also 
               see our", tags$b("video tutorials"), "on how to analyse RNA-seq and
               ChIP-seq data as well as on how to functionally annotate gene sets and genomic loci. Our
               code is freely available at", tags$b("Github."), "Please cite our work if you find it useful in your research."),
               
               tags$div(align = "justify", "Below you can find the phylogenetic relationship between the different microalgae species supported in ALGAEFUN with MARACAS. 
                        Our web tool also supports the generation of phylogenetic trees and the identification of potential orthologous genes using the functionality ",
                        tags$b("Funtree"), " available from the navigation bar on the left."),
               tags$br(),tags$br(),
               # 
               tags$div(align ="center",img(src='phylogeny.png', align = "center", width=600))
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
      
      conditionalPanel(condition = "input.navigation_bar == 'funtree'",
                       tags$div(align="justify", tags$b("Funtree"), "allows researchers to explore orthologous genes across different evolutionally 
                       distant species given a target gene. This tool is based on our phylogenomic analysis carried out using", 
                                tags$b(tags$a(href="https://github.com/davidemms/OrthoFinder", "orthofinder")),
                                "To perform the analysis, please follow the next instructions:",
                                tags$ol(
                                  tags$li("In the lower panel choose as many species as you wish to include in the tree.",),
                                  tags$li("Then write the ID of your gene of interest in the textbox, which has to correspond to one of the organisms selected
                                  in the previous step. If not, please select it before continuing. An example for the ",tags$i("Arabidopsis thaliana  CCA1"), 
                                  " gene (AT2G46830, Circadian Clock Associated 1) appears in the text box."),
                                  tags$li("Finally, click the", tags$b("Have Fun!"), "button to construct the tree. Each organism appears with an specific
                                          color detailed in legend and the target gene is highlighted in red. Also, a text box allows user to copy the names
                                          of all the detected orthologs in a simple way."))),
                       checkboxGroupInput(inputId = "selected_organisms",
                                          selected = c("mp"), 
                                          choiceNames = c("Porphyra umbilicalis","Phaeodactylum tricornutum","Nannochloropsis gaditana",
                                                          "Cyanophora paradoxa",  "Ostreococcus tauri", "Bathycoccus prasinos",
                                                          "Micromonas pusilla", "Ulva mutabilis", "Chromochloris zofingiensis",
                                                          "Raphidocelis subcapitata", "Dunaliella salina", "Chlamydomonas reinhardtii", 
                                                          "Volvox carteri", "Chlorokybus atmophyticus", "Mesostigma viride",
                                                          "Klebsormidium nitens",  "Mesotaenium endlicherianum", "Spirogloea muscicola",
                                                          "Marchantia polymorpha", "Sphangum magellanicum",  "Physcomitrium patens",
                                                          "Ceratodon purpureus", "Anthoceros agrestis",  "Selaginella moellendorffii", 
                                                          "Ceratopteris richardii", "Azolla filiculoides", "Salvinia cucullata",
                                                          "Cycas panzhihuaensis", "Thuja plicata", "Triticum aestivum", 
                                                          "Aegilops tauschii", "Oryza sativa", "Sorghum bicolor",
                                                          "Zea mays", "Solanum lycopersicum", "Arabidopsis thaliana"), 
                                          choiceValues=c("pu","pt","ng",
                                                         "cyano","ot","bp",
                                                         "mi","um", "cz",
                                                         "rs", "ds", "cr",
                                                         "vc", "ca", "mv",
                                                         "kn", "me", "sp",
                                                         "mp",  "smag", "pp",
                                                         "cp", "aa", "sm",
                                                         "cri", "af", "sc",
                                                         "cyc", "tp", "ta",
                                                         "aegi", "os", "sb",
                                                         "zm", "sl","a
                                          label= "Select the organisms to show in tree"),
                       textInput(inputId = "geneInt",value = "AT2G46830",label = NULL, placeholder = "AT2G46830"),
                       actionButton(inputId = "funtree_button",label = "Have fun!", icon("send")),
                       
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'chip'",
                       tags$div(align="justify", tags$b("AlgaeFUN"), "allows researchers to perform", tags$b("annotation analysis 
                                of genomic loci or regions."), "These are typically generated from", tags$b("ChIP-seq"), "studies 
                                of the genome-wide distribution of", tags$b("epigenetic marks or transcription factor binding sites."),
                                "Our tool", tags$b("MARACAS"), "can be used to perform this type of analysis. The set of marked genes 
                                can be obtained as well as the distribution of the genomic loci overlapping specific genes parts. Also, individual marked genes 
                                and the average signal level around the TSS (Transcription Start Site) and TES (Transcription
                                End Site) over the complete set of marked genes can be visualized.", " See our", tags$b("video tutorial"),
                                "for details or follow the next steps to perform your analysis:",
                                tags$ol(
                                  tags$li("In the left panel choose your ", tags$b("microalgae")," of interest, the gene", 
                                          tags$b("promoter length"), "and the",  tags$b("gene parts"), "that will be considered
                                          when determining the marked genes."),
                                  tags$li("Insert in the text box your ", tags$b("set of genomic regions"), " as a table consisting 
                                          of three tab-separated columns representing the chromosome, the start and end position of 
                                          the regions. An example can be loaded by clicking on the ",  tags$b("Example"), " button. 
                                          Click on", tags$b("Clear"), " button to remove the loaded gene set. Alternatively, using the",
                                          tags$b("Browse..."), "button, the genomic regions can be uploaded from a file in BED format as 
                                          described previously containing as least three columns. This file can be obained using our tool", 
                                          tags$b("MARACAS.")),
                                  tags$li("Optionally, users can upload the genome wide signal level of a epigenetic mark or transcription 
                                          factor binding in a BigWig file. This file can be obained using our tool", tags$b("MARACAS.")),
                                  tags$li("Click on the ", tags$b("Have Fun"), " button to perform the specified analysis. The
                                          results will be shown in the different tabs below.")
                                )
                       )),
      conditionalPanel(condition = "input.navigation_bar == 'maracas'",
                       tags$div(align = "justify", tags$b("MARACAS (MicroAlgae RnA-seq and Chip-seq AnalysiS) ")," 
                       is an automatic computational pipeline specifically designed for the analysis of ",  
                       tags$b("microalgae RNA-seq and ChIP-seq data"),". MARACAS starts processing raw fastq files 
                       and generates lists of differentially expressed genes for RNA-seq data and lists of genomic 
                       loci in bed format for ChIP-seq data. BigWig files containing normalized mapping signal are 
                       also generated. The analysis are performed according to user specified parameters. Reports 
                       in html and pdf format are produced for easy exploration of the results.", tags$br(),
                       "Differential expressed gene lists and genomic loci lists can be further analyzed using this 
                       online tool AlgaeFUN for functional analysis.", tags$br(),
                       tags$b("MARACAS"), " supports a wide range of microalge including:",tags$br(),
                       tags$ul(
                         tags$li("Ostreococcus tauri"),
                         tags$li("Micromonas pusilla CCMP1545"),
                         tags$li("Bathycoccus prasinos"),
                         tags$li("Coccomyxa subellipsoidea"),
                         tags$li("Chlamydomonas reinhardtii"),
                         tags$li("Volvox Carteri"),
                         tags$li("Dunaliella salina"),
                         tags$li("Haematococcus lacustris"),
                         tags$li("Chlomochloris zofingiensis"),
                         tags$li("Klebsormidium nitens"),
                         tags$li("Mesotaenium endlicherianum"),
                         tags$li("Spirogloea muscicola"),
                         tags$li("Phaeodactylum tricornutum"),
                         tags$li("Nannochloropsis gaditana"),
                       ),
                       tags$b("MARACAS"), " can be executed in a sequential mode in a laptop or server and 
                       in a distributed/parallel mode in a computer cluster.", tags$br()),
                       tags$div(align="center",
                       tags$a(href="https://github.com/fran-romero-campero/MARACAS", target="_blank",
                              tags$h4(tags$b("Click here for instructions on how to download, install and execute MARACAS.")))
                       )
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'github'",
                       tags$div(align = "justify", tags$b("AlgaeFUN,"), "is entirely developed using 
        the R package", tags$b( tags$a(href="https://shiny.rstudio.com/", "shiny.")), "The 
        source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                tags$b("GitHub."), "If you experience any problem using AlgaeFUN please create an", 
                                tags$b(tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN/issues","issue")), 
                                "in GitHub and we will address it."),
                       tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN",
                                                                     "AlgaeFUN at GitHub")))),
                       tags$br(),
                       tags$br(),
                       
                       tags$div(align = "justify", tags$b("MARACAS,"), "is developed using bash scripting and several 
        bioconductor R packages. The source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                tags$b("GitHub."), "If you experience any problem using AlgaeFUN please create an", 
                                tags$b(tags$a(href="https://github.com/fran-romero-campero/MARACAS/issues","issue")), 
                                "in GitHub and we will address it."),
                       tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/MARACAS",
                                                                     "MARACAS at GitHub")))),
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'citation'",
                       tags$div(align = "justify", "We are strongly committed to", tags$b("open access software"), 
                                "and", tags$b("open science."),"Following our philosophy we have deposited our GitHub code 
                       into", tags$a(href="https://zenodo.org/record/4754516#.YJxLPSaxUws", target="_blank",tags$b("Zenodo")), ", a
                       general-purpose open-access repository developed under the", 
                                tags$a(href="https://www.openaire.eu/", target="_blank", tags$b("European OpenAIRE program.")), "If you use AlgaeFUN with
                                MARACAS in your research please cite us:",
                                tags$br(),
                                tags$br(),
                                tags$div(tags$a(href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04639-5", target="_blank",
                                tags$b("Romero-Losada, A.B., Arvanitidou, C., de los Reyes, P., García-González, M. & Romero-Campero, F.J. 
                                ALGAEFUN with MARACAS, microALGAE FUNctional enrichment tool for MicroAlgae RnA-seq and Chip-seq AnalysiS. 
                                BMC Bioinformatics 23, 113 (2022). https://doi.org/10.1186/s12859-022-04639-5")))),
                       tags$br(),
                       tags$br(),
                       tags$div("If you want to contact us send an email to:", img(src='email.png',width=100)),
                       tags$br()
                       
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'tutorials'",
                       tags$div(align="center",uiOutput("video_tutorial")),
                       tags$div(align = "justify", 
                                tags$br(),
                                tags$br(),
                                tags$div(tags$h4(tags$b("Above you can find a video tutorial on how to use the different tools implemented 
                                in AlgaeFUN with MARACAS."))))
                       
      ),
      
    ),
    column(
      width = 2,
      img(src='logo_ibvf.jpg', align = "center", width=100),
      img(src='logo_us.png', align = "center", width=100),
      tags$br(),tags$br(),tags$br(),
      img(src='logo_csic.jpg', align = "center", width=100),
      tags$br(),tags$br(),
      tags$div(align="center",width=60,
               HTML("<script type=\"text/javascript\" src=\"//rf.revolvermaps.com/0/0/8.js?i=5jamj0c2y0z&amp;m=7&amp;c=ff0000&amp;cr1=ffffff&amp;f=arial&amp;l=33\" async=\"async\"></script>"))
    )
  ),

  
  tags$br(),tags$br(),
  
  #Interface where the user can choose his/her preferencies, separated by columns
  fluidRow(
      column(width = 4,

      conditionalPanel(condition = "input.navigation_bar == 'genes' || input.navigation_bar == 'chip'",
        #Choose the target microalgae
        selectInput(inputId = "microalgae", label="Choose your favourite microalgae", 

                  choices=c("Ostreococcus tauri" = "otauri", 
                            "Micromonas pusilla CCMP1545" = "mpusilla",
                            "Bathycoccus prasinos" = "bprasinos",
                            "Coccomyxa subellipsoidea" = "csubellipsoidea",
                            "Chlamydomonas reinhardtii" = "creinhardtii", 
                            "Volvox carteri" = "vcarteri",
                            "Dunaliella salina" = "dsalina",
                            "Haematococcus lacustris" = "hlacustris",
                            "Chromochloris zofingiensis" = "czofingiensis",
                            "Klebsormidium nitens" = "knitens",
                            "Mesotaenium endlicherianum" = "mendlicherianum",
                            "Spirogloea muscicola" = "smuscicola",
                            "Phaeodactylum tricornutum" = "ptricornutum",
                            "Nannochloropsis gaditana" = "ngaditana"
                            ))),

     conditionalPanel(condition = "input.navigation_bar == 'genes'",    
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
      conditionalPanel(condition = "input.navigation_bar == 'genes'",
                       
        numericInput(inputId = "pvalue", 
                     label= "Which will be your chosen p-value?", 
                     value= 0.05)
      ),
      
       
      #Choose a promoter length
      conditionalPanel(condition = "input.navigation_bar == 'chip'",
                       sliderInput(inputId = "promoter_length", 
                                    label= "Choose the distance in base pairs around the 
                                    Transcriptional Start Site defining gene promoters", 
                                    min=100, max=2000,value=1000,step=100)),
      
      #Choose genomic features
      conditionalPanel(condition = "input.navigation_bar == 'chip'",
                       checkboxGroupInput(inputId = "selected_genomic_features",selected = c("Promoter","5' UTR","Exon","Intron","3' UTR"), 
                                   label= "A gene will be associated to an input genomic locus
                                   when it overlaps one of the following gene features:",
                                   choices = c("Promoter","5' UTR","Exon","Intron","3' UTR")))
      
      ),
      
      column( width = 8,
        #UI for functional enrichment over a gene set obtained for instance from an RNAseq data analysis
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
        #UI for functional annotation of genomic loci obtained for instance from a ChIPseq data analysis 
        textAreaInput(inputId = "genomic_regions", 
                      label= "Insert a set of genomic regions", 
                      width="200%", height = "200px", 
                      placeholder = "Insert set of genomic regions",
                      value= ""),
        actionButton(inputId = "example_genomic_regions",label = "Example"),
        actionButton(inputId = "clear_genomic_regions",label = "Clear"),
        fileInput(inputId = "genomic_regions_file",label = "Choose File with the Genomic Regions to Upload",width = "100%"),
        fileInput(inputId = "bw_file",label = "Choose BigWig File to Upload for Profile Representations: (Optional)", width= "100%"),
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
      #KEGG pathway maps for gene set enrichment analysis

      conditionalPanel(condition = "(input.navigation_bar == 'genes')",
         tabsetPanel(type ="tabs",
            tabPanel(tags$b("GO ENRICHMENT"),
                  shinyjs::useShinyjs(),
                  hidden(div(id='loading.enrichment.go',h3('Please be patient, computing GO enrichment ...'))), 
                  hidden(div(id='ready.enrichment.go',h3('Your GO enrichment is ready!'))), 
                  htmlOutput(outputId = "gene_sanity_go"),
                  htmlOutput(outputId = "wrong_genes_go"),
                  tags$br(),
                  htmlOutput(outputId = "intro_go"),
                  tags$br(),
                  tags$br(),
                  tabsetPanel(type = "tabs",
                     tabPanel(tags$b("GO Enrichment Table"),
                              tags$br(),
                              htmlOutput(outputId = "textGOTable"),
                              tags$br(), tags$br(),
                              dataTableOutput(outputId = "output_go_table"),
                              uiOutput(outputId = "download_ui_for_go_table")#,
                              #htmlOutput(outputId = "revigo")
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
                    # tabPanel(tags$b("GO Emap"),
                    #          tags$br(),
                    #          htmlOutput(outputId = "emapplot_text"),
                    #          tags$br(),
                    #          div(style= "text-align: center;",
                    #              plotOutput(outputId = "emap.plot",inline=TRUE)),
                    #          tags$br(),
                    #          tags$br(), tags$br()),
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
                     tags$br(),
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
                     )#closes tabset panel for KEGG results
                     
            ),#closes tabs panel for kegg enrichment
            tabPanel(tags$b("COMMON IDs CONVERSION"),
                     shinyjs::useShinyjs(),
                     hidden(div(id='loading.conversion',h3('Please be patient, computing common IDs conversion ...'))), 
                     hidden(div(id='ready.conversion',h3('Your common IDs conversion is ready!'))), 
                     htmlOutput(outputId = "gene_sanity_conversion"),
                     htmlOutput(outputId = "wrong_genes_conversion"),
                     tags$br(),
                     htmlOutput(outputId = "intro_conversion"),
                     tags$br(),
                     tags$br(),
                     htmlOutput(outputId = "textconversionTable"),
                     tags$br(), tags$br(),
                     dataTableOutput(outputId = "output_conversion_table"),
                     uiOutput(outputId = "download_ui_for_conversion_table")
            ) # close tabPanel for common IDs conversion
            )#closes tab panel for GO, KEGG enrichment and common id conversion
      ), #closes conditional panel for gene option

      ## Main panel containing the results organized in different tabs for gene 
      ## genomic loci functional annotation
      
      conditionalPanel(condition = "input.navigation_bar == 'funtree'",
                       tabsetPanel(type = "tabs",
                          tabPanel(tags$b("Genes Tree"),
                               tags$br(), 
                               verbatimTextOutput(outputId = "treeTips"), 
                               tags$br(),
                               uiOutput(outputId = "download_tips"),
                               tags$br(),
                               uiOutput(outputId = "download_tree"),
                               tags$br(),
                               plotOutput("treePlot")
                               
                               )
                       )),
      
      conditionalPanel(condition = "input.navigation_bar == 'chip'",
                       hidden(div(id='loading.chip',h3('Please be patient, computing genomic loci analysis ...'))), 
                       hidden(div(id='ready.chip',h3('Your genomic loci analysis is ready!'))),
                       tabsetPanel(type = "tabs",
                           tabPanel(tags$b("Marked Genes Table"),
                               tags$br(), 
                               htmlOutput(outputId = "textTableAnnotatedGenes"),
                               tags$br(), 
                               dataTableOutput(outputId = "output_gene_chip_table"),
                               uiOutput(outputId = "download_gene_chip_table"),
                               tags$br(), tags$br()),
                           tabPanel(tags$b("Overlapping Gene Parts Distribution"),
                               tags$br(),
                               htmlOutput(outputId = "piechart.text"),
                               div(style= "text-align: center;",
                                   plotOutput(outputId = "annotation.pie.chart",inline=TRUE))), 
                           tabPanel(tags$b("Distance to TSS Visualization"),
                               tags$br(),
                               htmlOutput(outputId = "tss.distance.text"),      
                               tags$br(), 
                               plotOutput(outputId = "distance.to.tss",inline=TRUE),
                               tags$br(), tags$br()),
                           tabPanel(tags$b("Individual Genes Visualization"),
                                    tags$br(), 
                                    htmlOutput(outputId = "gene.signal.text"),
                                    tags$br(),
                                    uiOutput(outputId = "annotated_genes"),
                                    tags$br(), tags$br()),
                           tabPanel(tags$b("TSS/TES signal visualization"),
                               hidden(div(id='loading.tss.signal',h3('Please be patient, computing singal around TSS ...'))), 
                               hidden(div(id='ready.tss.signal',h3('Your signal around TSS is computed!'))),
                               tags$br(),
                               htmlOutput(outputId = "tss.signal.text"),      
                               tags$br(), 
                               plotOutput(outputId = "tss_signal"),
                               tags$br(), tags$br())
                           )
      )
      )
))

server <- shinyServer(function(input, output, session) {
  
  ## video tutorial
  output$video_tutorial <- renderUI({
    HTML("<iframe width=\"560\" height=\"315\" src=\"https://www.youtube.com/embed/ZCWrqOxrdJM\" title=\"YouTube video player\" frameborder=\"0\" allow=\"accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture\" allowfullscreen></iframe>")
  })  
  
  ## Clear content of gene set text area and previous results
  observeEvent(input$clear_gene_set, {
    updateTextAreaInput(session=session, inputId = "genes",value = "")
    
    shinyjs::hideElement(id = 'loading.enrichment.go')
    shinyjs::hideElement(id = 'ready.enrichment.go')
    output$intro_go <- renderText(expr = "")
    output$gene_sanity_go <- renderText(expr = "")
    output$wrong_genes_go  <- renderText(expr = "")
    output$textGOTable <- renderText(expr = "")
    output$output_go_table <- renderDataTable(expr = NULL)
    output$download_ui_for_go_table<- renderUI(expr = NULL)
    #output$revigo<- renderUI(expr = NULL)
    output$go_graph <- renderText(expr = "")
    output$go.plot <- renderPlot(expr = NULL)
    output$barplot_text <- renderText("")
    output$bar.plot <- renderPlot(expr = NULL)
    output$dotplot_text <- renderText("")
    output$dot.plot <- renderPlot(expr = NULL)
    #output$emapplot_text <- renderText("")
    #output$emap.plot <- renderPlot(expr = NULL)
    output$cnetplot_text <- renderText("")
    output$cnet.plot <- renderPlot(expr = NULL)
    
    shinyjs::hideElement(id = 'ready.enrichment.kegg')
    shinyjs::hideElement(id = 'loading.enrichment.kegg')
    output$intro_kegg <- renderText(expr = "")
    output$textKEGGTable <- renderText(expr = "")
    output$output_pathway_table <- renderDataTable(expr = NULL)
    output$textKEGGImage <- renderText(expr = "")
    output$kegg_selectize <- renderUI(expr = NULL)
    output$kegg_image <- renderImage(expr = NULL,deleteFile = T)
    output$text_module_kegg <- renderText(expr = "")
    output$output_module_table <- renderDataTable(expr = NULL)
    
    shinyjs::hideElement(id = 'ready.conversion')
    shinyjs::hideElement(id = 'loading.conversion')
    output$intro_conversion <- renderText(expr = "")
    output$textconversionTable <- renderText(expr = "")
    output$output_conversion_table <- renderDataTable(expr = NULL)
    
  })
  
  ## Clear content of genomic regions set text area
  observeEvent(input$clear_genomic_regions, {
    updateTextAreaInput(session=session, inputId = "genomic_regions",value = "")
    output$textTableAnnotatedGenes <- renderText(expr = "")
    output$output_gene_chip_table <- renderDataTable(expr = NULL)
    output$download_gene_chip_table <- renderUI(expr = NULL)
    output$piechart.text <- renderText(expr = "")
    output$annotation.pie.chart <- renderPlot(expr = NULL)
    output$tss.distance.text <- renderText(expr = "")
    output$distance.to.tss <- renderPlot(expr = NULL)
    output$gene.signal.text <- renderText(expr = "")
    output$annotated_genes <- renderUI(expr = NULL)
    output$tss.signal.text <- renderText(expr = "")
    output$tss_signal <- renderPlot(expr = NULL)
    shinyjs::hideElement(id = 'loading.tss.signal')
    shinyjs::hideElement(id = 'ready.tss.signal')
    shinyjs::hideElement(id = 'loading.chip')
    shinyjs::hideElement(id = 'ready.chip')
  })

  ## Clear content of universe set text area
  observeEvent(input$clear_universe_set, {
    updateTextAreaInput(session=session, inputId = "background",value = "")
  })
  
  ## Add an example of gene set to text area
  observeEvent(input$example_genes, {
    example.file <- paste(c("example_files/example_",input$microalgae,".txt"),collapse="")
    example.genes <- read.table(file = example.file,header = F,as.is = T,comment.char="",quote="\"")[[1]]
    updateTextAreaInput(session=session, inputId = "genes",value = paste(example.genes,collapse="\n"))
  })
 
  ## Add an example of genomic regions to text area
  observeEvent(input$example_genomic_regions, {
    example.file <- paste(c("example_files/example_genomic_regions_",input$microalgae,".txt"),collapse="")
    example.genomic.regions <- read.table(file = example.file,header = F,as.is = T)
    example.text <- NULL
    for(i in 1:nrow(example.genomic.regions))
    {
      example.text <- paste(example.text,paste(example.genomic.regions[i,],collapse = "\t"),sep="\n")
    }
    #print(example.text)
    updateTextAreaInput(session=session, inputId = "genomic_regions",value = example.text)
    bw.file <- paste(c("example_files/example_",input$microalgae,".bw"),collapse="")
    if(file.exists(bw.file))
    {
      selected.bigwig.files <- bw.file
    }
  })
  
  ## Actions to perform after click the go button
  observeEvent(input$go.button , {
    shinyjs::showElement(id = 'loading.enrichment.go')
    shinyjs::hideElement(id = 'ready.enrichment.go')
    
    # Remove previous results
    output$intro_go <- renderText(expr = "")
    output$gene_sanity_go <- renderText(expr = "")
    output$wrong_genes_go  <- renderText(expr = "")
    output$textGOTable <- renderText(expr = "")
    output$output_go_table <- renderDataTable("")
    output$download_ui_for_go_table<- renderUI(expr = NULL)
    #output$revigo<- renderUI(expr = NULL)
    output$go_graph <- renderText(expr = "")
    output$go.plot <- renderPlot(expr = NULL)
    output$barplot_text <- renderText("")
    output$bar.plot <- renderPlot(expr = NULL)
    output$dotplot_text <- renderText("")
    output$dot.plot <- renderPlot(expr = NULL)
    #output$emapplot_text <- renderText("")
    #output$emap.plot <- renderPlot(expr = NULL)
    output$cnetplot_text <- renderText("")
    output$cnet.plot <- renderPlot(expr = NULL)
    
    output$intro_kegg <- renderText(expr = "")
    output$output_pathway_table <- renderDataTable(expr = NULL)
    output$textKEGGImage <- renderText(expr = "")
    output$kegg_selectize <- renderUI(expr = NULL)
    output$kegg_image <- renderImage(expr = NULL,deleteFile = T)
    output$text_module_kegg <- renderText(expr = "")
    output$output_module_table <- renderDataTable(expr = NULL)
    
    shinyjs::hideElement(id = 'ready.conversion')
    shinyjs::hideElement(id = 'loading.conversion')
    output$intro_conversion <- renderText(expr = "")
    output$textconversionTable <- renderText(expr = "")
    output$output_conversion_table <- renderDataTable(expr = NULL)
    
    # Load libraries
    library(clusterProfiler)
    library(pathview)
      
    ## Select org.Db 
    if(input$microalgae == "otauri")
    {
      library(org.Otauri.eg.db)
      org.db <- org.Otauri.eg.db
      microalgae.genes <- read.table(file = "universe/otauri_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- ostta.gene.link
    } else if (input$microalgae == "mpusilla")
    {
      library(org.MpusillaCCMP1545.eg.db)
      org.db <- org.MpusillaCCMP1545.eg.db
      microalgae.genes <- read.table(file = "universe/mpusilla_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- mpusilla.gene.link
    } else if (input$microalgae == "bprasinos")
    {
      library(org.Bprasinos.eg.db)
      org.db <- org.Bprasinos.eg.db
      microalgae.genes <- read.table(file = "universe/bathy_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- bathy.gene.link
    } else if (input$microalgae == "csubellipsoidea")
    {
      library(org.Csubellipsoidea.eg.db)
      org.db <- org.Csubellipsoidea.eg.db
      microalgae.genes <- read.table(file = "universe/cocsu_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- csubellipsoidea.gene.link
    } else if (input$microalgae == "creinhardtii")
    {
      library(org.Creinhardtii.eg.db)
      org.db <- org.Creinhardtii.eg.db
      microalgae.genes <- read.table(file = "universe/cre_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- chlamy.gene.link
    } else if (input$microalgae == "vcarteri")
    {
      library(org.Vcarteri.eg.db)
      org.db <- org.Vcarteri.eg.db
      microalgae.genes <- read.table(file = "universe/vocar_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- vcarteri.gene.link
    } else if (input$microalgae == "dsalina")
    {
      library(org.Dsalina.eg.db)
      org.db <- org.Dsalina.eg.db
      microalgae.genes <- read.table(file = "universe/dusal_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- dsalina.gene.link
    } else if (input$microalgae == "hlacustris")
    {
      library(org.Hlacustris.eg.db)
      org.db <- org.Hlacustris.eg.db
      microalgae.genes <- read.table(file = "universe/hlacustris_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- ncbi.gene.link
    } else if (input$microalgae == "czofingiensis")
    {
      library(org.Czofingiensis.eg.db)
      org.db <- org.Czofingiensis.eg.db
      microalgae.genes <- read.table(file = "universe/zofi_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- czofingiensis.gene.link
    } else if (input$microalgae == "knitens")
    {
      library(org.Knitens.eg.db)
      org.db <- org.Knitens.eg.db
      microalgae.genes <- read.table(file = "universe/klebsor_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- no.gene.link
    } else if (input$microalgae == "mendlicherianum")
    {
      library(org.Mendlicherianum.eg.db)
      org.db <- org.Mendlicherianum.eg.db
      microalgae.genes <- read.table(file = "universe/mesotaenium_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- no.gene.link
    } else if (input$microalgae == "smuscicola")
    {
      library(org.Smuscicola.eg.db)
      org.db <- org.Smuscicola.eg.db
      microalgae.genes <- read.table(file = "universe/smuscicola_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- no.gene.link
    } else if (input$microalgae == "ptricornutum")
    {
      library(org.Ptricornutum.eg.db)
      org.db <- org.Ptricornutum.eg.db
      microalgae.genes <- read.table(file = "universe/phatri_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- phaeodactylum.gene.link
    } else if (input$microalgae == "ngaditana")
    {
      library(org.Ngaditana.eg.db)
      org.db <- org.Ngaditana.eg.db
      microalgae.genes <- read.table(file = "universe/naga_universe.txt",as.is = T, comment.char = "", quote="\"")[[1]]
      gene.link.function <- ngaditana.gene.link
    }  
    
    ## Extract genes from text box or uploaded file
    if(is.null(input$gene_set_file))
    {
      target.genes <- as.vector(unlist(strsplit(input$genes, split="\n",
                                                fixed = TRUE)[1]))
    } else
    {
      target.genes <- read.table(file=input$gene_set_file$datapath, header = F,as.is = TRUE,comment.char = "")[[1]]
    }

    ## Select gene universe
    if(input$input_mode == "No")
    {
      gene.universe <- microalgae.genes #unique(select(org.db,columns = c("GO"),keys=keys(org.db,keytype = "GID"))[["GID"]])
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
      output$gene_sanity_conversion <-renderText(expr = gene.sanity.text)
      output$wrong_genes_conversion <- renderText(expr = paste(wrong.genes,collapse="\n"))
    }

    ## GO term enirchment analysis
    if(input$analysis == "kegg")
    {
      output$intro_go <- renderText(expr = "<p style=\"color:blue\"><b> You have chosen only to perform a KEGG enrichment analysis.
                                    Please check the content of the KEGG ENRICHMENT tab.</b></p>")
    } else if((input$analysis == "go" || input$analysis == "both") && (length(target.genes) == 0))
    {
      output$intro_go <- renderText(expr = "<p style=\"color:red\"><b> You forgot to input a gene set. Please, paste a gene set in the text box above or
                                    select a text file containing your gene set of interest.</b></p>")
    } else if((input$analysis == "go" || input$analysis == "both") && (length(wrong.genes) == 0) && (length(target.genes) > 0))
    {
      ## Intro text for GO enrichment
      go.intro.text <- paste(c("The tabs below present the results from the <b>GO enrichment analysis</b> 
                                      performed over the input target genes (<i>",
                               paste(target.genes[1:3],collapse=" "),
                               "</i> ...) from the microalgae <b> <i>", microalgae.names[input$microalgae],
                               " </i> </b> with ", universe.text),collapse="") 
      
      output$intro_go <- renderText(expr = go.intro.text)

      ## Perform GO enrichment
      enrich.go <- enrichGO(gene          = target.genes,
                            universe      = gene.universe,
                            OrgDb         = org.db,
                            ont           = input$ontology,
                            pAdjustMethod = "BH",
                            pvalueCutoff  = input$pvalue,
                            qvalueCutoff = 1,
                            readable      = TRUE,
                            keyType = "GID")
      
      
      ## Generate ouput table
      enrich.go.result <- as.data.frame(enrich.go@result)
      enrich.go.result <- subset(enrich.go.result, pvalue < input$pvalue)

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
          tagList(downloadButton(outputId= "downloadGOTable", "Download GO Enrichment Table"),tags$br(),tags$br())
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
        # revigo.data <- paste(revigo.data <- apply(go.result.table[,c("GO ID", "q-value")], 1, paste, collapse = " "), collapse="\n")
        # 
        # url1 <- tags$a("here", href="#", onclick="document.revigoForm.submit();")
        # url2 <- tags$form(
        #   name="revigoForm", action="http://revigo.irb.hr/", method="post", target="_blank",
        #   tags$textarea(name="inputGoList", rows="1", cols="8", class="revigoText",
        #                 style="visibility: hidden", revigo.data)
        # )
        # 
        # output$revigo<- renderUI(
        #   tagList("The enriched GO terms above may be redundant. Visualize these results in REViGO in order to remove redundancy. Click", url1, url2)
        # )

        go.graph.text <- "The following acyclic graph represents the GO term enrichment
        in the target gene set. Each node stands for a GO term. The color of each node
        indicates the level of significance from grey, non-significant, to intense red,
        highly significant. An arrow is drawn from GO term A to GO term B when A is a more
        general GO term than B or B is more specific than A. Right click on the image to download it."
        
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
        of significance from blue, less significant, to red, more significant. Right click on the image to download it.")
        
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
        of significance from blue, less significant, to red, more significant. Right click on the image to download it.")
        
        ## Dotplot
        output$dot.plot <- renderPlot(
          width     = 870,
          height    = 600,
          res       = 120,
          expr = {
            dotplot(enrich.go)
          })

#         output$emapplot_text <- renderText("The following figure consists of an enrichment map where nodes represent enriched GO terms. The
#         size of a node is proportional to the number of genes annotated with the corresponding GO term in the target set.
# The node colors represent the level of significance from less signficant in blue to more significant in red. Edges are drawn
# between two nodes when the corresponding GO terms are semantically related. Right click on the image to download it.")
#         
#         ##EMAP plot
#         output$emap.plot <- renderPlot(
#           width     = 870,
#           height    = 600,
#           res       = 120,
#           expr = {
#             emapplot(enrich.go)
#           })
        
        output$cnetplot_text <- renderText("The following figure corresponds to a gene-concept network. The beige
nodes represents GO terms and the grey nodes genes. An edge is drawn from a gene to a GO term when the gene is annotated
with the corresponding gene. The size of nodes representing GO terms is proportional to the number of genes annotated
with the corresponding GO term. Right click on the image to download it.")
        
        
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
        output$textGOTable <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")
        output$go_graph <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")
        output$barplot_text <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")
        output$dotplot_text <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")
        #output$emapplot_text <- renderText(expr = "<b>No GO term enrichment detected 
        #                                 in the input gene set.</b>")        
        output$cnetplot_text <- renderText(expr = "<b>No GO term enrichment detected 
                                         in the input gene set.</b>")        
      }
      
      shinyjs::hideElement(id = 'loading.enrichment.go')
      shinyjs::showElement(id = 'ready.enrichment.go')
    }
    
    ## KEGG pathways enrichment analysis
    if(input$analysis == "go")
    {
      output$intro_kegg <- renderText(expr = "<p style=\"color:blue\"><b> You have chosen only to perform a GO enrichment analysis.
                                    Please check the content of the GO ENRICHMENT tab.</b></p>")
    } else if((input$analysis == "kegg" || input$analysis == "both") && (length(target.genes) == 0))
    {
      output$intro_kegg <- renderText(expr = "<p style=\"color:red\"><b> You forgot to input a gene set. Please, paste a gene set in the text box above or
                                    select a text file containing your gene set of interest.</b></p>")
    } else if( (input$analysis == "kegg"  || input$analysis == "both") && (length(wrong.genes) == 0))
    {
      shinyjs::showElement(id = 'loading.enrichment.kegg')
      shinyjs::hideElement(id = 'ready.enrichment.kegg')
      
      ## Update target genes and universe depending on the microalgae
      if(input$microalgae == "otauri")
      {
        target.genes <- paste0("OT_",target.genes)
        gene.universe <- paste0("OT_",gene.universe)
        organism.id <- "ota"
      } else if (input$microalgae == "mpusilla")
      {
        mpusilla.ko <- AnnotationDbi::select(org.MpusillaCCMP1545.eg.db,columns = c("KO"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID"))
        ko.universe <- mpusilla.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(mpusilla.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05)) #input$pvalue))
        
        for(i in 1:nrow(pathway.enrichment))
        {
          current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
          
          current.genes <- c()
          for(j in 1:length(current.Ks))
          {
            current.genes <- c(current.genes,subset(mpusilla.ko, KO == current.Ks[j])$GID)
          }
          
          pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
        }
      } else if(input$microalgae == "bprasinos")
      {
        gene.universe <- AnnotationDbi::select(org.Bprasinos.eg.db,columns = c("GID"),keys=keys(org.Bprasinos.eg.db,keytype = "GID"))[[1]]

        organism.id <- "bpg"
      } else if(input$microalgae == "vcarteri")
      {
        vocar.volcadraft.map <- AnnotationDbi::select(org.Vcarteri.eg.db,columns = c("VOLCADRAFT"),keys=keys(org.Vcarteri.eg.db,keytype = "GID"))
        vocar.ids <- vocar.volcadraft.map$GID
        volcadraft.ids <- vocar.volcadraft.map$VOLCADRAFT
        names(volcadraft.ids) <- vocar.ids
        names(vocar.ids) <- volcadraft.ids
        
        target.genes <- volcadraft.ids[target.genes]
        names(target.genes) <- NULL
        
        gene.universe <- volcadraft.ids[gene.universe]
        names(gene.universe) <- NULL
        
        organism.id <- "vcn"
      } else if (input$microalgae == "csubellipsoidea")
      {
        csubellipsoidea.ko <- AnnotationDbi::select(org.Csubellipsoidea.eg.db,columns = c("KO"),keys=keys(org.Csubellipsoidea.eg.db,keytype = "GID"))
        ko.universe <- csubellipsoidea.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(csubellipsoidea.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05)) #input$pvalue))
        
        for(i in 1:nrow(pathway.enrichment))
        {
          current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
          
          current.genes <- c()
          for(j in 1:length(current.Ks))
          {
            current.genes <- c(current.genes,subset(csubellipsoidea.ko, KO == current.Ks[j])$GID)
          }
          
          pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
        }
      } else if(input$microalgae == "creinhardtii")
      {
        cre.chlredraft.map <- AnnotationDbi::select(org.Creinhardtii.eg.db,columns = c("CHLREDRAFT"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID"))
        cre.ids <- cre.chlredraft.map$GID
        chlredraft.ids <- cre.chlredraft.map$CHLREDRAFT
        names(chlredraft.ids) <- cre.ids
        names(cre.ids) <- chlredraft.ids
        
        target.genes <- chlredraft.ids[target.genes]
        names(target.genes) <- NULL
        
        gene.universe <- chlredraft.ids[gene.universe]
        names(gene.universe) <- NULL
        
        organism.id <- "cre"
      } else if(input$microalgae == "dsalina")
      {
        dsalina.ko <- AnnotationDbi::select(org.Dsalina.eg.db,columns = c("KO"),keys=keys(org.Dsalina.eg.db,keytype = "GID"))
        ko.universe <- dsalina.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(dsalina.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05)) #input$pvalue))
        
        for(i in 1:nrow(pathway.enrichment))
        {
          current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
          
          current.genes <- c()
          for(j in 1:length(current.Ks))
          {
            current.genes <- c(current.genes,subset(dsalina.ko, KO == current.Ks[j])$GID)
          }
          
          pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
        }
      } else if(input$microalgae == "ptricornutum")
      {
        phatri.draft.map <- AnnotationDbi::select(org.Ptricornutum.eg.db,columns = c("PHATRIDRAFT"),keys=keys(org.Ptricornutum.eg.db,keytype = "GID"))
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
        naga.draft.map <- AnnotationDbi::select(org.Ngaditana.eg.db,columns = c("NAGADRAFT"),keys=keys(org.Ngaditana.eg.db,keytype = "GID"))
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
        knitens.ko <- AnnotationDbi::select(org.Knitens.eg.db,columns = c("KO"),keys=keys(org.Knitens.eg.db,keytype = "GID"))
        ko.universe <- knitens.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(knitens.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05)) #input$pvalue))
        
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
      } else if(input$microalgae == "mendlicherianum")
      {
        mendlicherianum.ko <- AnnotationDbi::select(org.Mendlicherianum.eg.db,columns = c("KO"),keys=keys(org.Mendlicherianum.eg.db,keytype = "GID"))
        ko.universe <- mendlicherianum.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(mendlicherianum.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05)) #input$pvalue))
        
        for(i in 1:nrow(pathway.enrichment))
        {
          current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
          
          current.genes <- c()
          for(j in 1:length(current.Ks))
          {
            current.genes <- c(current.genes,subset(mendlicherianum.ko, KO == current.Ks[j])$GID)
          }
          
          pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
        }
      } else if(input$microalgae == "smuscicola")
      {
        smuscicola.ko <- AnnotationDbi::select(org.Smuscicola.eg.db,columns = c("KO"),keys=keys(org.Smuscicola.eg.db,keytype = "GID"))
        ko.universe <- smuscicola.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(smuscicola.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05)) #input$pvalue))
        
        for(i in 1:nrow(pathway.enrichment))
        {
          current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
          
          current.genes <- c()
          for(j in 1:length(current.Ks))
          {
            current.genes <- c(current.genes,subset(smuscicola.ko, KO == current.Ks[j])$GID)
          }
          
          pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
        }
      } else if(input$microalgae == "hlacustris")
      {
        hlacustris.ko <- AnnotationDbi::select(org.Hlacustris.eg.db,columns = c("KO"),keys=keys(org.Hlacustris.eg.db,keytype = "GID"))
        ko.universe <- hlacustris.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(hlacustris.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05)) #input$pvalue))
        
        if(nrow(pathway.enrichment) > 1)
        {
          for(i in 1:nrow(pathway.enrichment))
          {
            current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
            
            current.genes <- c()
            for(j in 1:length(current.Ks))
            {
              current.genes <- c(current.genes,subset(hlacustris.ko, KO == current.Ks[j])$GID)
            }
            
            pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
          }
        }
      } else if(input$microalgae == "czofingiensis")
      {
        zofi.ko <- AnnotationDbi::select(org.Czofingiensis.eg.db,columns = c("KO"),keys=keys(org.Czofingiensis.eg.db,keytype = "GID"))
        ko.universe <- zofi.ko$KO
        ko.universe <- ko.universe[!is.na(ko.universe)]
        
        target.ko <- subset(zofi.ko,GID %in% target.genes)$KO
        target.ko <- target.ko[!is.na(target.ko)]
        
        pathway.enrichment <- as.data.frame(enrichKEGG(gene = target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05 )) #input$pvalue))
        
        for(i in 1:nrow(pathway.enrichment))
        {
          current.Ks <- strsplit(pathway.enrichment$geneID[i],split="/")[[1]]
          
          current.genes <- c()
          for(j in 1:length(current.Ks))
          {
            current.genes <- c(current.genes,subset(zofi.ko, KO == current.Ks[j])$GID)
          }
          
          pathway.enrichment$geneID[i] <- paste(intersect(unique(current.genes),target.genes),collapse="/")
        }
      }         
      
      ## Compute KEGG pathway enrichment
      if (input$microalgae != "hlacustris" && input$microalgae != "knitens" && 
          input$microalgae != "czofingiensis" && input$microalgae != "mpusilla" &&
          input$microalgae != "mendlicherianum" && input$microalgae != "smuscicola" &&
          input$microalgae != "dsalina" && input$microalgae != "csubellipsoidea")
      {
        pathway.enrichment <- enrichKEGG(gene = target.genes, organism = organism.id, keyType = "kegg",
                                         universe = gene.universe,pvalueCutoff = input$pvalue, qvalueCutoff = 1)
      }
      shinyjs::showElement(id = 'ready.enrichment.kegg')
      shinyjs::hideElement(id = 'loading.enrichment.kegg')
      
      pathway.enrichment.result <- as.data.frame(pathway.enrichment@result)
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
        
        ## Generate gene lists for enriched pathways
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
        } else if (input$microalgae == "hlacustris")
        {
          kegg.enriched.genes <- pathway.enrichment.result$geneID
          for(i in 1:length(kegg.enriched.genes))
          {
            kegg.enriched.genes[i] <- paste(strsplit(kegg.enriched.genes[i],split="/")[[1]],collapse=" ")
          }
        } else if (input$microalgae == "czofingiensis" | 
                  input$microalgae == "mpusilla" |
                  input$microalgae == "bprasinos" |
                  input$microalgae == "mendlicherianum" | 
                  input$microalgae == "smuscicola" |
                  input$microalgae == "dsalina" |
                  input$microalgae == "csubellipsoidea")
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
        # if(input$microalgae == "otauri")
        # {
        #   gene.link.function <- ostta.gene.link
        # } else if(input$microalgae == "creinhardtii" | input$microalgae == "vcarteri" | 
        #               input$microalgae == "csubellipsoidea" | input$microalgae == "dsalina" |
        #               input$microalgae == "mpusilla")
        # {
        #   gene.link.function <- phytozome.gene.link
        # } else if(input$microalgae == "bprasinos")
        # {
        #   gene.link.function <- bathy.gene.link
        # } else if(input$microalgae == "ptricornutum")
        # {
        #   gene.link.function <- phaeodactylum.gene.link
        # } else if(input$microalgae == "ngaditana")
        # {
        #   gene.link.function <- ngaditana.gene.link
        # } else if(input$microalgae == "knitens")
        # {
        #   gene.link.function <- knitens.gene.link
        # }else if(input$microalgae == "hlacustris")
        # {
        #   gene.link.function <- ncbi.gene.link
        # }else if(input$microalgae == "czofi")
        # {
        #   gene.link.function <- zofi.gene.link
        # }
        
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
        
        
        ## Figures for KEGG pathway enrichment analysis
        
        ## Prepare gene set for representation
        if(input$microalgae == "knitens" | input$microalgae == "hlacustris" | 
           input$microalgae == "czofingiensis" | input$microalgae == "mpusilla" |
           input$microalgae == "mendlicherianum" | input$microalgae == "smuscicola" |
           input$microalgae == "dsalina" | input$microalgae == "csubellipsoidea")
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
      } else
      {
        output$textKEGGTable <- renderText(expr = "<b>No significant KEGG 
                                           pathway enrichment detected in 
                                           the input gene set.")
      }

      if( input$microalgae == "knitens" | input$microalgae == "hlacustris" | 
          input$microalgae == "czofingiensis" | input$microalgae == "mpusilla" |
          input$microalgae == "mendlicherianum" | input$microalgae == "smuscicola" |
          input$microalgae == "dsalina"| input$microalgae == "csubellipsoidea")
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
        } else if (input$microalgae == "bprasinos")
        {
          modules.enriched.genes <- modules.enrichment.result$geneID
          for(i in 1:length(modules.enriched.genes))
          {
            modules.enriched.genes[i] <- paste(strsplit(modules.enriched.genes[i],split="/")[[1]],collapse=" ")
          }
        }
        else if (input$microalgae == "creinhardtii")
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
        } else if(input$microalgae == "knitens" | input$microalgae == "hlacustris" | 
                  input$microalgae == "czofingiensis" | input$microalgae == "mpusilla" |
                  input$microalgae == "smuscicola"| input$microalgae == "mendlicherianum" |
                  input$microalgae == "dsalina" | input$microalgae == "csubellipsoidea")
        {
          microalga.ko <- AnnotationDbi::select(org.db,columns = c("KO"),keys=keys(org.db,keytype = "GID"))
          modules.enriched.genes <- modules.enrichment.result$geneID
          for(i in 1:length(modules.enriched.genes))
          {
            current.Ks <- strsplit(modules.enriched.genes[i],split="/")[[1]]
            current.genes <- c()
            for(j in 1:length(current.Ks))
            {
              current.genes <- c(current.genes,subset(microalga.ko, KO == current.Ks[j])$GID)
            }
            modules.enriched.genes[i] <- paste(current.genes,collapse=" ")
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
        # if(input$microalgae == "otauri")
        # {
        #   gene.link.function <- ostta.gene.link
        # } else if(input$microalgae == "creinhardtii" | input$microalgae == "vcarteri"  | 
        #               input$microalgae == "cocsu" | input$microalgae == "dsalina" |
        #               input$microalgae == "mpusilla")
        # {
        #   gene.link.function <- phytozome.gene.link
        # } else if(input$microalgae == "ptricornutum")
        # {
        #   gene.link.function <- phaeodactylum.gene.link
        # } else if(input$microalgae == "ngaditana")
        # {
        #   gene.link.function <- ngaditana.gene.link
        # } else if(input$microalgae == "knitens")
        # {
        #   gene.link.function <- knitens.gene.link
        # }else if(input$microalgae == "hlacustris")
        # {
        #   gene.link.function <- ncbi.gene.link
        # }else if(input$microalgae == "zofi")
        # {
        #   gene.link.function <- zofi.gene.link
        # }
        
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
      } else 
      {
        output$text_module_kegg <- renderText(expr = "<b>No significant KEGG 
                                              module enrichment detected in the 
                                              input gene set.")
      }
    } 
    
    if((input$analysis == "go" || input$analysis == "kegg" || input$analysis == "both") && (length(target.genes) == 0))
    { 
      output$intro_conversion <- renderText(expr = "<p style=\"color:red\"><b> You forgot to input a gene set. Please, paste a gene set in the text box above or
                                    select a text file containing your gene set of interest.</b></p>")
      
      }else if ((input$analysis == "go" || input$analysis == "kegg" || input$analysis == "both") && (length(wrong.genes) == 0))
      {
        if(input$microalgae == "smuscicola" || input$microalgae == "mendlicherianum" || input$microalgae == "ngaditana" || input$microalgae == "ptricornutum")
        {
          conversion.intro.text <- paste(c("There is not data aviable to convert to common IDs the different gene names of the microalgae <b> <i>", 
                                           microalgae.names[input$microalgae],
                                           " </i> </b> "),collapse="") 
          
          output$intro_conversion <- renderText(expr = conversion.intro.text)
        }else
        {
        shinyjs::showElement(id = 'loading.conversion')
        shinyjs::hideElement(id = 'ready.conversion')
        
        ## Intro text for conversion
        conversion.intro.text <- paste(c("The table below present the results from the <b>common IDs conversion</b> 
                                      performed over the input target genes (<i>",
                                 paste(target.genes[1:3],collapse=" "),
                                 "</i> ...) from the microalgae <b> <i>", microalgae.names[input$microalgae],
                                 " </i> </b> "),collapse="") 
        
        output$intro_conversion <- renderText(expr = conversion.intro.text)
        
        if(input$microalgae == "otauri")
        {
          description_table <- read.table(file="common_ids/otauri_IDs.tsv", header = T)
          
          colnames(description_table) <- c("Genes", "Common ID or description")
          
        } else if (input$microalgae == "mpusilla")
        {
          description_table <- read.table(file="common_ids/mpusilla_IDs.tsv", header = T)
          
          colnames(description_table) <- c("Genes", "Common ID or description")
          
        } else if(input$microalgae == "bprasinos")
        {
          description_table <- read.csv(file="common_ids/bprasinos_IDs.csv", header = T, sep="\t")
          
          colnames(description_table) <- c("Genes", "Common ID or description")
          
        } else if(input$microalgae == "vcarteri")
        {
          description_table <- read.csv(file="common_ids/vcarteri_IDs.txt", header = T, sep="\t")
          
          colnames(description_table) <- c("Genes", "Common ID or description")
          
        } else if (input$microalgae == "csubellipsoidea")
        {
          description_table <- read.csv(file="common_ids/csubellipsoidea_IDs.txt", header = T, sep="\t")
          
          colnames(description_table) <- c("Genes", "Common ID or description")
          
        } else if(input$microalgae == "creinhardtii")
        {
          description_table <- read.csv(file="common_ids/creinhardtii_IDs.csv", header = F, sep="\t")
          
          colnames(description_table) <- c("Genes", "Common ID or description")
          
        } else if(input$microalgae == "dsalina")
        {
          description_table <- read.csv(file="common_ids/dsalina_IDs.txt", header = T, sep="\t")
          
          colnames(description_table) <- c("Genes", "Common ID or description")
          
        } else if(input$microalgae == "knitens")
        {
          description_table <- read.csv(file="common_ids/knitens_IDs.txt", header = T, sep="\t")
          
          colnames(description_table) <- c("Genes", "Common ID or description")
          
        } else if(input$microalgae == "hlacustris")
        {
          description_table <- read.csv(file="common_ids/hlacustris_IDs.csv", header = T, sep=",")
          
          colnames(description_table) <- c("Genes", "Common ID or description")
          
        } else if(input$microalgae == "czofingiensis")
        {
          description_table <- read.csv(file="common_ids/czofingiensis_IDs.txt", header = T, sep="\t")
          
          colnames(description_table) <- c("Genes", "Common ID or description")
        } 
        
          conversion_d <- description_table$`Common ID or description`
          names(conversion_d) <- description_table$Genes
          dataframe_conversion<-data.frame(target.genes, conversion_d[target.genes], stringsAsFactors = F)
          colnames(dataframe_conversion) <- c("Genes","Common ID or description")
          ## Add link to genes
          dataframe_conversion_links <- dataframe_conversion
          for(i in 1:length(target.genes))
          {
            dataframe_conversion_links$Genes[i] <- paste(sapply(X = target.genes[i],FUN = gene.link.function),collapse=" ")
          }
          
          
          output$output_conversion_table <- renderDataTable({
            dataframe_conversion_links
          },escape=FALSE,options =list(pageLength = 15))
        
        ## Generate UI to download go enrichment table and creating a downlodable table
        output$download_ui_for_conversion_table<- renderUI(
          tagList(downloadButton(outputId= "downloadconversionTable", "Download Common IDs Conversion Table"),tags$br(),tags$br())
        )
        
        ## Download result
        output$downloadconversionTable<- downloadHandler(
          filename= function() {
            paste("conversion_common_IDs",microalgae.names[input$microalgae] , ".tsv", sep="")
          },
          content= function(file) {
            write.table(x = dataframe_conversion,quote = F,sep = "\t",
                        file=file,row.names=FALSE,col.names=TRUE)
          })
        }#cierra el else
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
      } else if(input$microalgae == "knitens" | input$microalgae == "hlacustris" |
                input$microalgae == "czofingiensis" | input$microalgae == "mpusilla" |
                input$microalgae == "mendlicherianum" | input$microalgae == "smuscicola" |
                input$microalgae == "dsalina" | input$microalgae == "csubellipsoidea"
                    )
      {
        organism.id <- "ko"
      }
      
        output$kegg_image <- renderImage({
          pathview(gene.data = sort(genes.pathway,decreasing = TRUE),kegg.dir = "pathways",
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
    
    ## Remove previous results
    output$textTableAnnotatedGenes <- renderText(expr = "")
    output$output_gene_chip_table <- renderDataTable(expr = NULL)
    output$download_gene_chip_table <- renderUI(expr = NULL)
    output$piechart.text <- renderText(expr = "")
    output$annotation.pie.chart <- renderPlot(expr = NULL)
    output$tss.distance.text <- renderText(expr = "")
    output$distance.to.tss <- renderPlot(expr = NULL)
    output$gene.signal.text <- renderText(expr = "")
    output$annotated_genes <- renderUI(expr = NULL)
    output$tss.signal.text <- renderText(expr = "")
    output$tss_signal <- renderPlot(expr = NULL)
    shinyjs::hideElement(id = 'loading.tss.signal')
    shinyjs::hideElement(id = 'ready.tss.signal')
    shinyjs::hideElement(id = 'loading.chip')
    shinyjs::hideElement(id = 'ready.chip')
    
    shinyjs::showElement(id = 'loading.chip')
    shinyjs::hideElement(id = 'ready.chip')

    ## Load libraries
    library(ChIPseeker)
    library(ChIPpeakAnno)
    library(rtracklayer)
    library(Biostrings)
    library(seqinr)

    ## Select txdb 
    if(input$microalgae == "otauri")
    {
      gene.link.function <- ostta.gene.link
      library("org.Otauri.eg.db")
      library("TxDb.Otauri.JGI")
      txdb <- TxDb.Otauri.JGI
      org.db <- org.Otauri.eg.db
      microalgae.annotation <- read.table(file = "annotations/otauri_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/otauri_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "mpusilla")
    {
      gene.link.function <- mpusilla.gene.link
      library("org.MpusillaCCMP1545.eg.db")
      library("TxDb.MpusillaCCMP1545.Phytozome")
      org.db <- org.MpusillaCCMP1545.eg.db
      txdb <- TxDb.MpusillaCCMP1545.Phytozome
      microalgae.annotation <- read.table(file = "annotations/mpusilla_gene_annotation.tsv",sep="\t",header = T,as.is=T,comment.char = "", quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/mpusilla_gene_annotation_links.tsv",sep="\t",header = T,as.is=T,comment.char = "", quote="\"")
    } else if (input$microalgae == "bprasinos")
    {
      gene.link.function <- bathy.gene.link
      library("org.Bprasinos.eg.db")
      library("TxDb.Bprasinos.ORCAE")
      org.db <- org.Bprasinos.eg.db
      txdb <- TxDb.Bprasinos.ORCAE
      microalgae.annotation <- read.table(file = "annotations/bprasinos_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/bprasinos_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "csubellipsoidea")
    {
      gene.link.function <- csubellipsoidea.gene.link
      library("org.Csubellipsoidea.eg.db")
      library("TxDb.Csubellipsoidea.Phytozome")
      org.db <- org.Csubellipsoidea.eg.db
      txdb <- TxDb.Csubellipsoidea.Phytozome
      microalgae.annotation <- read.table(file = "annotations/csubellipsoidea_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/csubellipsoidea_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "creinhardtii")
    {
      gene.link.function <- chlamy.gene.link
      library("org.Creinhardtii.eg.db")
      library("TxDb.Creinhardtii.Phytozome")
      txdb <- TxDb.Creinhardtii.Phytozome
      org.db <- org.Creinhardtii.eg.db
      microalgae.annotation <- read.table(file = "annotations/creinhardtii_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/creinhardtii_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "vcarteri")
    {
      gene.link.function <- vcarteri.gene.link
      library("org.Vcarteri.eg.db")
      library("TxDb.Vcarteri.Phytozome")
      txdb <- TxDb.Vcarteri.Phytozome
      org.db <- org.Vcarteri.eg.db
      microalgae.annotation <- read.table(file = "annotations/vcarteri_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/vcarteri_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "dsalina")
    {
      gene.link.function <- dsalina.gene.link
      library("org.Dsalina.eg.db")
      library("TxDb.Dsalina.Phytozome")
      txdb <- TxDb.Dsalina.Phytozome
      org.db <- org.Dsalina.eg.db
      microalgae.annotation <- read.table(file = "annotations/dsalina_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/dsalina_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "hlacustris")
    {
      gene.link.function <- ncbi.gene.link
      library("org.Hlacustris.eg.db")
      library("TxDb.Hlacustris.NCBI")
      org.db <- org.Hlacustris.eg.db
      txdb <- TxDb.Hlacustris.NCBI
      microalgae.annotation <- read.table(file = "annotations/hlacustris_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/hlacustris_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "czofingiensis")
    {
      gene.link.function <- czofingiensis.gene.link
      library("org.Czofingiensis.eg.db")
      library("TxDb.Czofingiensis.Phytozome")
      org.db <- org.Czofingiensis.eg.db
      txdb <- TxDb.Czofingiensis.Phytozome
      microalgae.annotation <- read.table(file = "annotations/czofingiensis_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/czofingiensis_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "knitens")
    {
      gene.link.function <- no.gene.link
      library("org.Knitens.eg.db")
      library("TxDb.Knitens.Phycocosm")
      org.db <- org.Knitens.eg.db
      txdb <- TxDb.Knitens.Phycocosm
      microalgae.annotation <- read.table(file = "annotations/knitens_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/knitens_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "mendlicherianum")
    {
      gene.link.function <- no.gene.link
      library("org.Mendlicherianum.eg.db")
      library("TxDb.Mendlicherianum.pub")
      org.db <- org.Mendlicherianum.eg.db
      txdb <- TxDb.Mendlicherianum.pub
      microalgae.annotation <- read.table(file = "annotations/mendlicherianum_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/mendlicherianum_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "smuscicola")
    {
      gene.link.function <- no.gene.link
      library("org.Smuscicola.eg.db")
      library("TxDb.Smuscicola.pub")
      org.db <- org.Smuscicola.eg.db
      txdb <- TxDb.Smuscicola.pub
      microalgae.annotation <- read.table(file = "annotations/smuscicola_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/smuscicola_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "ptricornutum")
    {
      gene.link.function <- phaeodactylum.gene.link
      library("org.Ptricornutum.eg.db")
      library("TxDb.Ptricornutum.Ensembl.Protists")
      txdb <- TxDb.Ptricornutum.Ensembl.Protists
      org.db <- org.Ptricornutum.eg.db
      microalgae.annotation <- read.table(file = "annotations/ptricornutum_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/ptricornutum_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } else if (input$microalgae == "ngaditana")
    {
      gene.link.function <- ngaditana.gene.link
      library("org.Ngaditana.eg.db")
      library("TxDb.Ngaditana.JGI")
      txdb <- TxDb.Ngaditana.JGI
      org.db <- org.Ngaditana.eg.db
      microalgae.annotation <- read.table(file = "annotations/ngaditana_gene_annotation.tsv",sep="\t",header = T,as.is=T, quote="\"")
      microalgae.annotation.links <- read.table(file = "annotations/ngaditana_gene_annotation_links.tsv",sep="\t",header = T,as.is=T, quote="\"")
    } 
    
    
    ## Load reference genome for the chosen microalgae
    microalgae.genome.data <- read.fasta(file = paste(c("genomes/",input$microalgae,".fa"),collapse=""),
                                         seqtype = "DNA")
    microalgae.genome <- getSequence(microalgae.genome.data)
    names(microalgae.genome) <- getName(microalgae.genome.data)
    
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
    peakAnno <- annotatePeak(peak = genomic.regions, 
                             tssRegion=c(-input$promoter_length, input$promoter_length),
                             TxDb=txdb)

    ## Introductory text for pie chart
    output$piechart.text <- renderText(expr = "<b>The following piechart represents the distribution of the 
    input genomic loci overlapping different gene features. For example, if a section appears corresponding to 
    Promoter (90.75%), this means that 90.75% of the input genomic loci overlap a gene promoter:</b>")
    
    ## Plot pie chart with annotation
    output$annotation.pie.chart <- renderPlot(width = 940, height = 900, res = 120, {
        plotAnnoPie(peakAnno)
      })

    ## Introductory text for tss distance distribution
    output$tss.distance.text <- renderText(expr = "<b>The following graph represents the distance distribution 
    either upstream or downstream from genes TSS of the input genomic loci:</b>")

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
    genes.annotation.download <- subset(microalgae.annotation, Gene.ID %in% genes) 
    genes.annotation.links <- microalgae.annotation.links[genes,]  
    rownames(genes.annotation.links) <- NULL
    genes.annotation.links <- genes.annotation.links[!is.na(genes.annotation.links$Gene.ID),]

    ## Introductory text for target genes 
    annotated.genes.table.text <- "<b>The table below enumerates the potential gene targets associated with the input
    genomic loci. A gene is associated as a target of a genomic locus when it overlaps at least one of the selected
    gene features (i.e. promoter, exon, etc.). Each row represents a gene with its available annotation. Click on the 
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
    
    ## Generate UI to download gene table from genomic annotation and creating a downlodable table
    output$download_gene_chip_table<- renderUI(
      tagList(downloadButton(outputId= "downloadGeneTargets", "Download Potential Gene Targets"),tags$br(),tags$br())
    )

    ## Download result
    output$downloadGeneTargets <- downloadHandler(
      filename= function() {
        paste("annotated_genes_",microalgae.names[input$microalgae] , ".tsv", sep="")
      },
      content= function(file) {
        write.table(x = genes.annotation.download, quote=F, sep="\t",#genes.annotation.download,quote = F,sep = "\t",
                    file=file,row.names=FALSE,col.names=TRUE)
      })

    ## Extraction of the genomic features of the specified genes.
    genes.data <- subset(genes(txdb), gene_id %in% genes)

    ## Extract info from genes for profile representations
    genes.data.df <- as.data.frame(genes.data)
    exons.data <- as.data.frame(exons(txdb))
    cds.data <- as.data.frame(cds(txdb))
    
    output$annotated_genes <- renderUI({
      fluidRow(
        column (width = 3,  
                selectInput(inputId = "selected_annotated_gene", 
                            label="Choose a gene to inspect:",
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
        column(width = 9,
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
    
    observeEvent(input$individual_gene_mark, { 
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
      if(!is.null(input$bw_file))
      {
        selected.bigwig.files <- input$bw_file$data
        #selected.bed.files <- input$genomic_regions_file$data
        
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
        
        ## Colors to draw signal
        line.colors <- "blue"
        area.colors <- "lightblue"
       } #else
      # {
      #   chip.signal.mean <- rep(20,length(cord.x)) 
      #   ## Colors to draw signal
      #   line.colors <- "white"
      #   area.colors <- "white"
      # }
      
      ## Determine upper limit of the graph
      upper.lim <- 10#21
      
      ## Height to draw DNA strand
      gene.height <- -25
      cord.x <- 1:current.length
      
      ## Exon width to plot
      exon.width <- 1#2
      
      ## Cds width to plot
      cds.width <- 1.5#3
      
      ## Width of the rectangule representing the peak reagion
      peak.width <- 0.5#1
      
      ## Extract exons for target gene
      exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      ## Transform exon coordinates to current range
      min.pos <- min(exons.data.target.gene$start)

      exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$promoter_length
      exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$promoter_length

      ## Extract cds for target gene
      cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$promoter_length
      cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$promoter_length

      output$individual_gene_profile <- renderPlot({#width = 940, height = 700, res = 120, {
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

        if(nrow(cds.data.target.gene) > 0)
        {
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
        current.peaks <- as.data.frame(genomic.regions)#read.table(file=input$genomic_regions_file$data,header = F, as.is = T)
        peak.coordinates <- subset(current.peaks, seqnames == as.character(range.to.plot$seqnames) & start >= range.to.plot$start & end <= range.to.plot$end) 
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

        ## Draw profiles if bw file is provided
        if(!is.null(input$bw_file))
        {
          ## Compute base line for current TF
          current.base.line <- - 10
          
          ## Represent signal from the current TF
          lines(chip.signal.mean+current.base.line,type="l",col=line.colors,lwd=3)
          ## Determine polygon coordinates and represent it
          cord.y <- c(current.base.line,chip.signal.mean+current.base.line,current.base.line)
          cord.x <- 1:length(cord.y)
          polygon(cord.x,cord.y,col=area.colors)
        }
        
        ## Determine TFBS motifs to search for
        if(input$all.motifs)
        {
          selected.motifs.pwm <- motifs.pwm
        } else if(!is.null(input$selected.motifs))
        {
          selected.motifs.pwm <- motifs.pwm[input$selected.motifs]
        }
        
        if( input$all.motifs | !is.null(input$selected.motifs))
        {
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
                motif.pwm <- selected.motifs.pwm[[k]]
              
                hits.fw <- matchPWM(motif.pwm, peak.sequence, 
                                    min.score = paste0(input$min_score_pwm,"%"))
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
        
          # TF binding sites colors and symbol shapes
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
        }
      }) # close output$individual_gene_profile
    }) # observeEvent input$individual_gene_mark

    
    individual.gene.signal.text <- "<b> Here you can visualize an individual marked
    gene following these steps:
    <ol>
      <li>Choose or type your gene of interest from the left dropdown menu where
      autocompletion can be used.</li>
      <li>Optionally, choose or type the name of a DNA motif or transcription factor consensus
      binding sequence to be searched for in the genomic loci associated with your
      gene of interest. Alternatively tick the \"Select All Motifs\" box when you
      want to search for DNA motifs in our database. Here autocompletion is also activated.</li>
      <li>Choose a a minimun score for the selected DNA motif identification.</li>
      <li>Click on \"Go\" to visualize the selected genes with its associated input
      genomic loci, the identification of the selected DNA motifs according to the
      specified score and the signal profile over the selected gene when a BigWig
      file is provided.</li>
    </ol>
    The visualization is available on the right. The DNA sequence is represented by a
    black line. The widest blue rectangles represent exons. Untranslated regions at 5' and
    3' end of the transcript are depicted as thinner blue rectangles. Black lines also 
    represent introns. An arrow is used to mark the TSS (Transcriptional Start Site) 
    and the transcription sense. The red boxes located on top to the gene representation
    correspond to the specific input genomic loci associated with the selected gene. When
    a BigWig file is provided representing, for example, the signal of a epigenetic mark
    the signal profile over the selected gene is reprenting in lighblue. 
    </b>
    "
    
    ## Plot signal around tss and tes
    if(!is.null(input$bw_file))
    {
      # shinyjs::showElement(id = 'loading.tss.signal')
      # shinyjs::hideElement(id = 'ready.tss.signal')

      ## Extraction of the TSS 
      genes.tss <- resize(genes.data, width=1, fix='start')
      
      ## Centering around TSS with promoter length
      around.genes.tss <- genes.tss
      start(around.genes.tss) <- start(genes.tss) - input$promoter_length
      end(around.genes.tss) <- end(genes.tss) + input$promoter_length

      ## Importing bigWig file
      cvglists <- sapply(input$bw_file$data, import,
                         format="BigWig",
                         which=around.genes.tss,
                         as="RleList")

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
      
      output$tss.signal.text <- renderText(expr = "<b>The graph below represents
      the average signal around the TSS (Transcription Start Site) and TES (Transcription
      End Site) of marked genes with your input genomic loci.</b>")
      output$gene.signal.text <- renderText(expr = individual.gene.signal.text)
    } else
    {
      output$tss.signal.text <- renderText(expr = "<b> Average signal levels around
      genes TSS (Transcription Start Site) and TES (Trascription End Site) could not 
      be computed since NO BigWig file was selected.</b>")
      output$gene.signal.text <- renderText(expr = paste(individual.gene.signal.text,"<br><b> Average signal levels for 
      specific marked genes cannot be computed since NO BigWig file was selected.</b>"))
    }
    shinyjs::hideElement(id = 'loading.chip')
    shinyjs::showElement(id = 'ready.chip')
  })
  
  # Activate Funtree panel when selected
  observeEvent(input$funtree_button, {
    tree <- eventReactive(input$funtree_button,{
      
      gene.name.tree <- input$geneInt

      # Load table with orthogroups information
      ortho.table <- read.csv("Resolved_Gene_Trees/Orthogroups.tsv", header = T, sep = "\t", as.is = T,
                              fill = T, blank.lines.skip = F)
      # Find orthogroup of target gene
      found <- F
      file.name <- NULL
      i <- 1
      while (!(found) && (i <= nrow(ortho.table)))
      {
        object <- as.character(ortho.table[i,])
        gene.number <- grep(pattern = gene.name.tree, object)
        if (length(gene.number) != 0)
        {
          found <- T
          file.name <- ortho.table[i,1]
        }
        else
        {
          i <- i+1
        }
      }
      # Error if orthogroup not found
      validate(need(!is.null(file.name),"No results for this query due to not supported gene name or 
                    lack of orthologs in the selected organisms"))
      
      # Load gene tree file
      tree.name <- paste("Resolved_Gene_Trees",paste(file.name, "tree.txt", sep = "_"), sep="/")
      
      # Error if tree file not found
      validate(need(file.exists(tree.name),"No results for this query due to not supported gene name or 
                    lack of orthologs in the selected organisms"))
      
      tree <- read.tree(tree.name)
      return(tree)
    })
    
    # Tips to keep of each species with proper notation
    tips_to_keep.mp <- reactive({
      
      tree <- tree()
      # Selection of organisms
      organisms.list <- c(input$selected_organisms)
      
      # Selection of genes from the selected organism
      tips_to_keep.mp <- c()
      if ("mp" %in% organisms.list)
      {
        tips_to_keep.mp <- grep(pattern = "marchantia",tree$tip.label)
      }
      return(tips_to_keep.mp)
    })
    
    tips_to_keep.ot <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.ot <- c()
      if ("ot" %in% organisms.list)
      {
        tips_to_keep.ot <- grep(pattern = "ostreoco",tree$tip.label)
      }
      return(tips_to_keep.ot)
    })
    
    tips_to_keep.at <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.at <- c()
      if ("at" %in% organisms.list)
      {
        tips_to_keep.at <- grep(pattern = "arabidopsis",tree$tip.label)
      }
      return(tips_to_keep.at)
    })
    
    tips_to_keep.cp <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.cp <- c()
      if ("cp" %in% organisms.list)
      {
        tips_to_keep.cp <- grep(pattern = "ceratodon",tree$tip.label)
      }
      
      return(tips_to_keep.cp)
    })
    
    tips_to_keep.cr <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.cr <- c()
      if ("cr" %in% organisms.list)
      {
        tips_to_keep.cr <- grep(pattern = "chlamy",tree$tip.label)
      }
      return(tips_to_keep.cr)
    })
    
    tips_to_keep.cz <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.cz <- c()
      if ("cz" %in% organisms.list)
      {
        tips_to_keep.cz <- grep(pattern = "chromochloris",tree$tip.label)
      }
      return(tips_to_keep.cz)
    })
    
    tips_to_keep.kn <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.kn <- c()
      if ("kn" %in% organisms.list)
      {
        tips_to_keep.kn <- grep(pattern = "klebsormidium",tree$tip.label)
      }
      return(tips_to_keep.kn)
    })
    
    tips_to_keep.me <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.me <- c()
      if ("me" %in% organisms.list)
      {
        tips_to_keep.me <- grep(pattern = "mesotaenium",tree$tip.label)
      }
      return(tips_to_keep.me)
    })
    
    tips_to_keep.mi <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.mi <- c()
      if ("mi" %in% organisms.list)
      {
        tips_to_keep.mi <- grep(pattern = "micromonas",tree$tip.label)
      }
      return(tips_to_keep.mi)
    })
    
    tips_to_keep.pp <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.pp <- c()
      if ("pp" %in% organisms.list)
      {
        tips_to_keep.pp <- grep(pattern = "physcomitrium",tree$tip.label)
      }
      return(tips_to_keep.pp)
    })
    
    tips_to_keep.sl <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.sl <- c()
      if ("sl" %in% organisms.list)
      {
        tips_to_keep.sl <- grep(pattern = "solanum",tree$tip.label)
      }
      return(tips_to_keep.sl)
    })
    
    tips_to_keep.sm <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.sm <- c()
      if ("sm" %in% organisms.list)
      {
        tips_to_keep.sm <- grep(pattern = "selaginella",tree$tip.label)
      }
      return(tips_to_keep.sm)
    })
    
    tips_to_keep.sp <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.sp <- c()
      if ("sp" %in% organisms.list)
      {
        tips_to_keep.sp <- grep(pattern = "spirogloea",tree$tip.label)
      }
      return(tips_to_keep.sp)
    })
    
    tips_to_keep.ta <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.ta <- c()
      if ("ta" %in% organisms.list)
      {
        tips_to_keep.ta <- grep(pattern = "triticum",tree$tip.label)
      }
      return(tips_to_keep.ta)
    })
    
    tips_to_keep.vc <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.vc <- c()
      if ("vc" %in% organisms.list)
      {
        tips_to_keep.vc <- grep(pattern = "volvox",tree$tip.label)
      }
      return(tips_to_keep.vc)
    })
    
    tips_to_keep.bp <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.bp <- c()
      if ("bp" %in% organisms.list)
      {
        tips_to_keep.bp <- grep(pattern = "bathycoccus",tree$tip.label)
      }
      return(tips_to_keep.bp)
    })
    
    tips_to_keep.cri <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.cri <- c()
      if ("cri" %in% organisms.list)
      {
        tips_to_keep.cri <- grep(pattern = "ceratopteris",tree$tip.label)
      }
      return(tips_to_keep.cri)
    })
    
    tips_to_keep.ds <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.ds <- c()
      if ("ds" %in% organisms.list)
      {
        tips_to_keep.ds <- grep(pattern = "dunaliella",tree$tip.label)
      }
      return(tips_to_keep.ds)
    })
    
    tips_to_keep.os <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.os <- c()
      if ("os" %in% organisms.list)
      {
        tips_to_keep.os <- grep(pattern = "oryza",tree$tip.label)
      }
      return(tips_to_keep.os)
    })
    
    tips_to_keep.smag <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.smag <- c()
      if ("smag" %in% organisms.list)
      {
        tips_to_keep.smag <- grep(pattern = "sphagum",tree$tip.label)
      }
      return(tips_to_keep.smag)
    })
    
    tips_to_keep.tp <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.tp <- c()
      if ("tp" %in% organisms.list)
      {
        tips_to_keep.tp <- grep(pattern = "thuja",tree$tip.label)
      }
      return(tips_to_keep.tp)
    })
    
    tips_to_keep.aa <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.aa <- c()
      if ("aa" %in% organisms.list)
      {
        tips_to_keep.aa <- grep(pattern = "anthoceros",tree$tip.label)
      }
      return(tips_to_keep.aa)
    })
    
    tips_to_keep.um <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.um <- c()
      if ("um" %in% organisms.list)
      {
        tips_to_keep.um <- grep(pattern = "ulva",tree$tip.label)
      }
      return(tips_to_keep.um)
    })
    
    tips_to_keep.rs <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.rs <- c()
      if ("rs" %in% organisms.list)
      {
        tips_to_keep.rs <- grep(pattern = "raphidocelis",tree$tip.label)
      }
      return(tips_to_keep.rs)
    })
    
    tips_to_keep.cyc <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.cyc <- c()
      if ("cyc" %in% organisms.list)
      {
        tips_to_keep.cyc <- grep(pattern = "cycas",tree$tip.label)
      }
      return(tips_to_keep.cyc)
    })
    
    tips_to_keep.pu <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.pu <- c()
      if ("pu" %in% organisms.list)
      {
        tips_to_keep.pu <- grep(pattern = "porphyra",tree$tip.label)
      }
      return(tips_to_keep.pu)
    })
    
    tips_to_keep.pt <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.pt <- c()
      if ("pt" %in% organisms.list)
      {
        tips_to_keep.pt <- grep(pattern = "phaeodactylum",tree$tip.label)
      }
      return(tips_to_keep.pt)
    })
    
    tips_to_keep.ng <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.ng <- c()
      if ("ng" %in% organisms.list)
      {
        tips_to_keep.ng <- grep(pattern = "gaditana",tree$tip.label)
      }
      return(tips_to_keep.ng)
    })
    
    tips_to_keep.cyano <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.cyano <- c()
      if ("cyano" %in% organisms.list)
      {
        tips_to_keep.cyano <- grep(pattern = "cyanophora",tree$tip.label)
      }
      return(tips_to_keep.cyano)
    })
  
    tips_to_keep.ca <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.ca <- c()
      if ("ca" %in% organisms.list)
      {
        tips_to_keep.ca <- grep(pattern = "chlorokybus",tree$tip.label)
      }
      return(tips_to_keep.ca)
    })
    
    tips_to_keep.mv <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.mv <- c()
      if ("mv" %in% organisms.list)
      {
        tips_to_keep.mv <- grep(pattern = "mesostigma",tree$tip.label)
      }
      return(tips_to_keep.mv)
    })
    
    tips_to_keep.af <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.af <- c()
      if ("af" %in% organisms.list)
      {
        tips_to_keep.af <- grep(pattern = "azolla",tree$tip.label)
      }
      return(tips_to_keep.af)
    })
    
    tips_to_keep.sc <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.sc <- c()
      if ("sc" %in% organisms.list)
      {
        tips_to_keep.sc <- grep(pattern = "salvinia",tree$tip.label)
      }
      return(tips_to_keep.sc)
    })
    
    tips_to_keep.aegi <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.aegi <- c()
      if ("aegi" %in% organisms.list)
      {
        tips_to_keep.aegi <- grep(pattern = "aegilops",tree$tip.label)
      }
      return(tips_to_keep.aegi)
    })
    
    tips_to_keep.sb <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.sb <- c()
      if ("sb" %in% organisms.list)
      {
        tips_to_keep.sb <- grep(pattern = "sorghum",tree$tip.label)
      }
      return(tips_to_keep.sb)
    })
    
    tips_to_keep.zm <- reactive({
      
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      tips_to_keep.zm <- c()
      if ("zm" %in% organisms.list)
      {
        tips_to_keep.zm <- grep(pattern = "mays",tree$tip.label)
      }
      return(tips_to_keep.zm)
    })
    
 
    
  # Create complete gene tree with the proper name for each gene
    tree_adj <- reactive({
      tree <- tree()
      organisms.list <- c(input$selected_organisms)
      
      if ("mp" %in% organisms.list)
      {
        tips_to_keep.mp <- grep(pattern = "marchantia",tree$tip.label) 
        if (length(tips_to_keep.mp) != 0)
        {
          mp.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.mp]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.mp] <- mp.v
        }
      }
      
      if ("ot" %in% organisms.list)
      {
        tips_to_keep.ot <- grep(pattern = "ostreoco",tree$tip.label) 
        if (length(tips_to_keep.ot) != 0)
        {
          ost.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.ot]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.ot] <- ost.v
        }
      }
      
      if ("at" %in% organisms.list)
      {
        tips_to_keep.at <- grep(pattern = "arabidopsis",tree$tip.label) 
        if (length(tips_to_keep.at) != 0)
        {
          arabi.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.at]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.at] <- arabi.v
        }
      }
      
      if ("cp" %in% organisms.list)
      {
        tips_to_keep.cp <- grep(pattern = "ceratodon",tree$tip.label) 
        if (length(tips_to_keep.cp) != 0)
        {
          cer.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.cp]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.cp] <- cer.v
        }
      }
      
      if ("cr" %in% organisms.list)
      {
        tips_to_keep.cr <- grep(pattern = "chlamy",tree$tip.label)
        if (length(tips_to_keep.cr) != 0)
        {
          chlamy.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.cr]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.cr] <- chlamy.v
        }
      }
      
      if ("cz" %in% organisms.list)
      {
        tips_to_keep.cz <- grep(pattern = "chromochloris",tree$tip.label) 
        if (length(tips_to_keep.cz) != 0)
        {
          chromo.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.cz]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.cz] <- chromo.v
        }
      }
      
      if ("kn" %in% organisms.list)
      {
        tips_to_keep.kn <- grep(pattern = "klebsormidium",tree$tip.label) 
        if (length(tips_to_keep.kn) != 0)
        {
          klebs.v1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.kn]), "_"), function(x) x[[6]])
          klebs.v2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.kn]), "_"), function(x) x[[7]])
          klebs.v <- paste(klebs.v1, klebs.v2, sep = "_")
          tree$tip.label[tips_to_keep.kn] <- klebs.v
        }
      }
      
      if ("me" %in% organisms.list)
      {
        tips_to_keep.me <- grep(pattern = "mesotaenium",tree$tip.label) 
        if (length(tips_to_keep.me) != 0)
        {
          meso.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.me]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.me] <- meso.v
        }
      }
      
      if ("mi" %in% organisms.list)
      {
        tips_to_keep.mi <- grep(pattern = "micromonas",tree$tip.label) 
        if (length(tips_to_keep.mi) != 0)
        {
          micro.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.mi]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.mi] <- micro.v
        }
      }
      
      if ("pp" %in% organisms.list)
      {
        tips_to_keep.pp <- grep(pattern = "physcomitrium",tree$tip.label) 
        if (length(tips_to_keep.pp) != 0)
        {
          phys.v1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.pp]), "_"), function(x) x[[6]])
          phys.v2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.pp]), "_"), function(x) x[[7]])
          phys.v <- paste(phys.v1, phys.v2, sep = "_")
          tree$tip.label[tips_to_keep.pp] <- phys.v
        }
      }
      
      if ("sl" %in% organisms.list)
      {
        tips_to_keep.sl <- grep(pattern = "solanum",tree$tip.label) 
        if (length(tips_to_keep.sl) != 0)
        {
          sola.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.sl]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.sl] <- sola.v
        }
      }
      
      if ("sm" %in% organisms.list)
      {
        tips_to_keep.sm <- grep(pattern = "selaginella",tree$tip.label) 
        if (length(tips_to_keep.sm) != 0)
        {
          sel.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.sm]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.sm] <- sel.v
        }
      }
      
      if ("sp" %in% organisms.list)
      {
        tips_to_keep.sp <- grep(pattern = "spirogloea",tree$tip.label) 
        if (length(tips_to_keep.sp) != 0)
        {
          spiro.v <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.sp]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.sp] <- spiro.v
        }
      }
      
      if ("ta" %in% organisms.list)
      {
        tips_to_keep.ta <- grep(pattern = "triticum",tree$tip.label) 
        if (length(tips_to_keep.ta) != 0)
        {
          tri.v1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.ta]), "_"), function(x) x[[6]])
          tri.v2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.ta]), "_"), function(x) x[[7]])
          tri.v3 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.ta]), "_"), function(x) x[[8]])
          tri.v <- paste(tri.v1,tri.v2,tri.v3, sep = "_")
          tree$tip.label[tips_to_keep.ta] <- tri.v
        }
      }
      
      if ("vc" %in% organisms.list)
      {
        tips_to_keep.vc <- grep(pattern = "volvox",tree$tip.label)
        if (length(tips_to_keep.vc) != 0)
        {
          vc.vec <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.vc]), "_"), function(x) x[[6]])
          vc1 <- sapply(strsplit(as.character(vc.vec), "\\."), function(x) x[[1]])
          vc2 <- sapply(strsplit(as.character(vc.vec), "\\."), function(x) x[[2]])
          vc.v <- paste(vc1, vc2, sep = ".")
          tree$tip.label[tips_to_keep.vc] <- vc.v
        }
      }
      
      if ("bp" %in% organisms.list)
      {
        tips_to_keep.bp <- grep(pattern = "bathycoccus",tree$tip.label)
        if (length(tips_to_keep.bp) != 0)
        {
          bp.vec <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.bp]), "_"), function(x) x[[6]])
          bp.v <- sapply(strsplit(as.character(bp.vec), "leng"), function(x) x[[1]])
          tree$tip.label[tips_to_keep.bp] <- bp.v
        }
      }
      
      if ("cri" %in% organisms.list)
      {
        tips_to_keep.cri <- grep(pattern = "ceratopteris",tree$tip.label)
        if (length(tips_to_keep.cri) != 0)
        {
          cri.vec <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.cri]), "_"), function(x) x[[6]])
          cri1 <- sapply(strsplit(as.character(cri.vec), "\\."), function(x) x[[1]])
          cri2 <- sapply(strsplit(as.character(cri.vec), "\\."), function(x) x[[2]])
          cri.v <- paste(cri1, cri2, sep = ".")
          tree$tip.label[tips_to_keep.cri] <- cri.v
        }
      }
      
      if ("ds" %in% organisms.list)
      {
        tips_to_keep.ds <- grep(pattern = "dunaliella",tree$tip.label)
        if (length(tips_to_keep.ds) != 0)
        {
          ds.vec <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.ds]), "_"), function(x) x[[6]])
          ds1 <- sapply(strsplit(as.character(ds.vec), "\\."), function(x) x[[1]])
          ds2 <- sapply(strsplit(as.character(ds.vec), "\\."), function(x) x[[2]])
          ds.v <- paste(ds1, ds2, sep = ".")
          tree$tip.label[tips_to_keep.ds] <- ds.v
        }
      }
      
      if ("os" %in% organisms.list)
      {
        tips_to_keep.os <- grep(pattern = "oryza",tree$tip.label)
        if (length(tips_to_keep.os) != 0)
        {
          os.vec <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.os]), "_"), function(x) x[[7]])
          os.v <- sapply(strsplit(as.character(os.vec), "\\."), function(x) x[[1]])
          tree$tip.label[tips_to_keep.os] <- os.v
        }
      }
      
      if ("smag" %in% organisms.list)
      {
        tips_to_keep.smag <- grep(pattern = "sphagum",tree$tip.label)
        if (length(tips_to_keep.smag) != 0)
        {
          smag.vec <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.smag]), "_"), function(x) x[[6]])
          smag.v <- sapply(strsplit(as.character(smag.vec), "\\."), function(x) x[[1]])
          tree$tip.label[tips_to_keep.smag] <- smag.v
        }
      }
      
      if ("tp" %in% organisms.list)
      {
        tips_to_keep.tp <- grep(pattern = "thuja",tree$tip.label)
        if (length(tips_to_keep.tp) != 0)
        {
          tp.vec <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.tp]), "_"), function(x) x[[6]])
          tp1 <- sapply(strsplit(as.character(tp.vec), "\\."), function(x) x[[1]])
          tp2 <- sapply(strsplit(as.character(tp.vec), "\\."), function(x) x[[2]])
          tp.v <- paste(tp1, tp2, sep = ".")
          tree$tip.label[tips_to_keep.tp] <- tp.v
        }
      }
      
      if ("aa" %in% organisms.list)
      {
        tips_to_keep.aa <- grep(pattern = "anthoceros",tree$tip.label)
        if (length(tips_to_keep.aa) != 0)
        {
          aa.vec <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.aa]), "_"), function(x) x[[7]])
          aa1 <- sapply(strsplit(as.character(aa.vec), "\\."), function(x) x[[3]])
          aa2 <- sapply(strsplit(as.character(aa.vec), "\\."), function(x) x[[4]])
          aa.v <- paste(aa1, aa2, sep = ".")
          tree$tip.label[tips_to_keep.aa] <- aa.v
        }
      }
      
      if ("um" %in% organisms.list)
      {
        tips_to_keep.um <- grep(pattern = "ulva",tree$tip.label)
        if (length(tips_to_keep.um) != 0)
        {
          um.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.um]), "_"), function(x) x[[6]])
          um.vec2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.um]), "_"), function(x) x[[7]])
          um2 <- sapply(strsplit(as.character(um.vec2), "\\."), function(x) x[[1]])
          um.v <- paste(um.vec1, um2, sep = "_")
          tree$tip.label[tips_to_keep.um] <- um.v
        }
      }
      
      if ("rs" %in% organisms.list)
      {
        tips_to_keep.rs <- grep(pattern = "raphidocelis",tree$tip.label)
        if (length(tips_to_keep.rs) != 0)
        {
          rs.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.rs]), "_"), function(x) x[[6]])
          rs.vec2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.rs]), "_"), function(x) x[[7]])
          rs.v <- paste(rs.vec1, rs.vec2, sep = "_")
          tree$tip.label[tips_to_keep.rs] <- rs.v
        }
      }
      
      if ("cyc" %in% organisms.list)
      {
        tips_to_keep.cyc <- grep(pattern = "cycas",tree$tip.label)
        if (length(tips_to_keep.cyc) != 0)
        {
          cyc.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.cyc]), "_"), function(x) x[[6]])
          cyc.vec2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.cyc]), "_"), function(x) x[[7]])
          cyc.v <- paste(cyc.vec1, cyc.vec2, sep = "_")
          tree$tip.label[tips_to_keep.cyc] <- cyc.v
        }
      }
      
      if ("pu" %in% organisms.list)
      {
        tips_to_keep.pu <- grep(pattern = "porphyra",tree$tip.label)
        if (length(tips_to_keep.pu) != 0)
        {
          pu.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.pu]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.pu] <- pu.vec1
        }
      }
      
      if ("pt" %in% organisms.list)
      {
        tips_to_keep.pt <- grep(pattern = "phaeodactylum",tree$tip.label)
        if (length(tips_to_keep.pt) != 0)
        {
          pt.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.pt]), "_"), function(x) x[[6]])
          pt.vec2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.pt]), "_"), function(x) x[[7]])
          pt.v <- paste(pt.vec1, pt.vec2, sep = "_")
          tree$tip.label[tips_to_keep.pt] <- pt.v
        }
      }
      
      if ("ng" %in% organisms.list)
      {
        tips_to_keep.ng <- grep(pattern = "gaditana",tree$tip.label)
        if (length(tips_to_keep.ng) != 0)
        {
          ng.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.ng]), "_"), function(x) x[[6]])
          ng.vec2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.ng]), "_"), function(x) x[[7]])
          ng.v <- paste(ng.vec1, ng.vec2, sep = "_")
          tree$tip.label[tips_to_keep.ng] <- ng.v
        }
      }
      
      if ("cyano" %in% organisms.list)
      {
        tips_to_keep.cyano <- grep(pattern = "cyanophora",tree$tip.label)
        if (length(tips_to_keep.cyano) != 0)
        {
          cyano.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.cyano]), "_"), function(x) x[[6]])
          cyano.vec2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.cyano]), "_"), function(x) x[[7]])
          cyano.v <- paste(cyano.vec1, cyano.vec2, sep = "_")
          tree$tip.label[tips_to_keep.cyano] <- cyano.v
        }
      }
      
      if ("ca" %in% organisms.list)
      {
        tips_to_keep.ca <- grep(pattern = "chlorokybus",tree$tip.label)
        if (length(tips_to_keep.ca) != 0)
        {
          ca.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.ca]), "_"), function(x) x[[6]])
          ca.vec2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.ca]), "_"), function(x) x[[7]])
          ca.v <- paste(ca.vec1, ca.vec2, sep = "_")
          tree$tip.label[tips_to_keep.ca] <- ca.v
        }
      }
      
      if ("mv" %in% organisms.list)
      {
        tips_to_keep.mv <- grep(pattern = "mesostigma",tree$tip.label)
        if (length(tips_to_keep.mv) != 0)
        {
          mv.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.mv]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.mv] <- mv.vec1
        }
      }
      
      if ("af" %in% organisms.list)
      {
        tips_to_keep.af <- grep(pattern = "azolla",tree$tip.label)
        if (length(tips_to_keep.af) != 0)
        {
          af.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.af]), "_"), function(x) x[[6]])
          af.vec2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.af]), "_"), function(x) x[[7]])
          af.v <- paste(af.vec1, af.vec2, sep = "_")
          tree$tip.label[tips_to_keep.af] <- af.v
        }
      }
      
      if ("sc" %in% organisms.list)
      {
        tips_to_keep.sc <- grep(pattern = "salvinia",tree$tip.label)
        if (length(tips_to_keep.sc) != 0)
        {
          sc.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.sc]), "_"), function(x) x[[6]])
          sc.vec2 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.sc]), "_"), function(x) x[[8]])
          sc.v <- paste(sc.vec1, sc.vec2, sep = "_")
          tree$tip.label[tips_to_keep.sc] <- sc.v
        }
      }
      
      if ("aegi" %in% organisms.list)
      {
        tips_to_keep.aegi <- grep(pattern = "aegilops",tree$tip.label)
        if (length(tips_to_keep.aegi) != 0)
        {
          aegi.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.aegi]), "_"), function(x) x[[6]])
          aegi.v <- sapply(strsplit(as.character(aegi.vec1), "[.]"), function(x) x[[1]])
          tree$tip.label[tips_to_keep.aegi] <- aegi.v
        }
      }
      
      if ("sb" %in% organisms.list)
      {
        tips_to_keep.sb <- grep(pattern = "sorghum",tree$tip.label)
        if (length(tips_to_keep.sb) != 0)
        {
          sb.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.sb]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.sb] <- sb.vec1
        }
      }
      
      if ("zm" %in% organisms.list)
      {
        tips_to_keep.zm <- grep(pattern = "mays",tree$tip.label)
        if (length(tips_to_keep.zm) != 0)
        {
          zm.vec1 <- sapply(strsplit(as.character(tree$tip.label[tips_to_keep.zm]), "_"), function(x) x[[6]])
          tree$tip.label[tips_to_keep.zm] <- zm.vec1
        }
      }
      
      return(tree)
    })
    # Generate reduced tree when the corresponding button is activated
    tree_reduced <- eventReactive(input$funtree_button,{
      
      tree <- tree_adj()
      # Define tips to keep (selected organisms) and generate the reduced tree
      tips_to_keep.global <- c(tips_to_keep.mp(), tips_to_keep.ot(), tips_to_keep.at(), tips_to_keep.cp(),
                               tips_to_keep.cr(), tips_to_keep.cz(), tips_to_keep.kn(), tips_to_keep.me(),
                               tips_to_keep.mi(), tips_to_keep.pp(), tips_to_keep.sl(), tips_to_keep.sm(),
                               tips_to_keep.sp(), tips_to_keep.ta(), tips_to_keep.vc(), tips_to_keep.bp(),
                               tips_to_keep.cri(), tips_to_keep.ds(), tips_to_keep.os(), tips_to_keep.smag(),
                               tips_to_keep.tp(), tips_to_keep.aa(), tips_to_keep.um(), tips_to_keep.rs(),
                               tips_to_keep.cyc(), tips_to_keep.pu(), tips_to_keep.pt(), tips_to_keep.ng(),
                               tips_to_keep.cyano(), tips_to_keep.ca(), tips_to_keep.mv(), tips_to_keep.af(),
                               tips_to_keep.sc(), tips_to_keep.aegi(), tips_to_keep.sb(), tips_to_keep.zm())
      validate(need(length(tips_to_keep.global) > 1,"Unable to construct tree with a single tip, please select more
                    organisms"))
     
      tips_to_drop <- setdiff(1:length(tree$tip.label), tips_to_keep.global)
      tree_reduced <- drop.tip(tree, tips_to_drop)
      return(tree_reduced)
    })
    
    output$treeTips <- renderPrint({
      #validate(need(!is.null(tree_reduced()$tip.label),"No results for this query due to not supported gene name or 
      #              lack of homologs in the selected organisms"))
      cat(tree_reduced()$tip.label)
    })
    # Create UI for download buttons
    output$download_tips<- renderUI(
      tagList(downloadButton(outputId= "downloadTips", "Download Tree Tips"))
    )
    
    ## Download result
    output$downloadTips <- downloadHandler(
      filename= function() {
        paste("tips", ".tsv", sep="")
      },
      content= function(file) {
        write.table(x = tree_reduced()$tip.label, quote=F, sep="\t",
                    file=file,row.names=FALSE,col.names=FALSE)
      })
    
    tree_plot <-  eventReactive(input$funtree_button, {
      
      # Define previous variables
      tree_reduced <- tree_reduced()
      gene.name.tree <- input$geneInt
      tree <- tree_adj()
      
      tips_to_keep.mp <- tips_to_keep.mp()
      tips_to_keep.ot <- tips_to_keep.ot()
      tips_to_keep.at <- tips_to_keep.at()
      tips_to_keep.cp <- tips_to_keep.cp()
      tips_to_keep.cr <- tips_to_keep.cr()
      tips_to_keep.cz <- tips_to_keep.cz()
      tips_to_keep.kn <- tips_to_keep.kn()
      tips_to_keep.me <- tips_to_keep.me()
      tips_to_keep.mi <- tips_to_keep.mi()
      tips_to_keep.pp <- tips_to_keep.pp()
      tips_to_keep.sl <- tips_to_keep.sl()
      tips_to_keep.sm <- tips_to_keep.sm()
      tips_to_keep.sp <- tips_to_keep.sp()
      tips_to_keep.ta <- tips_to_keep.ta()
      tips_to_keep.vc <- tips_to_keep.vc()
      tips_to_keep.bp <- tips_to_keep.bp()
      tips_to_keep.cri <- tips_to_keep.cri()
      tips_to_keep.ds <- tips_to_keep.ds()
      tips_to_keep.os <- tips_to_keep.os()
      tips_to_keep.smag <- tips_to_keep.smag()
      tips_to_keep.tp <- tips_to_keep.tp()
      tips_to_keep.aa <- tips_to_keep.aa()
      tips_to_keep.um <- tips_to_keep.um()
      tips_to_keep.rs <- tips_to_keep.rs()
      tips_to_keep.cyc <- tips_to_keep.cyc()
      tips_to_keep.pu <- tips_to_keep.pu()
      tips_to_keep.pt <- tips_to_keep.pt()
      tips_to_keep.ng <- tips_to_keep.ng()
      tips_to_keep.cyano <- tips_to_keep.cyano()
      tips_to_keep.ca <- tips_to_keep.ca()
      tips_to_keep.mv <- tips_to_keep.mv()
      tips_to_keep.af <- tips_to_keep.af()
      tips_to_keep.sc <- tips_to_keep.sc()
      tips_to_keep.aegi <- tips_to_keep.aegi()
      tips_to_keep.sb <- tips_to_keep.sb()
      tips_to_keep.zm <- tips_to_keep.zm()
      
      if (length(tree_reduced$tip.label) < 2)
      {
        cat("")
      }
      else 
      {
        # Highlight the target gene
        high.gene <<- tree_reduced$tip.label[grep(pattern = gene.name.tree, tree_reduced$tip.label)]
        
        # Color asignment per species
        col.factor <- c()
        org.factor <- c()
        
        for (i in 1:length(tree_reduced$tip.label))
        {
          if (tree_reduced$tip.label[i] %in% high.gene)
          {
            col.factor <- c(col.factor,"#CD0000")
            org.factor <- c(org.factor,"Gen of interest")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.mp])
          {
            col.factor <- c(col.factor,"#006400")
            org.factor <- c(org.factor,"Marchantia")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.ot])
          {
            col.factor <- c(col.factor,"#00008B")
            org.factor <- c(org.factor,"Ostreococcus")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.at])
          {
            col.factor <- c(col.factor,"#CD661D")
            org.factor <- c(org.factor,"Arabidopsis")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.cp])
          {
            col.factor <- c(col.factor,"#458B74")
            org.factor <- c(org.factor,"Ceratodon")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.cr])
          {
            col.factor <- c(col.factor,"#8B7355")
            org.factor <- c(org.factor,"Chlamydomonas")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.cz])
          {
            col.factor <- c(col.factor,"#458B00")
            org.factor <- c(org.factor,"Chromochloris")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.kn])
          {
            col.factor <- c(col.factor,"#CD1076")
            org.factor <- c(org.factor,"Klebsormidium")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.me])
          {
            col.factor <- c(col.factor,"#8B8878")
            org.factor <- c(org.factor,"Mesotaenium")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.mi])
          {
            col.factor <- c(col.factor,"#666666")
            org.factor <- c(org.factor,"Micromonas")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.pp])
          {
            col.factor <- c(col.factor,"#B8860B")
            org.factor <- c(org.factor,"Physcomitrium")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.sl])
          {
            col.factor <- c(col.factor,"#8B008B")
            org.factor <- c(org.factor,"Solanum")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.sm])
          {
            col.factor <- c(col.factor,"#6E8B3D")
            org.factor <- c(org.factor,"Selaginella")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.sp])
          {
            col.factor <- c(col.factor,"#79CDCD")
            org.factor <- c(org.factor,"Spirogloea")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.ta])
          {
            col.factor <- c(col.factor,"#CDCD00")
            org.factor <- c(org.factor,"Triticum")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.vc])
          {
            col.factor <- c(col.factor,"#16317d")
            org.factor <- c(org.factor,"Volvox")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.bp])
          {
            col.factor <- c(col.factor,"#007e2f")
            org.factor <- c(org.factor,"Bathycoccus")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.cri])
          {
            col.factor <- c(col.factor,"#ffcd12")
            org.factor <- c(org.factor,"Ceratopteris")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.ds])
          {
            col.factor <- c(col.factor,"#b86092")
            org.factor <- c(org.factor,"Dunaliella")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.os])
          {
            col.factor <- c(col.factor,"#721b3e")
            org.factor <- c(org.factor,"Oryza")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.smag])
          {
            col.factor <- c(col.factor,"#00b7a7")
            org.factor <- c(org.factor,"Sphangum")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.tp])
          {
            col.factor <- c(col.factor,"#67000d")
            org.factor <- c(org.factor,"Thuja")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.aa])
          {
            col.factor <- c(col.factor,"#5b2c6f")
            org.factor <- c(org.factor,"Anthoceros")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.um])
          {
            col.factor <- c(col.factor,"#15e71b")
            org.factor <- c(org.factor,"Ulva")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.rs])
          {
            col.factor <- c(col.factor,"#e67e22")
            org.factor <- c(org.factor,"Raphidocelis")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.cyc])
          {
            col.factor <- c(col.factor,"#873600")
            org.factor <- c(org.factor,"Cycas")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.pu])
          {
            col.factor <- c(col.factor,"#dc1c0f")
            org.factor <- c(org.factor,"Porphyra")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.pt])
          {
            col.factor <- c(col.factor,"#a04000")
            org.factor <- c(org.factor,"Phaeodactylum")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.ng])
          {
            col.factor <- c(col.factor,"#935116")
            org.factor <- c(org.factor,"Nannochloropsis")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.cyano])
          {
            col.factor <- c(col.factor,"#2874a6")
            org.factor <- c(org.factor,"Cyanophora")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.ca])
          {
            col.factor <- c(col.factor,"#0b5345")
            org.factor <- c(org.factor,"Chlorokybus")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.mv])
          {
            col.factor <- c(col.factor,"#283747")
            org.factor <- c(org.factor,"Mesostigma")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.af])
          {
            col.factor <- c(col.factor,"#145a32")
            org.factor <- c(org.factor,"Azolla")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.sc])
          {
            col.factor <- c(col.factor,"#3339e6")
            org.factor <- c(org.factor,"Salvinia")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.aegi])
          {
            col.factor <- c(col.factor,"#e6338f")
            org.factor <- c(org.factor,"Aegilops")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.sb])
          {
            col.factor <- c(col.factor,"#cd016a")
            org.factor <- c(org.factor,"Sorghum")
          }
          else if (tree_reduced$tip.label[i] %in% tree$tip.label[tips_to_keep.zm])
          {
            col.factor <- c(col.factor,"#666909")
            org.factor <- c(org.factor,"Zea")
          }
          
        }
        
        #Matrix with labels and colors and transform to dplyr format
        data.tree <- data.frame(node = 1:length(tree_reduced$tip.label), label = tree_reduced$tip.label,
                                col = col.factor, org = org.factor)
        
        d2 <- dplyr::mutate(data.tree, lab = data.tree$label,
                            color = data.tree$col,
                            organism = data.tree$org,
                            name = glue("<i style='color:{color}'> {lab} </i>"))
        
        
        tree_plot <- ggtree(tree_reduced) %<+% d2 + geom_tiplab() + theme(legend.position =) +
          xlim(0, max(tree_reduced$edge.length)*3) + geom_tiplab(aes(label = label, color = organism)) +
          scale_color_manual(values = unique(d2$col), breaks = unique(d2$org)) +
          geom_highlight(mapping=aes(subset = label %in% high.gene,
                                     node = node,
                                     fill = as.factor(node)), extend = 0.8) + 
          labs(fill = "Node of interest")
      
      }})
    # Create tree output
    output$download_tree<- renderUI(
      tagList(downloadButton(outputId= "downloadTree", "Download Tree"))
    )

    image_height <- function(){300 + 11*length(tree_reduced()$tip.label)}
    image_width <- function(){200 + 400*max(tree_reduced()$edge.length)}
    ## Download result
    output$downloadTree <- downloadHandler(
      filename= function() {
        paste("tree", ".png", sep="")
      },
      content= function(file) {
        png(file, height = image_height(), width = image_width())
        plot(tree_plot())
        dev.off()
      })

    output$treePlot <- renderPlot({
      plot(tree_plot())
    }, height = image_height(),
    width = image_height()
    )
  
  }
  )
  
  
})

# Run the application 
shinyApp(ui = ui, server = server)