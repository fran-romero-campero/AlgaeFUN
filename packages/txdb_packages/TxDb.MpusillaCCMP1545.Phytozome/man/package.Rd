\name{TxDb.MpusillaCCMP1545.Phytozome}
\docType{package}

\alias{TxDb.MpusillaCCMP1545.Phytozome-package}
\alias{TxDb.MpusillaCCMP1545.Phytozome}
\alias{TxDb.MpusillaCCMP1545.Phytozome}


\title{Annotation package for TxDb object(s)}

\description{  
  This package loads one or more TxDb objects.  Such TxDb
  objects are an R interface to prefabricated databases contained by
  this package.

  The names of any objects exposed by this package indicate the origin and
  resources exposed.  So for example TxDb.Hsapiens.UCSC.hg19.knownGene
  would be a TxDb object, of Homo sapiens data from UCSC build
  hg19 based on the knownGene Track.
}

\note{
  This data package was made from resources at Phytozome on
  2021-02-16 06:49:23 +0100 (Tue, 16 Feb 2021) and based on the Phytozome
}

\author{Christina Arvanitidou, Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.MpusillaCCMP1545.Phytozome)
## list the contents that are loaded into memory
ls('package:TxDb.MpusillaCCMP1545.Phytozome')
## show the db object that is loaded by calling it's name
TxDb.MpusillaCCMP1545.Phytozome

}

\keyword{package}
\keyword{data}
