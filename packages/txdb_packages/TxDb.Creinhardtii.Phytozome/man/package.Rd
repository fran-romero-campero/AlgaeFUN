\name{TxDb.Creinhardtii.Phytozome}
\docType{package}

\alias{TxDb.Creinhardtii.Phytozome-package}
\alias{TxDb.Creinhardtii.Phytozome}
\alias{TxDb.Creinhardtii.Phytozome}


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
  2018-12-14 12:00:49 +0100 (Fri, 14 Dec 2018) and based on the Phytozome
}

\author{Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Creinhardtii.Phytozome)
## list the contents that are loaded into memory
ls('package:TxDb.Creinhardtii.Phytozome')
## show the db object that is loaded by calling it's name
TxDb.Creinhardtii.Phytozome

}

\keyword{package}
\keyword{data}
