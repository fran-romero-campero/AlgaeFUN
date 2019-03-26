\name{TxDb.Otauri.JGI}
\docType{package}

\alias{TxDb.Otauri.JGI-package}
\alias{TxDb.Otauri.JGI}
\alias{TxDb.Otauri.JGI}


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
  This data package was made from resources at JGI on
  2019-01-20 12:45:44 +0100 (Sun, 20 Jan 2019) and based on the JGI
}

\author{Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Otauri.JGI)
## list the contents that are loaded into memory
ls('package:TxDb.Otauri.JGI')
## show the db object that is loaded by calling it's name
TxDb.Otauri.JGI

}

\keyword{package}
\keyword{data}
