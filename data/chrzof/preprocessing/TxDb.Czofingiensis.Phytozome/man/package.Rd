\name{TxDb.Czofingiensis.Phytozome}
\docType{package}

\alias{TxDb.Czofingiensis.Phytozome-package}
\alias{TxDb.Czofingiensis.Phytozome}
\alias{TxDb.Czofingiensis.Phytozome}


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
  2019-04-23 11:52:34 +0200 (Tue, 23 Apr 2019) and based on the Phytozome
}

\author{Pedro de los Reyes, Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Czofingiensis.Phytozome)
## list the contents that are loaded into memory
ls('package:TxDb.Czofingiensis.Phytozome')
## show the db object that is loaded by calling it's name
TxDb.Czofingiensis.Phytozome

}

\keyword{package}
\keyword{data}
