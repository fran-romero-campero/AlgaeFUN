\name{TxDb.Dsalina.Phytozome}
\docType{package}

\alias{TxDb.Dsalina.Phytozome-package}
\alias{TxDb.Dsalina.Phytozome}
\alias{TxDb.Dsalina.Phytozome}


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
  2020-05-27 16:44:32 +0200 (Wed, 27 May 2020) and based on the Phytozome
}

\author{Ana B. Romero-Losada, Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Dsalina.Phytozome)
## list the contents that are loaded into memory
ls('package:TxDb.Dsalina.Phytozome')
## show the db object that is loaded by calling it's name
TxDb.Dsalina.Phytozome

}

\keyword{package}
\keyword{data}
