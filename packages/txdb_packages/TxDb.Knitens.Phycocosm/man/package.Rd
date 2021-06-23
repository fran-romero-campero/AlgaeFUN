\name{TxDb.Knitens.Phycocosm}
\docType{package}

\alias{TxDb.Knitens.Phycocosm-package}
\alias{TxDb.Knitens.Phycocosm}
\alias{TxDb.Knitens.Phycocosm}


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
  This data package was made from resources at Phycocosm on
  2021-06-23 16:00:12 +0200 (Wed, 23 Jun 2021) and based on the Phycocosm
}

\author{Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Knitens.Phycocosm)
## list the contents that are loaded into memory
ls('package:TxDb.Knitens.Phycocosm')
## show the db object that is loaded by calling it's name
TxDb.Knitens.Phycocosm

}

\keyword{package}
\keyword{data}
