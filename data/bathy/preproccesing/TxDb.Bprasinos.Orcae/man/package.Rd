\name{TxDb.Bprasinos.Orcae}
\docType{package}

\alias{TxDb.Bprasinos.Orcae-package}
\alias{TxDb.Bprasinos.Orcae}
\alias{TxDb.Bprasinos.Orcae}


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
  This data package was made from resources at Orcae on
  2019-04-25 14:06:24 +0200 (Thu, 25 Apr 2019) and based on the Orcae
}

\author{Ana B. Romero-Losada, Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Bprasinos.Orcae)
## list the contents that are loaded into memory
ls('package:TxDb.Bprasinos.Orcae')
## show the db object that is loaded by calling it's name
TxDb.Bprasinos.Orcae

}

\keyword{package}
\keyword{data}
