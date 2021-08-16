\name{TxDb.Bprasinos.ORCAE}
\docType{package}

\alias{TxDb.Bprasinos.ORCAE-package}
\alias{TxDb.Bprasinos.ORCAE}
\alias{TxDb.Bprasinos.ORCAE}


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
  This data package was made from resources at ORCAE on
  2021-08-17 00:18:19 +0200 (Tue, 17 Aug 2021) and based on the ORCAE
}

\author{Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Bprasinos.ORCAE)
## list the contents that are loaded into memory
ls('package:TxDb.Bprasinos.ORCAE')
## show the db object that is loaded by calling it's name
TxDb.Bprasinos.ORCAE

}

\keyword{package}
\keyword{data}
