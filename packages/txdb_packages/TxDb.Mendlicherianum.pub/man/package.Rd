\name{TxDb.Mendlicherianum.pub}
\docType{package}

\alias{TxDb.Mendlicherianum.pub-package}
\alias{TxDb.Mendlicherianum.pub}
\alias{TxDb.Mendlicherianum.pub}


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
  This data package was made from resources at pub on
  2021-06-08 08:01:55 +0200 (Tue, 08 Jun 2021) and based on the pub
}

\author{Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Mendlicherianum.pub)
## list the contents that are loaded into memory
ls('package:TxDb.Mendlicherianum.pub')
## show the db object that is loaded by calling it's name
TxDb.Mendlicherianum.pub

}

\keyword{package}
\keyword{data}
