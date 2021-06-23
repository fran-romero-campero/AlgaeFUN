\name{TxDb.Ngaditana.JGI}
\docType{package}

\alias{TxDb.Ngaditana.JGI-package}
\alias{TxDb.Ngaditana.JGI}
\alias{TxDb.Ngaditana.JGI}


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
  2021-06-21 07:03:07 +0200 (Mon, 21 Jun 2021) and based on the JGI
}

\author{Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Ngaditana.JGI)
## list the contents that are loaded into memory
ls('package:TxDb.Ngaditana.JGI')
## show the db object that is loaded by calling it's name
TxDb.Ngaditana.JGI

}

\keyword{package}
\keyword{data}
