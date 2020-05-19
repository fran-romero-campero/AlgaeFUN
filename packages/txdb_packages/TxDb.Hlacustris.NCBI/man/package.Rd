\name{TxDb.Hlacustris.NCBI}
\docType{package}

\alias{TxDb.Hlacustris.NCBI-package}
\alias{TxDb.Hlacustris.NCBI}
\alias{TxDb.Hlacustris.NCBI}


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
  This data package was made from resources at NCBI on
  2020-05-19 13:17:35 +0200 (Tue, 19 May 2020) and based on the NCBI
}

\author{Ana B. Romero-Losada, Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Hlacustris.NCBI)
## list the contents that are loaded into memory
ls('package:TxDb.Hlacustris.NCBI')
## show the db object that is loaded by calling it's name
TxDb.Hlacustris.NCBI

}

\keyword{package}
\keyword{data}
