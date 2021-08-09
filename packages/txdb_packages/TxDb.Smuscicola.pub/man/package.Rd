\name{TxDb.Smuscicola.pub}
\docType{package}

\alias{TxDb.Smuscicola.pub-package}
\alias{TxDb.Smuscicola.pub}
\alias{TxDb.Smuscicola.pub}


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
  2021-06-22 11:54:08 +0200 (Tue, 22 Jun 2021) and based on the pub
}

\author{Christina Arvanitidou, Francisco J. Romero-Campero <fran@us.es> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Smuscicola.pub)
## list the contents that are loaded into memory
ls('package:TxDb.Smuscicola.pub')
## show the db object that is loaded by calling it's name
TxDb.Smuscicola.pub

}

\keyword{package}
\keyword{data}
