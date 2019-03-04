datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Dsalina.eg <- function() showQCData("org.Dsalina.eg", datacache)
org.Dsalina.eg_dbconn <- function() dbconn(datacache)
org.Dsalina.eg_dbfile <- function() dbfile(datacache)
org.Dsalina.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Dsalina.eg_dbInfo <- function() dbInfo(datacache)

org.Dsalina.egORGANISM <- "Dunaliella salina"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Dsalina.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Dsalina.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Dsalina.eg_dbconn())
}

