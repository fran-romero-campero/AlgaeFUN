datacache <- new.env(hash=TRUE, parent=emptyenv())

org.MpusillaCCMP1545.eg <- function() showQCData("org.MpusillaCCMP1545.eg", datacache)
org.MpusillaCCMP1545.eg_dbconn <- function() dbconn(datacache)
org.MpusillaCCMP1545.eg_dbfile <- function() dbfile(datacache)
org.MpusillaCCMP1545.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.MpusillaCCMP1545.eg_dbInfo <- function() dbInfo(datacache)

org.MpusillaCCMP1545.egORGANISM <- "Micromonas pusillaCCMP1545"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.MpusillaCCMP1545.eg.sqlite", package=pkgname, lib.loc=libname)
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
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.MpusillaCCMP1545.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.MpusillaCCMP1545.eg_dbconn())
}

