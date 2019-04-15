datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Vcarteri.eg <- function() showQCData("org.Vcarteri.eg", datacache)
org.Vcarteri.eg_dbconn <- function() dbconn(datacache)
org.Vcarteri.eg_dbfile <- function() dbfile(datacache)
org.Vcarteri.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Vcarteri.eg_dbInfo <- function() dbInfo(datacache)

org.Vcarteri.egORGANISM <- "Volvox carteri"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Vcarteri.eg.sqlite", package=pkgname, lib.loc=libname)
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
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Vcarteri.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Vcarteri.eg_dbconn())
}

