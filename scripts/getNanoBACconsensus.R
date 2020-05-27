#!/usr/bin/env Rscript

###-----------------
### Timing
###-----------------
message("Processing started:", date(),"\n")
ptm<-proc.time()

###-----------------
### Libraries
###-----------------
suppressMessages(require(NanoBAC))

###-----------------
### Arguments
###-----------------
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
                             defaults=list(
                                           msfFile=NULL,
                                           outputfile=NULL,
                                           consName = NULL
                                        ))

if (is.null(args$outputfile)) {
    stop("outputfile is NULL")
}

if (is.null(args$msfFile)) {
    stop("msfFile is NULL")
}

if (!file.exists(args$msfFile)) {
    stop(args$msfFile, " file not found")
}

if (file.exists(args$outputfile)) {
    warning(args$outputfile, " file exists. Overwriting")
}

if (is.null(args$consName)) {
    consName <- "NanoBAC"
} else {
    consName <- args$consName
}

cons <- NanoBAC::consensusFromMSF(args$msfFile,
                                  threshold = 0.5,
                                  removegaps = TRUE,
                                  consname = consName)

Biostrings::writeXStringSet(cons, 
                            args$outputfile)

###-----------------
### Timing
###-----------------
message("Processing Ended:", date(),"\n")
message("Timing:\n")
proc.time()-ptm

#Done
