#!/usr/bin/env Rscript

###-----------------
### Timing
message("Processing started:", date(),"\n")
ptm<-proc.time()


###-----------------
### Arguments
###-----------------
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
                             defaults=list(
                                        blastvec = NULL,
                                        blastGeneA = NULL,
                                        blastGeneB = NULL,
                                        pafHost = NULL,
                                        minHost_mapQ = 10L,
                                        readLength = NULL,
                                        vectorSequence = NULL,
                                        minaln = 1L,
                                        MinDVDsides = 10e3L,
                                        outputdir = NULL,
                                        outFileName = NULL))

# Verbose:
message("Running AnnotateBACreads with the following arguments:\n")
message(paste(paste(names(args), 
                    ":", 
                    format(args), 
                    sep = " "), 
              collapse="\n"))

# Output file

if (is.null(args$outputdir)) {
  outdir <- dirname(args$blastvec)
} else {
  outdir <- args$outputdir
}

if (is.null(args$outFileName)) {
  outn <- paste0("ReadTypeTable_",basename(args$readLength), ".rds")
} else {
  outn <- basename(args$outFileName)
}

if (file.exists(file.path(outdir, outn))) {
    warning("Overwriting ", file.path(outdir, outn))
}

###-----------------
### Script
###-----------------
# Run AnnotateBACreads
anno <- NanoBAC::AnnotateBACreads(
                    blastvec = args$blastvec,
                    blastGeneA = args$blastGeneA,
                    blastGeneB = args$blastGeneB,
                    pafHost = args$pafHost,
                    minHost_mapQ = args$minHost_mapQ,
                    readLength = args$readLength,
                    vectorSequence = args$vectorSequence,
                    minaln = args$minaln,
                    MinDVDsides = args$MinDVDsides)

# Save the result to an rds file
saveRDS(anno, file.path(outdir, outn))


### Timing
###-----------------
message("Processing Ended:", date(),"\n")
message("Timing:\n")
proc.time()-ptm

# Done
