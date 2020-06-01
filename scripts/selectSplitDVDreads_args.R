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
                                        ReadClass = NULL,
                                        blastvec = NULL,
                                        FastaFile = NULL,
                                        MinDNAlength = 3e4L,
                                        WithGeneA = TRUE,
                                        WithGeneB = TRUE,
                                        outRDS = NULL,
                                        outReadNames = NULL,
                                        outFasta = NULL))

# Verbose:
message("Running selectSplitDVDreads with the following arguments:\n")
message(paste(paste(names(args), 
                    ":", 
                    format(args), 
                    sep = " "), 
              collapse="\n"))


# Arguments & output files
args$WithGeneA <- as.logical(args$WithGeneA)
args$WithGeneB <- as.logical(args$WithGeneB)

## ReadClass
if (is.null(args$ReadClass)) {
    stop("Please provide a ReadClass file")
}

if (is.null(args$blastvec)) {
    stop("Please provide a blastvec file")
}

if (is.null(args$FastaFile)) {
    stop("Please provide a FastaFile file")
}


if (is.null(args$outRDS)) {
  stop("Please provide a path to output the rds file (outRDS argument)")
} else {
  outRDS <- args$outRDS
}

if (file.exists(args$outRDS)) {
    warning("Overwriting ", args$outRDS)
}

if (is.null(args$outReadNames)) {
  stop("Please provide a path to output the list of DVD reads (outReadNames argument)")
} else {
  outReadNames <- args$outReadNames
}

if (file.exists(args$outReadNames)) {
    warning("Overwriting ", args$outReadNames)
}

if (is.null(args$outFasta)) {
  stop("Please provide a path to output the sequence of split DVD reads (outReadNames argument)")
} else {
  outFasta <- args$outFasta
}

if (file.exists(args$outFasta)) {
    warning("Overwriting ", args$outFasta)
}

## Silently load libraries
suppressMessages(require(Biostrings))
suppressMessages(require(GenomicRanges))
suppressMessages(require(NanoBAC))

###-----------------
### Script
###-----------------
# Select ans split DVD reads
selDVD <- NanoBAC::splitDVDreads(
                         ReadClass = args$ReadClass,
                         blastvec = args$blastvec,
                         FastaFile = args$FastaFile,
                         MinDNAlength = args$MinDNAlength,
                         WithGeneA = args$WithGeneA,
                         WithGeneB = args$WithGeneB
                         )

if (is.null(selDVD)) {
    
    ## When using a ligation protocol, it is frequent to have zero DVD reads
    ## Thus, we throw a warning but no error if there is no DVD reads
    warning("Zero DVD reads selected")
    ## Create empty files as result
    system(paste("touch", outRDS, outReadNames, outFasta, sep = " "))
    ## Quit R without error
    quit(save = "no", status = 0)

} else {
    
    ## save reads
    Biostrings::writeXStringSet(selDVD$ReadSequence, outFasta)
    ## save read names
    write.table(GenomeInfoDb::seqlevels(selDVD$ReadDefinition),
                outReadNames,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    ## save read definitions (as a GRanges)
    saveRDS(selDVD$ReadDefinition, outRDS)

    ## Print some stats about the reads found
    Nreads <- length(GenomeInfoDb::seqlevels(selDVD$ReadDefinition))
    NreadFragments <- length(selDVD$ReadDefinition)
    avgReadSize <- round(mean(GenomicRanges::width(selDVD$ReadDefinition), na.rm=TRUE), 1)
    sdReadSize <- round(sd(GenomicRanges::width(selDVD$ReadDefinition), na.rm=TRUE), 1)
    avgInsertSize <- round(mean(GenomicRanges::mcols(selDVD$ReadDefinition)$DNAlength, na.rm=TRUE), 1)
    sdInsertSize <- round(sd(GenomicRanges::mcols(selDVD$ReadDefinition)$DNAlength, na.rm=TRUE), 1)
    message("\nNumber of DVD reads selected: ", Nreads)
    message("Number of DV/VD fragments after splitting and size selection: ", NreadFragments) 
    message("Read Length (mean +/- sd): ", avgReadSize, " +/- ", sdReadSize)
    message("Insert Length (mean +/- sd): ", avgInsertSize, " +/- ", sdInsertSize)

}

### Timing
###-----------------
message("Processing Ended:", date(),"\n")
message("Timing:\n")
proc.time()-ptm

# Done
