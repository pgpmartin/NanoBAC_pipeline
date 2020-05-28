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
                                        minLength = 3e4L,
                                        WithGeneA = TRUE,
                                        WithGeneB = TRUE,
                                        isHostAlign = FALSE,
                                        outRDS = NULL,
                                        outReadNames = NULL))

# Verbose:
message("Running selectlongVDreads with the following arguments:\n")
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

if (is.null(args$outRDS)) {
  stop("Please provide a path to output the rds file (outRDS argument)")
} else {
  outRDS <- args$outRDS
}

if (file.exists(args$outRDS)) {
    warning("Overwriting ", args$outRDS)
}

if (is.null(args$outReadNames)) {
  stop("Please provide a path to output the list of VDV reads (outReadNames argument)")
} else {
  outReadNames <- args$outReadNames
}

if (file.exists(args$outReadNames)) {
    warning("Overwriting ", args$outReadNames)
}


###-----------------
### Script
###-----------------
# Select long VD reads
selLongVD <- NanoBAC::FilterBACreads(
                         ReadClass = args$ReadClass,
                         readtype = "VD",
                         MinReadLength = args$minLength,
                         alnGeneA = args$WithGeneA,
                         alnGeneB = args$WithGeneB,
                         isHostAlign = args$isHostAlign
                         )

# Throw an error if 0 reads are selected
if (is.null(selLongVD)) {
    stop("Zero long VD reads selected")
} else {
    
    # save RDS file
    saveRDS(selLongVD, outRDS)

    # Save the list of long VD reads
    write.table(as.character(selLongVD$ReadName),
                outReadNames,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)

    # Print some stats about the reads found
    Nreads <- nrow(selLongVD)
    avgReadSize <- mean(selLongVD$ReadLength, na.rm=TRUE)
    sdReadSize <- sd(selLongVD$ReadLength, na.rm=TRUE)
    avgInsertSize <- mean(selLongVD$LongestDNA, na.rm=TRUE)
    sdInsertSize <- sd(selLongVD$LongestDNA, na.rm=TRUE)
    message("\nNumber of long VD reads selected: ", Nreads)
    message("Read Length (mean +/- sd): ", avgReadSize, " +/- ", sdReadSize)
    message("Insert Length (mean +/- sd): ", avgInsertSize, " +/- ", sdInsertSize)
}

### Timing
###-----------------
message("Processing Ended:", date(),"\n")
message("Timing:\n")
proc.time()-ptm

# Done
