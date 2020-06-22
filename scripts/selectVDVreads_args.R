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
                                        SizeTolerance = 0.05,
                                        WithGeneA = NULL,
                                        WithGeneB = NULL,
                                        ignoredReads = NULL,
                                        MaxClusters = 10L,
                                        makePlot = TRUE,
                                        plotVar = "InsertLength",
                                        plotRes = 150L,
                                        plotHeight = 4,
                                        plotWidth = 4,
                                        plotFormat = c("png", "pdf"),
                                        outPlot = NULL,
                                        outRDS = NULL,
                                        outVDVnames = NULL))

# parse the ignreReads argument if not null
if (toupper(args$ignoredReads) == "NULL") {
    args$ignoredReads = NULL
}
if (!is.null(args$ignoredReads) && toupper(args$ignoredReads) == "NA") {
    args$ignoredReads = NULL
}
if (!is.null(args$ignoredReads)) {
    args$ignoredReads <- gsub(" ", "", strsplit(args$ignoredReads, ",")[[1]])
    message("Number of ignored reads: ", length(args$ignoredReads))
}

# Verbose:
message("Running selectVDVreads with the following arguments:\n")
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


## Plot file(s)
args$plotVar <- match.arg(args$plotVar,
                          choices = c("ReadLength", "InsertLength"),
                          several.ok = FALSE)

args$plotFormat <- match.arg(args$plotFormat,
                        choices = c("png", "pdf"),
                        several.ok = TRUE)


## Output objects
if (args$makePlot) {

    if (is.null(args$outPlot)) {
        stop("Please provide a path to output the plot (outPlot argument)")
    } else {
        outPlotBase <- gsub("\\.(jpg|jpeg|png|pdf|svg)$", "", 
                            args$outPlot)
    }
}

if (is.null(args$outRDS)) {
  stop("Please provide a path to output the rds file (outRDS argument)")
} else {
  outRDS <- args$outRDS
}

if (file.exists(args$outRDS)) {
    warning("Overwriting ", args$outRDS)
}

if (is.null(args$outVDVnames)) {
  stop("Please provide a path to output the list of VDV reads (outVDVnames argument)")
} else {
  outVDVnames <- args$outVDVnames
}

if (file.exists(args$outVDVnames)) {
    warning("Overwriting ", args$outVDVnames)
}


###-----------------
### Script
###-----------------
# Run selectVDVreads

pdf(NULL)
selvdv <- NanoBAC::selectVDVreads(
                         ReadClass = args$ReadClass,
                         SizeTolerance = args$SizeTolerance,
                         WithGeneA = args$WithGeneA,
                         WithGeneB = args$WithGeneB,
                         ignoredReads = args$ignoredReads,
                         MaxClusters = args$MaxClusters,
                         makePlot = args$makePlot,
                         plotVar = args$plotVar
                         )
dev.off()

# save plot
if (args$makePlot) {

    require(ggplot2)
    
    if ("png" %in% args$plotFormat) {
        png(paste0(outPlotBase, ".png"),
            res = args$plotRes,
            width = args$plotWidth * args$plotRes,
            height = args$plotHeight * args$plotRes)
        print(selvdv$VDVlengthPlot)
        dev.off()
    }

    if ("pdf" %in% args$plotFormat) {
        pdf(paste0(outPlotBase, ".pdf"),
            width = args$plotWidth,
            height = args$plotHeight,
            useDingbats = FALSE)
        print(selvdv$VDVlengthPlot)
        dev.off()
    }
}

# save RDS file
saveRDS(selvdv, outRDS)

# Save the list of VDV reads
write.table(as.character(selvdv$VDVreads$ReadName[selvdv$VDVreads$Selected]),
            outVDVnames,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# Print the estimated BAC size (insert length)
message("\nThe estimated length of the DNA insert (in bp) is: ", 
        selvdv$InsertSizeEstimate$BACsize, "\n")

### Timing
###-----------------
message("Processing Ended:", date(),"\n")
message("Timing:\n")
proc.time()-ptm

# Done
