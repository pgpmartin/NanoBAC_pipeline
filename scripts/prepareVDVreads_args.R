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
                                        inputVDVfasta = NULL,
                                        outputVDVfasta = NULL,
                                        outputVDVjunctions = NULL,
                                        FilterOutputfasta = TRUE,
                                        blastvec = NULL,
                                        RestrictionSite = "G^AATTC",
                                        VectorSequence = NULL,
                                        UnalignedVectorLength = 1000L,
                                        SideSeqSearch = 10L,
                                        Ncores = 1
                                        ))





# Verbose:
message("Running prepareVDVreads with the following arguments:\n")
message(paste(paste(names(args), 
                    ":", 
                    format(args), 
                    sep = " "), 
              collapse="\n"))

# Arguments
if (is.null(args$inputVDVfasta)) {
    stop("inputVDVfasta argument is NULL")
}

if (!file.exists(args$inputVDVfasta)) {
    stop("file not found: ", args$inputVDVfasta)
}

if (is.null(args$outputVDVfasta)) {
    stop("outputVDVfasta argument is NULL")
}

if (is.null(args$outputVDVjunctions)) {
    stop("outputVDVjunctions argument is NULL")
}

if (is.null(args$blastvec)) {
    stop("blastvec argument is NULL")
}

if (!file.exists(args$blastvec)) {
    stop("file not found: ", args$blastvec)
}

if (is.null(args$VectorSequence)) {
    stop("VectorSequence argument is NULL")
}

if (!file.exists(args$VectorSequence)) {
    stop("file not found: ", args$VectorSequence)
}

if (file.exists(args$outputVDVfasta)) {
    warning(args$outputVDVfasta, " already exists. Overwritting")
}

if (file.exists(args$outputVDVjunctions)) {
    warning(args$outputVDVjunctions, " already exists. Overwritting")
}


###-----------------
### Libraries
###-----------------
suppressMessages(require(Biostrings))
suppressMessages(require(NanoBAC))
suppressMessages(require(parallel))

###-----------------
### Script
###-----------------

## import VDV read sequences:
    readsDNA <- Biostrings::readDNAStringSet(args$inputVDVfasta)

    if (length(readsDNA) == 0) {
        stop("inputVDVfasta is empty")
    }

    if (length(unique(names(readsDNA))) != length(readsDNA)) {
        stop("VDV read names are not unique")
    }

## get read sizes
    ReadLength <- data.frame("ReadName" = names(readsDNA),
                             "ReadLength" = width(readsDNA))

## import and filter vector alignment (BLAST) results:
    ## import blast table
    ReadVecALN <- NanoBAC::readBlast(args$blastvec)
    ### Keep only alignments for VDV reads
    ReadVecALN <- ReadVecALN[ReadVecALN$SubjectACC %in% names(readsDNA),]
    ### Filter irrelevant alignments (those overlapping on >50% of their length)
    ReadVecALN <- ReadVecALN[NanoBAC::SelectSingularBlastALN(ReadVecALN, 
                                        rl = ReadLength, 
                                        threshold = 0.5), ]

## Import vector sequence
    VecSeq <- Biostrings::readDNAStringSet(args$VectorSequence)[[1]]
    VecSize <- length(VecSeq)
    
## Find the junctions for all reads
    VDVjunc <- parallel::mclapply(names(readsDNA),
                    NanoBAC::findVDVjunctions,
                    ReadDNA = readsDNA,
                    ReadVecAlign = ReadVecALN,
                    RestrictionSite = args$RestrictionSite,
                    VectorSequence = VecSeq,
                    UnalignedVectorLength = args$UnalignedVectorLength,
                    SideSeqSearch = args$SideSeqSearch,
                    replaceVectorSequence = TRUE,
                    mc.cores = args$Ncores)

# Format a table with info about the junctions
    juncTable <- data.frame(
                       ReadName = unlist(lapply(VDVjunc, 
                                                `[[`, "ReadName"), 
                                         use.names = FALSE),
                       Strand = unlist(lapply(VDVjunc, 
                                                `[[`, "Strand"), 
                                        use.names = FALSE),
                       InsertStart = unlist(lapply(VDVjunc, 
                                                `[[`, "InsertStart"), 
                                        use.names = FALSE),
                       InsertEnd = unlist(lapply(VDVjunc, 
                                                `[[`, "InsertEnd"), 
                                        use.names = FALSE)
                            )

# Save junction table
    saveRDS(juncTable, 
            args$outputVDVjunctions)

# Format a DNAStringSet with the "corrected" sequences
    newReads <- DNAStringSet(lapply(VDVjunc, `[[`, "correctedRead"))
    names(newReads) <- juncTable$ReadName

#Count the reads with/without junctions
    message("\nInitial number of reads: ", 
            nrow(juncTable))
    StartIsNotOK <- is.na(juncTable[,"InsertStart"]) # start of read
    message("Reads missing a junction at the beginning: ", 
            sum(StartIsNotOK))
    EndIsNotOK <- is.na(juncTable[,"InsertEnd"]) # end of read
    message("Reads missing a junction at the end: ", 
            sum(EndIsNotOK))
    isReadOK <- !StartIsNotOK & !EndIsNotOK #both ends OK
    message("Reads missing a junction on any side: ", 
            sum(StartIsNotOK | EndIsNotOK))
    message("Reads missing a junction on both sides: ", 
            sum(StartIsNotOK & EndIsNotOK))
    message("Reads with both junctions found: ", 
            sum(isReadOK))

# if required, filter the sequences where junctions were not found
    if (args$FilterOutputfasta) {
        message("Filtering reads. Number of VDV reads returned: ", 
                sum(isReadOK))
        newReads <- newReads[isReadOK]
    } else {
        message("No filtering. Number of VDV reads returned: ", 
                nrow(juncTable))
    }

    if (args$FilterOutputfasta) {
    ## Additional statistics on first/last insert bases when filtering is done
        message("\nBase frequency at the first insert base:")
        print(
        colSums(Biostrings::alphabetFrequency(
                                Biostrings::subseq(newReads, 
                                                    VecSize + 1, 
                                                    VecSize + 1), 
                                baseOnly=TRUE))
            )
        message("Base frequency at the last insert base:")
        print(
        colSums(Biostrings::alphabetFrequency(
                                Biostrings::subseq(newReads, 
                                    Biostrings::width(newReads)-VecSize, 
                                    Biostrings::width(newReads)-VecSize), 
                                baseOnly=TRUE))
            )
    }

# Save the corrected sequences
    Biostrings::writeXStringSet(newReads, 
                                args$outputVDVfasta, 
                                format = "fasta")


### Timing
###-----------------
message("\nProcessing Ended:", date())
message("Timing:")
proc.time()-ptm

# Done
