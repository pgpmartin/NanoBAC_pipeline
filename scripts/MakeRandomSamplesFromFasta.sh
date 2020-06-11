#!/bin/bash

# seqtk must be in your PATH to use this script
# In order to use the default values for directories, the fasta file is expected to be in ${baseDIR}/SelectedReads/VDVprepared
# Otherwise, provide at least baseDIR and possibly outDIR

#--------
# Variables
#--------

  ## Test if fastafile is unset
    if [ -z "$fastafile" ]
    then
      printf "%s\n" "The variable fastafile is not set" >&2
      exit 1
    fi

  ## Test if sampleName is unset
    if [ -z "$sampleName" ]
    then
      printf "%s\n" "The variable sampleName is not set" >&2
      exit 1
    fi

  ## Test if outDIR is unset
    if [ -z "$outDIR" ]
    then
      printf "%s\n" "The variable outDIR is not set" >&2
      exit 1
    fi

  ## Test if NumRndSample is unset
    if [ -z "$NumRndSample" ]
    then
      printf "%s\n" "The variable NumRndSample is not set" >&2
      exit 1
    fi

  ## Require at least 11 VDV reads to get a reliable BAC sequence
    ### Get the number of total sequences:
    TotNumSeq=$(grep "^>" ${fastafile} | wc -l)
    
    if (( TotNumSeq < 11 ))
    then
      printf "%s\n" "${sampleName} has only ${TotNumSeq} VDV sequences" >&2
      exit 1
    fi

  ## Number of sequences per random set (defaults to 11)
    NumSeqPerSample=${NumSeqPerSample:-11}

  ## Return the values of all variables
    mkdir -p ${outDIR}
    printf "%s\n" "${sampleName} has ${TotNumSeq} VDV sequences"
    printf "%s\n" "Sample Name: ${sampleName}"
    printf "%s\n" "Input fasta file: ${fastafile}"
    printf "%s\n" "Output directory: ${outDIR}"
    printf "%s\n" "Generating ${NumRndSample} random samples of size ${NumSeqPerSample}"


#--------
# Generate ${NumRndSample} random numbers
#--------

    count=1
    while [ "$count" -le $NumRndSample ]; do
        RNDnum[$count]=$RANDOM
        let "count += 1"
    done

#--------
# Create the random samples
#--------
printf "%s\n" "Generating random samples for ${sampleName}"
printf "\n%s\n" "Random seeds used:"

    for (( i=1; i<=$NumRndSample; i++ ))
    do
        seqtk sample \
            -s ${RNDnum[$i]} \
            ${fastafile} \
            ${NumSeqPerSample} > \
            ${outDIR}/${sampleName}_RS${i}.fa

        printf "%s\t%s\n" "Random Sample ${SampleIndex}" "${RNDnum[$i]}"

    done

printf "\n"

#Done

