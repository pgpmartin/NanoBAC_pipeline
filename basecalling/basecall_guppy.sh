#!/bin/bash

set -euo pipefail

##----------------
## Timing
##----------------
start=`date +%s`

##----------------
## Variables
##----------------

guppyContainer=${guppyContainer:-/N/project/NOR/SingIMG/guppy.sif}
#Container can be built with:
#  module load singularity/2.6.1
#  singularity build guppy.sif shub://pgpmartin/SingIMG:guppy_4.0.1
#inpath=${inpath}
outpath=${outpath:-"$inpath/guppy_basecalls"}
## if a folder $inpath/fast5s exist, guppy will do the basecalling in the fast5 files/folders that it contains
## otherwise guppy will do the basecalling on all the fast5 files located in $inpath/fast5_* folders

# Check that inpath is set and not empty
if [ -z "$inpath" ]; then
    printf "%s\n" "variable inpath is not set" >&2
    exit 1
fi

##----------------
#modules / softwares
##----------------
singversion=$(singularity --version 2>/dev/null)
errcode=$?

if [ "$errcode" -ne 0 ]; then
  module load singularity
  errLoadModule=$?
  if [ $errLoadModule -ne 0 ]; then
      printf "%s\n" "Could not load singularity module" >&2
      exit $errLoadModule
  fi
fi


##----------------
## Tests
##----------------

# Check that inpath exists
if [ ! -d "$inpath" ]; then
    printf "%s\n" "input directory $inpath does not exist" >&2
    exit 1
fi


# Check that guppy container exists
if [ ! -f "$guppyContainer" ]; then
    printf "%s\n" "$guppyContainer not found" >&2
    exit 1
fi

# Check that guppy_basecaller works an get its version
GuppyVersion=$(singularity exec --nv ${guppyContainer} guppy_basecaller --version)
errcode=$?
if [ "$errcode" -ne 0 ]; then
    printf "%s\n" "The guppy_basecaller function returned error $errcode" >&2
    exit $errcode
else
    GuppyVersion=$(echo $GuppyVersion | grep -o 'Version\ .*\,' | cut -c9- | sed 's/,$//')
fi


# If fatst5s folder is absent, check that fast5_* folders are present in inpath
if [ ! -d ${inpath}/fast5s ]; then
    numfast5=$(ls -1dA $inpath/fast5_*/ | wc -l)
    if [ "$numfast5" -eq 0 ]; then 
        printf "%s\n" "Input directory does not contain fast5_* folders" >&2
        exit 1
    fi
fi



# Check if outpath exists.  If it exists check that it is empty. If it doesn't exist, create it.
if [ -d "$outpath" ]; then
    if (ls -1qA "$outpath" 2>/dev/null | grep -q .); then
        printf "%s\n" "$outpath directory is not empty!" >&2
        exit 1
    fi
else
    mkdir -p ${outpath}
fi


##----------------
## Summary of variables
##----------------
printf "%s\n" "Using guppy version: ${GuppyVersion}"
printf "%s\n" "Input directory (containing fast5_* folders): ${inpath}"
printf "%s\n" "Output directory: ${outpath}"


##----------------
## Prepare fast5s folder if it does not exist
##----------------


if (ls -1qA ${inpath}/fast5s/* 2>/dev/null | grep -q .); then
    printf "%s\n%s\n" \
        "${inpath}/fast5s directory already exists and is not empty!" \
        "Using this directory for basecalling"
else
    mkdir -p "${inpath}/fast5s"
    for fn in `ls -d ${inpath}/fast5_*`; do
        ln -s $(realpath $fn) ${inpath}/fast5s/$(basename $fn)
    done
fi

##----------------
## Run guppy
##----------------

singularity exec \
  --nv \
  ${guppyContainer} \
  guppy_basecaller \
    -r \
    -c dna_r9.4.1_450bps_hac.cfg \
    -x 'cuda:0 cuda:1' \
    --num_callers 8 \
    --gpu_runners_per_device 4 \
    --chunks_per_runner 1600 \
    --input_path ${inpath}/fast5s \
    --save_path ${outpath}

##----------------
## Timing
##----------------
end=`date +%s`
printf "%s\n" "Duration: $((end-start)) seconds"
printf "%s\n" "Basecalling done!"
#Done
