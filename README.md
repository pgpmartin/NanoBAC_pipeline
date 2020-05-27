# NanoBAC_pipeline

## Preparing to run the pipeline  

First create the NanoBAC pipeline
```bash
conda env create -f NanoBAC.yml
```

Activate the environment using

Then install R or make sure R (v>3.6) is in your PATH (e.g. `echo $PATH` and `which R`).    
In R, install the NanoBAC package using the following command:
```r
devtools::install_github("pascmart/NanoBAC")
```

Create a project folder, for example:
```bash
projPath="/N/slate/pascmart/PROJECTS/TestBACseq"
mkdir -p ${projPath}
cd ${projPath}
```

Create a folder in which the conda environments will be built:  
```bash
mkdir -p ${projPath}/myenvs
```

Create a folder in which the scripts will be placed:
```bash
mkdir -p ${projPath}/scripts
```
Place the content of the `scripts/` folder in this folder and edit the config.json file to indicate the location of this script folder 


Create a folder for your first BAC with the subdirectories `raw/` and `log`  
```
mkdir -p ${projPath}/myFirstBAC/raw
mkdir -p ${projPath}/myFirstBAC/log
```
In the `raw/` subdirectory, place the fastq file with your Nanopore reads for your BAC sequencing experiment (format must be: `sampleName.fastq` where `sampleName` is replaced by the actual name of your sample)  


## Running the pipeline  

For PBS/Torque:
```bash
snakemake \
  -j 100 \
  --cluster-config cluster.json \
  --use-conda \
  --conda-prefix /N/slate/pascmart/PROJECTS/TestBACseq/myenvs \
  --cluster "qsub -l nodes=1:ppn={threads},walltime={cluster.time},mem={cluster.mem},vmem={cluster.mem} -o {cluster.logfolder} -e {cluster.logfolder}"
```

For SLURM:
```bash
snakemake \
  -j 100 \
  --cluster-config cluster.json \
  --use-conda \
  --conda-prefix /N/slate/pascmart/PROJECTS/TestBACseq/myenvs \
  --cluster "sbatch -o {cluster.logfolder} -e {cluster.logfolder} --mem-per-cpu={cluster.mem} --time={cluster.time} --ntasks=1 --cpus-per-task={threads}"
```
