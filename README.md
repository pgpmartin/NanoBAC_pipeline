# NanoBAC_pipeline

## Preparing to run the pipeline  

First create the NanoBAC conda environment
```bash
conda env create -f NanoBAC.yml
```

Activate the environment using

Then install R or make sure R (v>3.6) is in your PATH (e.g. `echo $PATH` and `which R`).    
In R, install the NanoBAC package using the following command:
```r
devtools::install_github("pgpmartin/NanoBAC")
```

Create a project folder, for example:
```bash
projPath="/DATA/BACproject"
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
Place the content of the `NanoBAC_pipeline/scripts/` folder in this folder and edit the `config.json` file to indicate the location of this script folder 


Create the subdirectories `data/raw/` and `log/`  
```
mkdir -p ${projPath}/data/raw
mkdir -p ${projPath}/log
```
In the `data/raw/` subdirectory, place the fastq file(s) with your Nanopore reads for your BAC(s) sequencing experiment (format must be: `sampleName.fastq` where `sampleName` is replaced by the actual name(s) of your sample(s))  

In the `projPath` directory, place the following files:  

  - `Snakefile`  
  - `config.json` (edit to your need)  
  - `cluster.json`  
  - `envs/` folder with all its content (`*.yaml` files)  
  
Your folder should look something like this:
```
BACproject
├── cluster.json
├── config.json
├── data
│   └── raw
│       ├── BAC01.fastq
│       ├── BAC02.fastq
│       └── BAC03.fastq
├── envs
│   ├── blast.yaml
│   ├── kalign2.yaml
│   ├── minimap2.yaml
│   ├── racon.yaml
│   ├── samtools.yaml
│   └── seqtk.yaml
├── log
├── myenvs
├── scripts
│   ├── AnnotateBACreads_args.R
│   ├── getNanoBACconsensus.R
│   ├── MakeRandomSamplesFromFasta.sh
│   ├── prepareVDVreads_args.R
│   ├── selectlongVDreads_args.R
│   ├── selectSplitDVDreads_args.R
│   └── selectVDVreads_args.R
└── Snakefile
```

Make sure the scripts are executable with `chmod u+x ${projPath}/scripts/*` (same for the SnakeFile)


## Running the pipeline  

In all cases activate your conda environment and make sure R is in your path with the NanoBAC package installed:
```bash
source activate NanoBAC
which R #check that R is installed (and have the NanoBAC R package installed too)
```


Run locally (or interactively on a cluster node) using 4 threads:
```bash
snakemake \
  -j 4 \
  --use-conda \
  --conda-prefix ${projPath}/myenvs
```


Run on a cluster using **PBS/Torque** scheduler:
```bash
snakemake \
  -j 100 \
  --cluster-config cluster.json \
  --use-conda \
  --conda-prefix ${projPath}/myenvs/myenvs \
  --cluster "qsub -l nodes=1:ppn={threads},walltime={cluster.time},mem={cluster.mem},vmem={cluster.mem} -o {cluster.logfolder} -e {cluster.logfolder}"
```


Run on a cluster using **SLURM** scheduler:
```bash
snakemake \
  -j 100 \
  --cluster-config cluster.json \
  --use-conda \
  --conda-prefix ${projPath}/myenvs \
  --cluster "sbatch -o {cluster.logfolder} -e {cluster.logfolder} --mem-per-cpu={cluster.mem} --time={cluster.time} --ntasks=1 --cpus-per-task={threads}"
```

If running Snakemake in an interactive session and submitting batch job, it may be necessary to add the option `--latency-wait 120 all` as explained [here](https://hpc.nih.gov/apps/snakemake.html).
