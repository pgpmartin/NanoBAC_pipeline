configfile:
    "config.json"

shell.executable("/bin/bash")

localrules: all, Rename_Reads, RenamedReads_fq2fa, ReadLengthTable, AlignVector, AlignGeneA,
            AlignGeneB, Annotate_Reads, Select_VDVreads, selRead_Get_fastq, selRead_fastq2fasta, Prepare_VDV_reads,
            DefineRandomSampleParam, RandomSampling_VDVreads, consensus_VDVrandomSets, consensus_fromVDVcons, Select_longVDreads, SelectSplit_DVDreads,
            merge_VD_DV_reads, Polishing_step01, Polishing_step02

# Get sample names
SAMPLENAME, = glob_wildcards("data/raw/{sample}.fastq")

#modules
import os

# create log directory if it does not exist
LOGDIR = "log"
if not(os.path.exists(LOGDIR)):
    os.mkdir(LOGDIR)

# get the full path of the script directory
SCRIPTDIR = os.path.abspath(config['scriptDIR'])

# function to count the prepared VDV reads from a fasta file and ouput 3 values:
# 1st value is the number of VDV reads selected
# 2nd value is the number of random sets
# 3rd value is the number of sequences per random sets
def VDVcount(fafile):
    nvdv = len([1 for line in open(fafile) if line.startswith(">")])
#    print("Number of VDV reads: " + str(nvdv))
    if nvdv <= 10:
        print("There are " + nvdv + " VDV reads (<11 so stopping here)")
        exit(1)
    elif nvdv < 20:
        res = (nvdv, 7, 6)
    elif nvdv < 50:
        res = (nvdv, 11, 11)
    elif nvdv < 100:
        res = (nvdv, 21, 11)
    else:
        res = (nvdv, 21, 21)

    return(res)



rule all:
    input:
        expand("tables/ReadNames/{sample}_ReadNameTable.tsv", sample=SAMPLENAME),
        expand("data/renamed/{sample}.fq.gz", sample=SAMPLENAME),
#        expand("data/renamed/{sample}.fa", sample=SAMPLENAME),
#        expand("tables/ReadLength/{sample}_ReadLength.tsv", sample=SAMPLENAME),
#        expand("tables/ReadClass/{sample}_ReadClass.rds", sample=SAMPLENAME),
#        expand("align/Blast/results/{sample}_vectorblast.res", sample=SAMPLENAME),
#        expand("align/Blast/results/{sample}_geneAblast.res", sample=SAMPLENAME),
#        expand("align/Blast/results/{sample}_geneBblast.res", sample=SAMPLENAME),
#        expand("align/minimap2/host/{sample}.paf", sample=SAMPLENAME),
        expand("align/minimap2/host/{sample}_alignment.stats", sample=SAMPLENAME),
#        expand("SelectedReads/VDV/{sample}_Sel_VDV.rds", sample=SAMPLENAME),
#        expand("SelectedReads/VDV/{sample}_Sel_VDV.fa", sample=SAMPLENAME),
        expand("SelectedReads/VDV/{sample}_Sel_VDV.fq.gz", sample=SAMPLENAME),
        expand("SelectedReads/VDVprepared/{sample}_SelectedPreparedVDVreads.fa", sample=SAMPLENAME),
#        expand("SelectedReads/VDVprepared/RandomSamples/{sample}_RS{NUM}.fa", sample=SAMPLENAME, NUM = RNDSETS),
#        expand("align/kalign/VDVrndSets/{sample}_RS{NUM}.msf", sample=SAMPLENAME, NUM=RNDSETS),
#        expand("align/kalign/VDVrndSets/{sample}_RS{NUM}_cons.fa", sample=SAMPLENAME, NUM=RNDSETS),
#        expand("align/kalign/ConsFromVDV/{sample}_consAlign.msf", sample=SAMPLENAME),
#        expand("assembly/{sample}_InitialConsensus.fa", sample=SAMPLENAME),
#        expand("SelectedReads/longVD/{sample}_Sel_longVD.rds", sample=SAMPLENAME),
#        expand("SelectedReads/longVD/{sample}_Sel_longVD.fa", sample=SAMPLENAME),
#        expand("SelectedReads/longVD/{sample}_Sel_longVD.fq.gz", sample=SAMPLENAME),
#        expand("SelectedReads/DVD/{sample}_SelSplit_DVD.fa", sample=SAMPLENAME),
#        expand("SelectedReads/{sample}_VDreadsForPolishing.fa", sample=SAMPLENAME),
        expand("assembly/final_assembly_{sample}.fa", sample=SAMPLENAME)


# Create a table with new names and old names
rule Read_NameMapping:
    message: "Storing old and new names for {input} in {output}."
    input:
        "data/raw/{sample}.fastq"
    output:
        "tables/ReadNames/{sample}_ReadNameTable.tsv"
    params:
        sampleName = "{sample}"
    threads: 1
    shell:
        """
        cat {input} | \
          sed -n '1~4p' | \
          awk 'BEGIN {{ cntr = 0 }} {{ cntr++ ; print \"{params.sampleName}R\"cntr, \"\\t\", $0 }}'  > \
          {output}
        """

# Create a fastq file with the new (simple) names
rule Rename_Reads:
    message: "Renaming reads for {wildcards.sample}"
    input:
        "data/raw/{sample}.fastq"
    params:
        sampleName = "{sample}"
    output:
        temp("data/renamed/{sample}.fq")
    shell:
        """
        cat {input} | \
          sed -n '1~4s/^@.*$/@zyxabcdefzyx/p;2~4p;3~4p;0~4p' | \
          awk 'BEGIN {{ cntr = 0 }} /^@zyxabcdefzyx$/ {{ cntr++ ; print \"@{params.sampleName}R\"cntr }} !/^@zyxabcdefzyx$/ {{ print $0 }}' > \
          {output}
        """


# Create a fasta file with the new (simple) names
rule RenamedReads_fq2fa:
    message: "Converting {input} to fasta format"
    input:
        "data/renamed/{sample}.fq"
    output:
        "data/renamed/{sample}.fa"
    shell:
        """
        sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
        """

# Compress fastq file
rule Compress_Renamed_FastQ:
    message: "Compressing {input}"
    input: "data/renamed/{sample}.fq"
    output: "data/renamed/{sample}.fq.gz"
    threads: 1
    shell: " gzip -c {input} > {output} "


# Get a table with read lengths
rule ReadLengthTable:
    message: "Extracting Read Length for {input}"
    input:
        "data/renamed/{sample}.fa"
    output:
        "tables/ReadLength/{sample}_ReadLength.tsv"
    shell:
        """
        cat {input} | \
          awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) \"\\t\"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' - > \
          {output}
        """


# Create Read database for Blast
#makeblastdb creates variable file names depending on the size of the database
#Using the solution found here: https://stackoverflow.com/questions/59288002/missingoutputexception-and-latency-wait-error-with-snakemake
rule CreateBlastDatabase:
    message: "Creating BLAST database for {input}"
    input:
        "data/renamed/{sample}.fa"
    output:
        done = touch("align/Blast/db/{sample}.makeblastdb.done")
    log:
        "log/{sample}_makeblastdb.log"
    threads: 1
    conda:
        "envs/blast.yaml"
#    envmodule:
#        "blast/2.9.0"
    shell:
        """
        makeblastdb \
          -in {input} \
          -out align/Blast/db/{wildcards.sample}_blastdb \
          -parse_seqids \
          -dbtype nucl \
          -logfile {log}
        """

# Align the vector
rule AlignVector:
    message: "Aligning vector sequence on the reads from {wildcards.sample}"
    input:
        query = config['VectorFasta'],
        db_done = "align/Blast/db/{sample}.makeblastdb.done"
    output:
        "align/Blast/results/{sample}_vectorblast.res"
    log:
        "log/{sample}_blastVector.log"
    conda: "envs/blast.yaml"
#    envmodules: "blast/2.9.0"
    threads: 1
    shell:
        """
        blastn \
          -task megablast \
          -query {input.query} \
          -db align/Blast/db/{wildcards.sample}_blastdb \
          -out {output} \
          -outfmt 6 \
          -num_alignments 100000000 \
        >> {log} 2>&1
        """

# Align to gene A
rule AlignGeneA:
    message: "Aligning GeneA on the reads from {wildcards.sample}"
    input:
        query = config['GeneAFasta'],
        db_done = "align/Blast/db/{sample}.makeblastdb.done"
    output:
        "align/Blast/results/{sample}_geneAblast.res"
    log:
        "log/{sample}_blastGeneA.log"
    conda: "envs/blast.yaml"
#    envmodules: "blast/2.9.0"
    threads: 1
    shell:
        """
        blastn \
          -task megablast \
          -query {input.query} \
          -db align/Blast/db/{wildcards.sample}_blastdb \
          -out {output} \
          -outfmt 6 \
          -num_alignments 100000000 \
        >> {log} 2>&1
        """

# Align to gene B
rule AlignGeneB:
    message: "Aligning GeneB on the reads from {wildcards.sample}"
    input:
        query = config['GeneBFasta'],
        db_done = "align/Blast/db/{sample}.makeblastdb.done"
    output:
        "align/Blast/results/{sample}_geneBblast.res"
    log:
        "log/{sample}_blastGeneB.log"
    conda: "envs/blast.yaml"
#    envmodules: "blast/2.9.0"
    threads: 1
    shell:
        """
        blastn \
          -task megablast \
          -query {input.query} \
          -db align/Blast/db/{wildcards.sample}_blastdb \
          -out {output} \
          -outfmt 6 \
          -num_alignments 100000000 \
        >> {log} 2>&1
        """

# Align to host genome
rule Align2Host:
    message: "Aligning reads from {wildcards.sample} to host genome"
    input:
        reads = "data/renamed/{sample}.fq",
        hostGenome = config['hostGenomeFasta']
    output:
        temp("align/minimap2/host/{sample}.sam")
    log:
        "log/{sample}_hostalign.log"
    conda: "envs/minimap2.yaml"
#    envmodules: "minimap2/2.17"
    threads: 4
    shell:
        """
        minimap2 \
          -ax map-ont \
          -t {threads} \
          -L \
          {input.hostGenome} \
          {input.reads} > \
          {output}
        """

# Convert to paf (only convert the primary and supplementary alignments with the -p argument )
rule Hostsam2paf:
    message: "Convert {input} to PAF"
    input:
        "align/minimap2/host/{sample}.sam"
    output:
        "align/minimap2/host/{sample}.paf"
    threads: 1
    conda: "envs/minimap2.yaml"
#    envmodules: "minimap2/2.17"
    shell:
        """
        paftools.js sam2paf -p {input} > {output}
        """

# Gather some stats
rule HostAlignStats:
    message: "Host alignment stats for {input}"
    input:
        "align/minimap2/host/{sample}.sam"
    params:
        sampleName = "{sample}"
    output:
        "align/minimap2/host/{sample}_alignment.stats"
    threads: 1
    conda: "envs/samtools.yaml"
#    envmodules: "samtools/1.9"
    shell:
        """
        TotalAlignments=$(samtools view -c {input})
        UniqueReads=$(grep -v "^@(HD\|SQ\|RG\|PG\|CO)" {input} | awk '{{ print $1 }}' | sort | uniq | wc -l)
        AlignedReads=$(samtools view -F 260 {input} | grep -v "^@(HD\|SQ\|RG\|PG\|CO)" | awk '{{ print $1 }}' | sort | uniq | wc -l)
        AlignedReadsHighQ=$(samtools view -F 260 -q 60 {input} | grep -v "^@(HD\|SQ\|RG\|PG\|CO)" | awk '{{ print $1 }}' | sort | uniq | wc -l)
        UnmappedReads=$(samtools view -f 4 {input} | grep -v "^@(HD\|SQ\|RG\|PG\|CO)" | awk '{{ print $1 }}' | sort | uniq | wc -l)
        SecondaryAlignments=$(samtools view -f 256 {input} | grep -v "^@(HD\|SQ\|RG\|PG\|CO)" | awk '{{ print $1 }}' | sort | uniq | wc -l)

        printf "%s\\t" "Sample" > {output}
        printf "%s\\t" "Total_Alignments" >> {output}
        printf "%s\\t" "Total_Reads" >> {output}
        printf "%s\\t" "Aligned_Reads" >> {output}
        printf "%s\\t" "Aligned_Reads_Q60" >> {output}
        printf "%s\\t" "Unmapped_Reads" >> {output}
        printf "%s\\n" "ReadsWithSecondaryAlignment" >> {output}

        printf "%s\\t" "{params.sampleName}" >> {output}
        printf "%s\\t" "${{TotalAlignments}}" >> {output}
        printf "%s\\t" "${{UniqueReads}}" >> {output}
        printf "%s\\t" "${{AlignedReads}}" >> {output}
        printf "%s\\t" "${{AlignedReadsHighQ}}" >> {output}
        printf "%s\\t" "${{UnmappedReads}}" >> {output}
        printf "%s\\n" "${{SecondaryAlignments}}" >> {output}
        """

# Filter SAM and convert to BAM
rule Filter_HostAlignment_SAM2BAM:
    message: "Remove unmapped reads and secondary alignment from {input}, sort and convert to BAM"
    input:
        "align/minimap2/host/{sample}.sam"
    output:
        "align/minimap2/host/{sample}.bam"
    threads: 4
    conda: "envs/samtools.yaml"
#    envmodules: "samtools/1.9"
    shell:
        """
        if [ "{threads}" -ge 2 ];
        then
          sortthreads="-@ $(( {threads} - 1 ))"
        else
          sortthreads=""
        fi

        samtools view \
          -bh \
          -F 260 \
          {input} | \
        samtools sort \
          ${sortthreads} \
          -o {output}
        """

# Annotate the reads
rule Annotate_Reads:
    message: "Annotating reads for {wildcards.sample}"
    input:
        blastvec = "align/Blast/results/{sample}_vectorblast.res",
        blastGeneA = "align/Blast/results/{sample}_geneAblast.res",
        blastGeneB = "align/Blast/results/{sample}_geneBblast.res",
        pafHost = "align/minimap2/host/{sample}.paf",
        readLength = "tables/ReadLength/{sample}_ReadLength.tsv"
    params:
        vectorSequence = config['VectorFasta'],
        minHost_mapQ = 10,
        minaln = 1,
        MinDVDsides = 10000
    output:
        "tables/ReadClass/{sample}_ReadClass.rds"
    log:
        "log/{sample}_AnnotateReads.log"
#    envmodules: "r/3.6.0"
    singularity:
        "shub://pgpmartin/SingIMG:r36_nanobac"
    threads: 1
    shell:
        """
        Rscript {SCRIPTDIR}/AnnotateBACreads_args.R \
          -blastvec={input.blastvec} \
          -blastGeneA={input.blastGeneA} \
          -blastGeneB={input.blastGeneB} \
          -pafHost={input.pafHost} \
          -minHost_mapQ={params.minHost_mapQ} \
          -readLength={input.readLength} \
          -vectorSequence={params.vectorSequence} \
          -minaln={params.minaln} \
          -MinDVDsides={params.MinDVDsides} \
          -outputdir=$(dirname {output}) \
          -outFileName=$(basename {output}) \
        >> {log} 2>&1
        """

# Select VDV Reads
rule Select_VDVreads:
    message: "Selecting VDV reads for {wildcards.sample}"
    input:
        "tables/ReadClass/{sample}_ReadClass.rds"
    params:
        SizeTolerance = 0.05,
        WithGeneA = True,
        WithGeneB = True,
        ignoredReads = "NULL",
        MaxClusters = 10,
        plotVar = "InsertLength",
        plotRes = 150,
        plotHeight = 6,
        plotWidth = 6,
        plotFormat = "png"
    output:
        outPlot = "Plots/{sample}_VDVreadSelection.png",
        outRDS = "SelectedReads/VDV/{sample}_ALL_VDV.rds",
        outVDVnames = "SelectedReads/VDV/{sample}_Sel_VDV_ReadNames.tsv"
    log:
        "log/{sample}_Sel_VDV_Selection.log"
#    envmodules: "r/3.6.0"
    singularity:
        "shub://pgpmartin/SingIMG:r36_nanobac"
    threads: 1
    shell:
        """
        Rscript {SCRIPTDIR}/selectVDVreads_args.R \
          -ReadClass={input} \
          -SizeTolerance={params.SizeTolerance} \
          -WithGeneA={params.WithGeneA} \
          -WithGeneB={params.WithGeneB} \
          -ignoredReads={params.ignoredReads} \
          -MaxClusters={params.MaxClusters} \
          -plotVar={params.plotVar} \
          -plotRes={params.plotRes} \
          -plotHeight={params.plotHeight} \
          -plotWidth={params.plotWidth} \
          -plotFormat={params.plotFormat} \
          -outPlot={output.outPlot} \
          -outRDS={output.outRDS} \
          -outVDVnames={output.outVDVnames} \
        >> {log} 2>&1
        """


# Obtain fastq of selected reads
rule selRead_Get_fastq:
    message: "Obtain fastq file of selected reads for {wildcards.readSelection}/{wildcards.sample}"
    input:
        rawfastq = "data/renamed/{sample}.fq",
        ReadSelection = "SelectedReads/{readSelection}/{sample}_Sel_{readSelection}_ReadNames.tsv"
    output:
        temp("SelectedReads/{readSelection}/{sample}_Sel_{readSelection}.fq")
    conda: "envs/seqtk.yaml"
#    envmodules: "seqtk/1.3"
    threads: 1
    shell:
        """
        seqtk subseq \
            {input.rawfastq} \
            {input.ReadSelection} > \
            {output}
        """


# Convert fastq to fasta
rule selRead_fastq2fasta:
    message: "Convert {input} to fasta"
    input:
        "SelectedReads/{readSelection}/{sample}_Sel_{readSelection}.fq"
    output:
        "SelectedReads/{readSelection}/{sample}_Sel_{readSelection}.fa"
    shell:
        """
        sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
        """


# gzip fastq for selected reads
rule selRead_gzip_fastq:
    message: "Compress {input}"
    input: "SelectedReads/{readSelection}/{sample}_Sel_{readSelection}.fq"
    output: "SelectedReads/{readSelection}/{sample}_Sel_{readSelection}.fq.gz"
    threads: 1
    shell: "gzip -c {input} > {output}"


# Prepare VDV reads for multiple sequence alignments (by replacing the vector sequence)
rule Prepare_VDV_reads:
    message: "Preparing {wildcards.sample} VDV reads for alignment"
    input:
        VDVfasta = "SelectedReads/VDV/{sample}_Sel_VDV.fa",
        blastvec = "align/Blast/results/{sample}_vectorblast.res"
    params:
        FilterOutputfasta = True,
        vectorSequence = config['VectorFasta'],
        RestrictionSite = config['RestrictionSite'],
        UnalignedVectorLength = 1000,
        SideSeqSearch = 10
    output:
        outVDVfasta = "SelectedReads/VDVprepared/{sample}_SelectedPreparedVDVreads.fa",
        VDVjunctions = "SelectedReads/VDVprepared/{sample}_VDVreads_InsertVectorJunctions.rds"
    log:
        "log/{sample}_VDVreadPreparation.log"
#    envmodules: "r/3.6.0"
    singularity:
        "shub://pgpmartin/SingIMG:r36_nanobac"
    threads: 4
    shell:
        """
        Rscript {SCRIPTDIR}/prepareVDVreads_args.R \
          -inputVDVfasta={input.VDVfasta} \
          -outputVDVfasta={output.outVDVfasta} \
          -outputVDVjunctions={output.VDVjunctions} \
          -FilterOutputfasta={params.FilterOutputfasta} \
          -blastvec={input.blastvec} \
          -RestrictionSite={params.RestrictionSite} \
          -VectorSequence={params.vectorSequence} \
          -UnalignedVectorLength={params.UnalignedVectorLength} \
          -SideSeqSearch={params.SideSeqSearch} \
          -Ncores={threads} \
        >> {log} 2>&1
        """

rule DefineRandomSampleParam:
    message: "Define the parameters of VDV reads random sampling for {wildcards.sample}"
    input:
        "SelectedReads/VDVprepared/{sample}_SelectedPreparedVDVreads.fa"
    output:
        NumVDVreads = "params/VDVRndSets/{sample}_NumVDVreads.txt",
        NumSeqPerSample = "params/VDVRndSets/{sample}_NumSeqPerSample.txt",
        NumRndSample = "params/VDVRndSets/{sample}_NumRndSample.txt"
    run:
        RndSetParam = VDVcount(input[0])
        with open(output.NumVDVreads, "w") as out:
            out.write(str(RndSetParam[0]))
        with open(output.NumRndSample, "w") as out:
            out.write(str(RndSetParam[1]))
        with open(output.NumSeqPerSample, "w") as out:
            out.write(str(RndSetParam[2]))


# Obtain random samples of VDV reads
checkpoint RandomSampling_VDVreads:
    message: "Obtaining random samples of prepared VDV reads for {wildcards.sample}"
    input:
        VDVreads = "SelectedReads/VDVprepared/{sample}_SelectedPreparedVDVreads.fa",
        NumVDVreads = "params/VDVRndSets/{sample}_NumVDVreads.txt",
        NumSeqPerSample = "params/VDVRndSets/{sample}_NumSeqPerSample.txt",
        NumRndSample = "params/VDVRndSets/{sample}_NumRndSample.txt"
    params:
        sampleName = "{sample}"
    output:
        directory("SelectedReads/VDVprepared/RandomSamples/{sample}")
    log:
        "log/{sample}_VDVreadsRandomSampling.log"
    conda: "envs/seqtk.yaml"
#    envmodules: "seqtk/1.3"
    threads: 1
    shell:
        """
        NUMVDV=$(cat {input.NumVDVreads})
        NUMSET=$(cat {input.NumRndSample})
        NUMSEQ=$(cat {input.NumSeqPerSample})

        printf "%s\n" "Number of prepared VDV reads: ${{NUMVDV}}" >> {log}
        printf "%s\n" "Number of random sets of VDV reads: ${{NUMSET}}" >> {log}
        printf "%s\n" "Number of VDV reads per random set: ${{NUMSEQ}}" >> {log}

        fastafile={input.VDVreads} \
        sampleName={params.sampleName} \
        outDIR={output} \
        NumRndSample=${{NUMSET}} \
        NumSeqPerSample=${{NUMSEQ}} \
        {SCRIPTDIR}/MakeRandomSamplesFromFasta.sh \
        >> {log} 2>&1
        """

# Align VDV random samples
rule Kalign_RandomSetsOfVDVreads:
    message: "Aligning random samples RS{wildcards.NUM} of VDV reads for {wildcards.sample}"
    input:
        "SelectedReads/VDVprepared/RandomSamples/{sample}/{sample}_RS{NUM}.fa"
    params:
        gapopen = 55,
        gapextension = 8.5,
        teminalgap = 4.5,
        matrixbonus = 0.2
    output:
        "align/kalign/VDVrndSets/{sample}_RS{NUM}.msf"
    log:
        "log/{sample}_kalign_RS{NUM}.log"
    conda: "envs/kalign2.yaml"
#    envmodules: "kalign/2.04"
    threads: 1
    shell:
        """
        kalign \
            -dna \
            -gapopen {params.gapopen} \
            -gpe {params.gapextension} \
            -tgpe {params.teminalgap} \
            -bonus {params.matrixbonus} \
            -infile {input} \
            -outfile {output} \
            -format msf
        >> {log} 2>&1
        """

# Obtain consensus from multiple alignment of VDV reads random samples
rule consensus_VDVrandomSets:
    message: "Consensus sequence from RS{wildcards.NUM} of VDV reads for {wildcards.sample}"
    input:
        "align/kalign/VDVrndSets/{sample}_RS{NUM}.msf"
    output:
        "align/kalign/VDVrndSets/{sample}_RS{NUM}_cons.fa"
    params:
        consName = "cons_{sample}_RS{NUM}"
    log:
        "log/{sample}_ConsensusVDVRndSet_RS{NUM}.log"
#    envmodules: "r/3.6.0"
    singularity:
        "shub://pgpmartin/SingIMG:r36_nanobac"
    threads: 1
    shell:
        """
        Rscript {SCRIPTDIR}/getNanoBACconsensus.R \
          -msfFile={input} \
          -outputfile={output} \
          -consName={params.consName} \
        >> {log} 2>&1
        """

# See https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
def ConsensusFromRndSets(wildcards):
    checkpoint_output = checkpoints.RandomSampling_VDVreads.get(**wildcards).output[0]
    return expand("align/kalign/VDVrndSets/{sample}_RS{NUM}_cons.fa",
            sample=wildcards.sample,
            NUM=glob_wildcards(os.path.join(checkpoint_output, "{sample}_RS{NUM}.fa")).NUM)

# Align consensus obtained from random sets of VDV reads
rule align_consensus_VDVrandomSets:
    message: "Align consensus from random sets of VDV reads for {wildcards.sample}"
    input:
        ConsensusFromRndSets
    params:
        gapopen = 55,
        gapextension = 8.5,
        teminalgap = 4.5,
        matrixbonus = 0.2
    output:
        "align/kalign/ConsFromVDV/{sample}_consAlign.msf"
    log:
        "log/{sample}_consAlign.log"
    conda: "envs/kalign2.yaml"
#    envmodules: "kalign/2.04"
    threads: 1
    shell:
        """
        cat {input} | \
        kalign \
            -dna \
            -gapopen {params.gapopen} \
            -gpe {params.gapextension} \
            -tgpe {params.teminalgap} \
            -bonus {params.matrixbonus} \
            -outfile {output} \
            -format msf
        >> {log} 2>&1
        """


# Obtain final consensus from multiple alignment
rule consensus_fromVDVcons:
    message: "Consensus from MSA of random sets of VDV reads for {wildcards.sample}"
    input:
        "align/kalign/ConsFromVDV/{sample}_consAlign.msf"
    output:
        "assembly/{sample}_InitialConsensus.fa"
    params:
        consName = "cons_{sample}"
    log:
        "log/{sample}_ConsensusFromVDVcons.log"
#    envmodules: "r/3.6.0"
    singularity:
        "shub://pgpmartin/SingIMG:r36_nanobac"
    threads: 1
    shell:
        """
        Rscript {SCRIPTDIR}/getNanoBACconsensus.R \
          -msfFile={input} \
          -outputfile={output} \
          -consName={params.consName} \
        >> {log} 2>&1
        """


# Select long VD reads for polishing
rule Select_longVDreads:
    message: "Selecting long VD reads for {wildcards.sample}"
    input:
        "tables/ReadClass/{sample}_ReadClass.rds"
    params:
        WithGeneA = True,
        WithGeneB = True,
        isHostAlign = False,
        minLength = config['minVDreadLength']
    output:
        outRDS = "SelectedReads/longVD/{sample}_Sel_longVD.rds",
        outReadNames = "SelectedReads/longVD/{sample}_Sel_longVD_ReadNames.tsv"
    log:
        "log/{sample}_Sel_longVD_Selection.log"
#    envmodules: "r/3.6.0"
    singularity:
        "shub://pgpmartin/SingIMG:r36_nanobac"
    threads: 1
    shell:
        """
        Rscript {SCRIPTDIR}/selectlongVDreads_args.R \
          -ReadClass={input} \
          -minLength={params.minLength} \
          -WithGeneA={params.WithGeneA} \
          -WithGeneB={params.WithGeneB} \
          -isHostAlign={params.isHostAlign} \
          -outRDS={output.outRDS} \
          -outReadNames={output.outReadNames} \
        >> {log} 2>&1
        """


# Select DVD reads (if any)
rule SelectSplit_DVDreads:
    message: "Selecting and splitting DVD reads for {wildcards.sample}"
    input:
        ReadClass = "tables/ReadClass/{sample}_ReadClass.rds",
        blastvec = "align/Blast/results/{sample}_vectorblast.res",
        FastaFile = "data/renamed/{sample}.fa"
    params:
        WithGeneA = True,
        WithGeneB = True,
        MinDNAlength = config['minDNAlength_DVDreads']
    output:
        outRDS = "SelectedReads/DVD/{sample}_SelSplit_DVD.rds",
        outReadNames = "SelectedReads/DVD/{sample}_SelSplit_DVD_ReadNames.tsv",
        outFasta = "SelectedReads/DVD/{sample}_SelSplit_DVD.fa"
    log:
        "log/{sample}_DVD_SelectAndSplit.log"
#    envmodules: "r/3.6.0"
    singularity:
        "shub://pgpmartin/SingIMG:r36_nanobac"
    threads: 1
    shell:
        """
        Rscript {SCRIPTDIR}/selectSplitDVDreads_args.R \
          -ReadClass={input.ReadClass} \
          -blastvec={input.blastvec} \
          -FastaFile={input.FastaFile} \
          -MinDNAlength={params.MinDNAlength} \
          -WithGeneA={params.WithGeneA} \
          -WithGeneB={params.WithGeneB} \
          -outRDS={output.outRDS} \
          -outReadNames={output.outReadNames} \
          -outFasta={output.outFasta} \
        >> {log} 2>&1
        """

# Merge DV/VD reads (from DVD reads) with long VD reads to use for polishing step
rule merge_VD_DV_reads:
    message: "Merging fasta files for VD/DV reads (used for polishing) for {wildcards.sample}"
    input:
        DVD = "SelectedReads/DVD/{sample}_SelSplit_DVD.fa",
        VD = "SelectedReads/longVD/{sample}_Sel_longVD.fa"
    output:
        "SelectedReads/{sample}_VDreadsForPolishing.fa"
    log:
        "log/{sample}_Sel_VDreadsForPolishing.log"
    shell:
        """
        cat {input.DVD} {input.VD} > {output}
        NumReads=$(grep "^>" {output} | wc -l)
        if (( NumReads == 0 )); then
            printf "%s\\n" "Zero reads for polishing step" > {log}
            echo "Zero reads for polishing step. Stopping Here!"
            exit 1
        else
            printf "%s\\n" "Number of reads for polishing step: ${{NumReads}}" > {log}
        fi
        """

# Align VDreads to consensus
## TODO: evaluate if winnowmap (https://www.biorxiv.org/content/10.1101/2020.02.11.943241v1) instead of minimap2 gives better results at this step
## TODO: check if/how minimap2 uses the base-level quality scores. If it makes a difference, input fastq instead of fasta
rule align_VD_to_cons:
    message: "Aligning VD reads to consensus for {wildcards.sample}"
    input:
        VDreads = "SelectedReads/{sample}_VDreadsForPolishing.fa",
        consensus = "assembly/{sample}_{consensusType}.fa"
    output:
        "align/minimap2/polishing_{consensusType}/{sample}.paf"
    log:
        "log/{sample}_polishing_{consensusType}.log"
    conda: "envs/minimap2.yaml"
#    envmodules: "minimap2/2.17"
    threads: 4
    shell:
        """
        minimap2 \
          -x map-ont \
          -t {threads} \
          {input.consensus} \
          {input.VDreads} > \
          {output}
        """

# Polish using racon
rule racon_Polish:
    message: "Polishing consensus for {wildcards.sample}_{wildcards.consensusType}"
    input:
        VDreads = "SelectedReads/{sample}_VDreadsForPolishing.fa",
        consensus = "assembly/{sample}_{consensusType}.fa",
        VDalign = "align/minimap2/polishing_{consensusType}/{sample}.paf"
    output:
        "assembly/polishing/{sample}_{consensusType}_polished.fa"
    conda: "envs/racon.yaml"
#    envmodules: "racon/1.4"
    threads: 4
    shell:
        """
        racon \
            -t {threads} \
            {input.VDreads} \
            {input.VDalign} \
            {input.consensus} > \
            {output}
        """

# Define polishing step 01
rule Polishing_step01:
    message: "moving {input}"
    input: "assembly/polishing/{sample}_InitialConsensus_polished.fa"
    output: "assembly/{sample}_Polished01.fa"
    shell: "mv {input} {output}"

# Define polishing step 02
# To stop the polishing loop, do not name the final assembly {sample}_{something}.fa
rule Polishing_step02:
    message: "moving {input}"
    input: "assembly/polishing/{sample}_Polished01_polished.fa"
    output: "assembly/final_assembly_{sample}.fa"
    shell:
        """
        mv {input} {output}
        if ls -1qA assembly/polishing/ | grep -q .
        then
          echo "assembly/polishing directory is not empty"
        else
          rmdir assembly/polishing
        fi
        """
#Done
