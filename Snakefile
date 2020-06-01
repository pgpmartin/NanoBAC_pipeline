configfile:
    "config.json"

shell.executable("/bin/bash")

localrules: all, Rename_Reads, RenamedReads_fq2fa, ReadLengthTable, AlignVector, AlignGeneA,
            AlignGeneB, Annotate_Reads, Select_VDVreads, selRead_Get_fastq, selRead_fastq2fasta, Prepare_VDV_reads,
            RandomSampling_VDVreads, consensus_VDVrandomSets, consensus_fromVDVcons, Select_longVDreads, SelectSplit_DVDreads

WORKDIR = config['workingDIR']
LOGDIR = WORKDIR+"/log"
SCRIPTDIR = config['scriptDIR']
SAMPLENAME, = glob_wildcards(WORKDIR+"/data/raw/{id}.fastq")
NUMRNDSET = config['VDVreadsNumberOfRandomSets']
RNDSETS = list(range(1, int(NUMRNDSET)+1))
#modules
#import os.path

rule all:
    input:
        expand("tables/ReadNames/{sample}_ReadNameTable.tsv", sample=SAMPLENAME),
        expand("data/renamed/{sample}.fq.gz", sample=SAMPLENAME),
        expand("data/renamed/{sample}.fa", sample=SAMPLENAME),
        expand("tables/ReadLength/{sample}_ReadLength.tsv", sample=SAMPLENAME),
        expand("tables/ReadClass/{sample}_ReadClass.rds", sample=SAMPLENAME),
        expand("align/Blast/results/{sample}_vectorblast.res", sample=SAMPLENAME),
        expand("align/Blast/results/{sample}_geneAblast.res", sample=SAMPLENAME),
        expand("align/Blast/results/{sample}_geneBblast.res", sample=SAMPLENAME),
        expand("align/minimap2/host/{sample}.paf", sample=SAMPLENAME),
        expand("align/minimap2/host/{sample}_alignment.stats", sample=SAMPLENAME),
        expand("SelectedReads/VDV/{sample}_Sel_VDV.rds", sample=SAMPLENAME),
        expand("SelectedReads/VDV/{sample}_Sel_VDV.fa", sample=SAMPLENAME),
        expand("SelectedReads/VDV/{sample}_Sel_VDV.fq.gz", sample=SAMPLENAME),
        expand("SelectedReads/VDVprepared/{sample}_SelectedPreparedVDVreads.fa", sample=SAMPLENAME),
        expand("SelectedReads/VDVprepared/RandomSamples/{sample}_RS{NUM}.fa", sample=SAMPLENAME, NUM = RNDSETS),
        expand("align/kalign/VDVrndSets/{sample}_RS{NUM}.msf", sample=SAMPLENAME, NUM=RNDSETS),
        expand("align/kalign/VDVrndSets/{sample}_RS{NUM}_cons.fa", sample=SAMPLENAME, NUM=RNDSETS),
        expand("align/kalign/ConsFromVDV/{sample}_consAlign.msf", sample=SAMPLENAME),
        expand("assembly/{sample}_InitialConsensus.fa", sample=SAMPLENAME),
        expand("SelectedReads/longVD/{sample}_Sel_longVD.rds", sample=SAMPLENAME),
        expand("SelectedReads/longVD/{sample}_Sel_longVD.fa", sample=SAMPLENAME),
        expand("SelectedReads/longVD/{sample}_Sel_longVD.fq.gz", sample=SAMPLENAME),
        expand("SelectedReads/DVD/{sample}_SelSplit_DVD.fa", sample=SAMPLENAME)


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
    message: "Renaming reads for {input}"
    input:
        "data/raw/{sample}.fastq"
    params:
        sampleName = "{sample}"
    output:
        "data/renamed/{sample}.fq"
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
    message: "Extracting Read Length"
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
rule CreateBlastDatabase:
    message: "Creating BLAST database for {input}"
    input:
        "data/renamed/{sample}.fa"
    output:
        "align/Blast/db/{sample}_blastdb.nhr"
    log:
        "log/{sample}_makeblastdb.log"
    threads: 1
    conda: "envs/blast.yaml"
    shell:
        """
        outbase={output}
        outbase=${{outbase%%.nhr}}
        makeblastdb \
          -in {input} \
          -out ${{outbase}} \
          -parse_seqids \
          -dbtype nucl \
          -logfile {log}
        """

# Align the vector
rule AlignVector:
    message: "Aligning vector sequence on the reads {input}"
    input:
        query = config['VectorFasta'],
        db = "align/Blast/db/{sample}_blastdb.nhr"
    output:
        "align/Blast/results/{sample}_vectorblast.res"
    log:
        "log/{sample}_blastVector.log"
    conda: "envs/blast.yaml"
    threads: 1
    shell:
        """
        dbbase={input.db}
        dbbase=${{dbbase%%.nhr}}
        blastn \
          -task megablast \
          -query {input.query} \
          -db ${{dbbase}} \
          -out {output} \
          -outfmt 6 \
          -num_alignments 100000000 \
        >> {log} 2>&1
        """

# Align to gene A
rule AlignGeneA:
    message: "Aligning GeneA on the reads {input}"
    input:
        query = config['GeneAFasta'],
        db = "align/Blast/db/{sample}_blastdb.nhr"
    output:
        "align/Blast/results/{sample}_geneAblast.res"
    log:
        "log/{sample}_blastGeneA.log"
    conda: "envs/blast.yaml"
    threads: 1
    shell:
        """
        dbbase={input.db}
        dbbase=${{dbbase%%.nhr}}
        blastn \
          -task megablast \
          -query {input.query} \
          -db ${{dbbase}} \
          -out {output} \
          -outfmt 6 \
          -num_alignments 100000000 \
        >> {log} 2>&1
        """

# Align to gene B
rule AlignGeneB:
    message: "Aligning GeneB on the reads {input}"
    input:
        query = config['GeneBFasta'],
        db = "align/Blast/db/{sample}_blastdb.nhr"
    output:
        "align/Blast/results/{sample}_geneBblast.res"
    log:
        "log/{sample}_blastGeneB.log"
    conda: "envs/blast.yaml"
    threads: 1
    shell:
        """
        dbbase={input.db}
        dbbase=${{dbbase%%.nhr}}
        blastn \
          -task megablast \
          -query {input.query} \
          -db ${{dbbase}} \
          -out {output} \
          -outfmt 6 \
          -num_alignments 100000000 \
        >> {log} 2>&1
        """

# Align to host genome
rule Align2Host:
    message: "Aligning reads to host genome ({input})"
    input:
        reads = "data/renamed/{sample}.fq",
        hostGenome = config['hostGenomeFasta']
    output:
        temp("align/minimap2/host/{sample}.sam")
    log:
        "log/{sample}_hostalign.log"
    conda: "envs/minimap2.yaml"
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
    message: "Convert SAM to PAF ({input})"
    input:
        "align/minimap2/host/{sample}.sam"
    output:
        "align/minimap2/host/{sample}.paf"
    threads: 1
    conda: "envs/minimap2.yaml"
    shell:
        """
        paftools.js sam2paf -p {input} > {output}
        """

# Gather some stats
rule HostAlignStats:
    message: "Host alignment stats"
    input:
        "align/minimap2/host/{sample}.sam"
    params:
        sampleName = "{sample}"
    output:
        "align/minimap2/host/{sample}_alignment.stats"
    threads: 1
    conda: "envs/samtools.yaml"
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
    message: "Annotating reads"
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
    message: "Selecting VDV reads"
    input:
        "tables/ReadClass/{sample}_ReadClass.rds"
    params:
        SizeTolerance = 0.05,
        WithGeneA = True,
        WithGeneB = True,
        MaxClusters = 10,
        plotVar = "InsertLength",
        plotRes = 150,
        plotHeight = 4,
        plotWidth = 4,
        plotFormat = "png"
    output:
        outPlot = "Plots/{sample}_VDVreadSelection.png",
        outRDS = "SelectedReads/VDV/{sample}_Sel_VDV.rds",
        outVDVnames = "SelectedReads/VDV/{sample}_Sel_VDV_ReadNames.tsv"
    log:
        "log/{sample}_Sel_VDV_Selection.log"
    threads: 1
    shell:
        """
        Rscript {SCRIPTDIR}/selectVDVreads_args.R \
          -ReadClass={input} \
          -SizeTolerance={params.SizeTolerance} \
          -WithGeneA={params.WithGeneA} \
          -WithGeneB={params.WithGeneB} \
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
    message: "Obtain fastq file of selected reads"
    input:
        rawfastq = "data/renamed/{sample}.fq",
        ReadSelection = "SelectedReads/{readSelection}/{sample}_Sel_{readSelection}_ReadNames.tsv"
    output:
        "SelectedReads/{readSelection}/{sample}_Sel_{readSelection}.fq"
    conda: "envs/seqtk.yaml"
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
    message: "Convert fastq to fasta for selected reads"
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
    message: "Compress fastq for selected reads"
    input: "SelectedReads/{readSelection}/{sample}_Sel_{readSelection}.fq"
    output: "SelectedReads/{readSelection}/{sample}_Sel_{readSelection}.fq.gz"
    threads: 1
    shell: "gzip -c {input} > {output}"


# Prepare VDV reads for multiple sequence alignments (by replacing the vector sequence)
rule Prepare_VDV_reads:
    message: "Preparing VDV reads for alignment"
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

# Obtain random samples of VDV reads
rule RandomSampling_VDVreads:
    message: "Obtaining random samples of prepared VDV reads"
    input:
        "SelectedReads/VDVprepared/{sample}_SelectedPreparedVDVreads.fa"
    params:
        sampleName = "{sample}",
        NumSeqPerSample = config['VDVreadsNumberOfSequencesPerRndSample'],
        NumRndSample = NUMRNDSET
    output:
        expand("SelectedReads/VDVprepared/RandomSamples/{{sample}}_RS{NUM}.fa", NUM = RNDSETS)
    log:
        "log/{sample}_VDVreadsRandomSampling.log"
    conda: "envs/seqtk.yaml"
    threads: 1
    shell:
        """
        fastafile={input} \
        sampleName={params.sampleName} \
        outDIR=$(dirname {output[0]}) \
        NumRndSample={params.NumRndSample} \
        NumSeqPerSample={params.NumSeqPerSample} \
        {SCRIPTDIR}/MakeRandomSamplesFromFasta.sh \
        >> {log} 2>&1
        """
## Another option for this rule is to use dynamic files. See https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
## it would allow to not throw an error at less than 20 VDV reads (not sure it's a good thing though??)
## However that may make downstream rules more complicated (not knowing the number of input files)
## dynamic-flag will be deprecated in Snakemake 6.0
## It is thus recommended to use checkpoints instead: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution

# Align VDV random samples
rule Kalign_RandomSetsOfVDVreads:
    message: "Aligning random samples of VDV reads"
    input:
        "SelectedReads/VDVprepared/RandomSamples/{sample}_RS{NUM}.fa"
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

# Obtain consensus from alignment of VDV reads random samples
rule consensus_VDVrandomSets:
    message: "Obtaining consensus sequence from multiple alignment of random sets of VDV reads"
    input:
        "align/kalign/VDVrndSets/{sample}_RS{NUM}.msf"
    output:
        "align/kalign/VDVrndSets/{sample}_RS{NUM}_cons.fa"
    params:
        consName = "cons_{sample}_RS{NUM}"
    log:
        "log/{sample}_ConsensusVDVRndSet_RS{NUM}.log"
    threads: 1
    shell:
        """
        Rscript {SCRIPTDIR}/getNanoBACconsensus.R \
          -msfFile={input} \
          -outputfile={output} \
          -consName={params.consName} \
        >> {log} 2>&1
        """

# Align consensus obtained from random sets of VDV reads
rule align_consensus_VDVrandomSets:
    message: "Align consensus from multiple alignment of random sets of VDV reads"
    input:
        expand("align/kalign/VDVrndSets/{{sample}}_RS{NUM}_cons.fa", NUM=RNDSETS)
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
    input:
        "align/kalign/ConsFromVDV/{sample}_consAlign.msf"
    output:
        "assembly/{sample}_InitialConsensus.fa"
    params:
        consName = "cons_{sample}"
    log:
        "log/{sample}_ConsensusFromVDVcons.log"
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
    message: "Selecting long VD reads"
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
    message: "Selecting and splitting DVD reads"
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

# Merge DVD reads with long VD reads

# Polish consensus round 1

# Polish consensus round 2
    #final assembly/
