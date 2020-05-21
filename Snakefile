configfile:
    "config.json"

shell.executable("/bin/bash")

# SAMPLENAME = config['sampleName']
WORKDIR = config['workingDIR']
LOGDIR = WORKDIR+"/log"
SCRIPTDIR = config['scriptDIR']
SAMPLENAME, = glob_wildcards(WORKDIR+"/data/raw/{id}.fastq")

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
        expand("align/minimap2/host/{sample}_alignment.stats", sample=SAMPLENAME)


# Create a table with new names and old names
rule Read_NameMapping:
    message: "Storing old and new names for {input} in {output}."
    input:
        "data/raw/{sample}.fastq"
    output:
        "tables/ReadNames/{sample}_ReadNameTable.tsv"
    params:
        sampleName = "{sample}"
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
    input:
        "data/renamed/{sample}.fq"
    output:
        "data/renamed/{sample}.fq.gz"
    shell:
        """
        gzip -c {input} > {output}
        """


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
          -num_alignments 100000000
        """

# Check for the presence of data for geneA
# rule check_GeneA_exists:
#         message: "Check that GeneA is used"
#         input:
#             geneA=config['GeneAFasta']
#         shell:
#             """
#             if [[ ! -z {geneA} ]]
#             then
#                 touch .GeneAPresent
#             else
#                 touch .GeneAAbsent # maybe touch align/Blast/results/{sample}_geneAblast.res
#             """

# Align to gene A
rule AlignGeneA:
    message: "Aligning GeneA on the reads {input}"
    input:
        query = config['GeneAFasta'],
        db = "align/Blast/db/{sample}_blastdb.nhr"
    output:
        "align/Blast/results/{sample}_geneAblast.res"
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
          -num_alignments 100000000
        """

# Align to gene B
rule AlignGeneB:
    message: "Aligning GeneB on the reads {input}"
    input:
        query = config['GeneBFasta'],
        db = "align/Blast/db/{sample}_blastdb.nhr"
    output:
        "align/Blast/results/{sample}_geneBblast.res"
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
          -num_alignments 100000000
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
    threads: 4
    shell:
        """
        minimap2 \
          -ax map-ont \
          -t {threads} \
          -L \
          {input.hostGenome} \
          {input.reads} > \
          {output} \
          2> {log}
        """

# Convert to paf (only convert the primary and supplementary alignments with the -p argument )
rule Hostsam2paf:
    message: "Convert SAM to PAF ({input})"
    input:
        "align/minimap2/host/{sample}.sam"
    output:
        "align/minimap2/host/{sample}.paf"
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
    shell:
        """
        samtools view \
          -bh \
          -F 260 \
          {input} | \
        samtools sort \
          -@ {threads} \
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
