configfile:
    "config.json"

# SAMPLENAME = config['sampleName']
WORKDIR = config['workingDIR']
LOGDIR = WORKDIR+"/log"
SAMPLENAME, = glob_wildcards(WORKDIR+"/data/raw/{id}.fastq")

#modules
#import os.path

rule all:
    input:
        expand("align/Blast/results/{sample}_vectorblast.res", sample=SAMPLENAME),
        expand("tables/ReadNames/ReadNameTable_{sample}.tsv", sample=SAMPLENAME),
        expand("tables/ReadLength/ReadLength_{sample}.tsv", sample=SAMPLENAME)


# Create a table with new names and old names
rule Read_NameMapping:
    message: "Storing old and new names for {input} in {output}."
    input:
        "data/raw/{sample}.fastq",
    output:
        "tables/ReadNames/ReadNameTable_{sample}.tsv"
    params:
        sampleName = config['sampleName']
    shell:
        """
        cat {input} | \
          sed -n '1~4p' | \
          awk 'BEGIN {{ cntr = 0 }} {{ cntr++ ; print \"{SAMPLENAME}R\"cntr, \"\\t\", $0 }}'  > \
          {output}
        """

# Create a fastq file with the new (simple) names
rule Rename_Reads:
    message: "Renaming reads for {input}"
    input:
        "data/raw/{sample}.fastq",
    output:
        "data/renamed/{sample}.fq"
    shell:
        """
        cat {input} | \
          sed -n '1~4s/^@.*$/@zyxabcdefzyx/p;2~4p;3~4p;0~4p' | \
          awk 'BEGIN {{ cntr = 0 }} /^@zyxabcdefzyx$/ {{ cntr++ ; print \"@{SAMPLENAME}R\"cntr }} !/^@zyxabcdefzyx$/ {{ print $0 }}' > \
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

# Get a table with read lengths
rule ReadLengthTable:
    message: "Extracting Read Length"
    input:
        "data/renamed/{sample}.fa"
    output:
        "tables/ReadLength/ReadLength_{sample}.tsv"
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
        shell:
            """
            mkdir -p {LOGDIR}
            outbase={output}
            outbase=${{outbase%%.nhr}}
            makeblastdb \
              -in {input} \
              -out ${{outbase}} \
              -parse_seqids \
              -dbtype nucl \
              -logfile {LOGDIR}/makeblastdb_{SAMPLENAME}.log
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
#                 touch .GeneAAbsent
#             """

# Classify the reads
