import os
import pandas as pd

RAWDIR = "data/raw_fastq"
RAWEXT = "fastq.gz"
RESULTSDIR = config["data/processed"]
QCDIR = os.path.join(RESULTSDIR, "01_QC")

# Sample dataframe
data_input = pd.read_csv("data/SraRunTable.csv")
# Define sample names
SAMPLES = data_input.Run.values

# What are mates called in this workflow?
MATES = ['1', '2']

rule all:
    input:
        expand(RAWDIR + "/{sample}_{mate}."+RAWEXT, sample = SAMPLES, mate = MATES)

rule download_srr:
    output: 
        r1=RAWDIR + "/{sample}_1."+RAWEXT,
        r2=RAWDIR + "/{sample}_2."+RAWEXT
    params:
        outdir=RAWDIR
    conda:
        "../envs/sra_tools.yaml"
    resources:
        mem_mb=10000,
        time_min=120,
    threads: 6
    shell:
        """
        fasterq-dump --skip-technical --split-3 {wildcards.sample} -O {params.outdir} -e {threads}
        pigz {params.outdir}/{wildcards.sample}_*
        """

rule fastp:
    conda: "envs/qc.yaml"
    resources:
        mem_mb=20000,
        time_min=360,
    threads: 12
    input:
        R1=QCDIR+"/raw_fastq/{sample}_1.fastq.gz",
        R2=QCDIR+"/raw_fastq/{sample}_2.fastq.gz",
    output:
        T1=QCDIR+"/qc_reads/{sample}_1.trim.fastq.gz",
        T2=QCDIR+"/qc_reads/{sample}_2.trim.fastq.gz",
        json=QCDIR+"/qc_reads/{sample}.fastp.json",
        html=QCDIR+"/qc_reads/{sample}.fastp.html",
    shell:
        """
        fastp \
            --in1 {input.R1} \
            --in2 {input.R2} \
            --out1 {output.T1} \
            --out2 {output.T2} \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
            --detect_adapter_for_pe
        """




