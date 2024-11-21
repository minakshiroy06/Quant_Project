import os
import pandas as pd

RAWDIR = "data/raw_fastq"
RAWEXT = "fastq.gz"

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






