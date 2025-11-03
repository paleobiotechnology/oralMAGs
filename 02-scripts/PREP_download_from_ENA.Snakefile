################################################################################
# Project: Oral MAGs
# Part: Data preparation
# Step: Download the sequencing data from ENA
#
# Dependent on:
#   - PREP_prepare_AncientMetagenomeDir_samplelist.Snakefile
#   - PREP_prepare_ancient_dental_calculus_samplelist.Snakefile 
#
# Alex Hübner
################################################################################

import json
import os
from pathlib import Path
import sys

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

#### Load sample data ##########################################################
ancdentalcalc = (pd.read_csv("01-resources/ancient_dental_calculus_samples.tsv", sep="\t",
                             usecols=['archive_data_accession', 'library_layout'])
                 .rename({'archive_data_accession': 'run_accession'}, axis=1)
                 .set_index(['run_accession'])
)
saliva = (pd.read_csv("01-resources/modern_saliva_samples.tsv", sep="\t",
                      usecols=['run_accession', 'library_layout'])
          .set_index(['run_accession'])
)
samples = ancdentalcalc.index.tolist() + saliva.index.tolist()
################################################################################

#### Auxilliary functions ######################################################
EXPECTED_SUFFICES = {1: ['0'],
                     2: ['1', '2'],
                     3: ['0', '1', '2']}

def return_url(wildcards):
    with open(checkpoints.fetch_ena_info.get(**wildcards).output[0], "rt") as jsonfile:
        ena_rep = json.load(jsonfile)
    return ena_rep[int(wildcards.i) - 1]['url']


def return_expected_suffices(wildcards):
    with open(checkpoints.fetch_ena_info.get(**wildcards).output[0], "rt") as jsonfile:
        ena_rep = json.load(jsonfile)
    return EXPECTED_SUFFICES[len(ena_rep)]

################################################################################

rule all:
    input:
        expand("03-data/raw_data/{err}.validated", err=['ERR10167375']) # samples)

checkpoint fetch_ena_info:
    output:
        "tmp/ffq/{err}.json"
    message: "Download the FTP URL and MD5sum from ENA: {wildcards.err}"
    conda: "ENVS_ffq.yaml"
    resources:
        mem_mb = 4000
    shell:
        "ffq --ftp {wildcards.err} > {output}"

rule download_fastq_file:
    input:
        "tmp/ffq/{err}.json"
    output:
        "03-data/raw_data/{err}_{i}.fastq.gz"
    message: "Download the FastQ file: {wildcards.err} for read {wildcards.i}"
    resources:
        mem_mb = 4000
    params:
        url = lambda wildcards: return_url(wildcards),
        cutdirs = lambda wildcards: return_url(wildcards).count("/") - 1
    shell:
        """
        wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent \
            --reject "index.html*" --cut-dirs={params.cutdirs} -O {output} {params.url}
        """

rule calculate_md5sum:
    input:
        "03-data/raw_data/{err}_{i}.fastq.gz"
    output:
        temp("03-data/raw_data/{err}_{i}.md5")
    message: "Calculate the md5sum: {wildcards.err} for read {wildcards.i}"
    resources:
        mem_mb = 4000
    shell:
        "md5sum {input} > {output}"

rule validate_md5sum:
    input:
        lambda wildcards: [f"03-data/raw_data/{wildcards.err}_{i}.md5" for i in return_expected_suffices(wildcards)]
    output:
        "03-data/raw_data/{err}.validated"
    message: "Validate md5sum: {wildcards.err}"
    resources:
        mem_mb = 4000
    run:
        with open(checkpoints.fetch_ena_info.get(**wildcards).output[0], "rt") as jsonfile:
            ena_rep = json.load(jsonfile)

        md5sums = [open(fn, 'rt').readline().split()[0]
                   for fn in input]

        if all([ena_rep[i]['md5'] == md5sums[i] for i in range(len(md5sums))]):
            Path(output[0]).touch()
        else:
            print("The md5sums don't match. The files have the md5sums\n" +
                  ", ".join(md5sums) +
                  "\nwhile ENA reports the following values\n" +
                  ", ".join([f['md5'] for f in ena_rep]))
            sys.exit(1)
