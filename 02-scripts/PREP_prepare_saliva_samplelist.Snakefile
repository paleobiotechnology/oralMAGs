################################################################################
# Project: Oral MAGs
# Part: Data preparation
# Step: Compile the sample list of published saliva samples from the studies
#       Clemente 2015, Lassalle 2018, and Pedro 2022
#
# Dependent on:
#
# Alex Hübner
################################################################################

import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
ENA_PROJECT_CODES = {
    'Clemente2015': 'PRJNA245336',
    'Lassalle2018': 'PRJEB14383',
    'Pedro2022': 'PRJEB54966',
}
################################################################################

rule all:
    input:
        "01-resources/modern_saliva_samples.tsv"

rule pull_seqdata_info:
    output:
        temp("scratch_tmp/oralmags/{study}_samplelist.tsv")
    message: "Pull the information on the sequencing data from ENA: {wildcards.study}"
    conda: "ENVS_amdirt.yaml"
    resources:
        mem_mb = 4000
    params:
        accession = lambda wildcards: ENA_PROJECT_CODES[wildcards.study]
    threads: 1
    script:
        "pyscripts/PREP_prepare_saliva_samplelist-pull_seqdata_info.py"

rule concat:
    input:
        expand("scratch_tmp/oralmags/{study}_samplelist.tsv", study=ENA_PROJECT_CODES)
    output:
        "01-resources/modern_saliva_samples.tsv"
    message: "Concatenate the sequencing data information from the three studies"
    resources:
        mem_mb = 4000
    threads: 1
    run:
        (pd.concat([pd.read_csv(fn, sep="\t")
                    for fn in input])
         .sort_values(['publication_year', 'secondary_sample_accession'])
         .to_csv(output[0], sep="\t", index=False)
        )

