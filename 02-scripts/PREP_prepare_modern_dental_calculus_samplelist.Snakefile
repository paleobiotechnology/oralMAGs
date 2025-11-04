################################################################################
# Project: Oral MAGs
# Part: Data preparation
# Step: Compile the sample list of modern dental samples from the study
#       Fellows Yates 2021
#
# Dependent on:
#
# Alex Hübner
################################################################################

import os

import pandas as pd

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
ENA_PROJECT_CODES = {
    'JAE': 'PRJEB31185',
    'VLC': 'PRJEB34569',
}
################################################################################

rule all:
    input:
        "01-resources/modern_dental_calculus_samples.tsv"

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
        "pyscripts/PREP_prepare_modern_dental_calculus_samplelist-pull_seqdata_info.py"

rule concat:
    input:
        expand("scratch_tmp/oralmags/{study}_samplelist.tsv", study=ENA_PROJECT_CODES)
    output:
        "01-resources/modern_dental_calculus_samples.tsv"
    message: "Concatenate the sequencing data information from the two ENA project codes"
    resources:
        mem_mb = 4000
    threads: 1
    run:
        (pd.concat([pd.read_csv(fn, sep="\t")
                    for fn in input])
         .sort_values(['publication_year', 'secondary_sample_accession'])
         .to_csv(output[0], sep="\t", index=False)
        )
