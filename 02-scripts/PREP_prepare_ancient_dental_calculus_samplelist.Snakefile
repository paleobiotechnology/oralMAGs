################################################################################
# Project: Oral MAGs
# Part: Data preparation
# Step: Compile the list of dental calculus samples from the Warinner lab that
#       were not present in v23.06.0 of AncientMetagenomeDir and merge with the
#       published samples
#
# Dependent on:
#   - PREP_prepare_AncientMetagenomeDir_samplelist.Snakefile
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
# Links to AMDir release containing the samples of Velsko et al. (2024)
AMDIR_SAMPLES_TSV = "https://raw.githubusercontent.com/SPAAM-community/AncientMetagenomeDir/v25.03.0/ancientmetagenome-hostassociated/samples/ancientmetagenome-hostassociated_samples.tsv"
AMDIR_LIBRARIES_TSV = "https://raw.githubusercontent.com/SPAAM-community/AncientMetagenomeDir/v25.03.0/ancientmetagenome-hostassociated/libraries/ancientmetagenome-hostassociated_libraries.tsv"
################################################################################

rule all:
    input:
        "01-resources/ancient_dental_calculus_samples.tsv"

rule amdir_velsko2024:
    output:
        "scratch_tmp/oralmags/AncientMetagenomeDir_v25.03.0_Velsko2024.tsv"
    message: "Extract all metagenomic dental calculus samples from Velsko et al. (2024) from AncientMetagenomeDir v25.03.0"
    resources:
        mem_mb = 4000
    threads: 1
    run:
        # Prepare sample table by taking care of multiple accession codes per sample
        samples = pd.read_csv(AMDIR_SAMPLES_TSV, sep="\t")
        samples['archive_accession'] = samples['archive_accession'].str.split(",")
        samples = (samples.drop(['archive_accession'], axis=1) \
            .merge(samples[['project_name', 'sample_name', 'archive_accession']]
                   .set_index(['project_name', 'sample_name'])
                   .explode('archive_accession')
                   .reset_index(),
                   how="left", on=['project_name', 'sample_name'])
            .rename({'archive_accession': 'archive_sample_accession'}, axis=1)
        )

        # Pull the libraries table
        libraries = pd.read_csv(AMDIR_LIBRARIES_TSV,
                                sep="\t",
                                usecols=['archive_project', 'archive_sample_accession',
                                         'library_name', 'strand_type', 'library_treatment',
                                         'instrument_model', 'library_layout',
                                         'library_strategy', 'read_count',
                                         'archive_data_accession', 'download_links',
                                         'download_md5s'])
    
        # Filter for dental calculus samples with 
        (samples.merge(libraries, how="left", on=["archive_sample_accession", "archive_project"]) \
            .query("project_name == 'Velsko2024'") \
            .to_csv(output[0], sep="\t", index=False, float_format="%.3f")
        )

rule join:
    input:
        "scratch_tmp/oralmags/AncientMetagenomeDir_v25.03.0_Velsko2024.tsv"
    output:
        "01-resources/ancient_dental_calculus_samples.tsv"
    message: "Generate a single table with all ancient dental calculus samples"
    resources:
        mem_mb = 4000
    params:
        amdir = "01-resources/AncientMetagenomeDir_v23.06.0_dentalcalculus_samples.tsv",
        warinnerlab = "01-resources/WarinnerLab_ancient_dentalcalculus_samples.tsv"
    threads: 1
    run:
        velsko2024 = pd.read_csv(input[0], sep="\t")
        amdir = pd.read_csv(params.amdir, sep="\t")
        warinnerlab = pd.read_csv(params.warinnerlab, sep="\t")

        (pd.concat([amdir, velsko2024, warinnerlab])
         .sort_values(['publication_year', 'archive_sample_accession'])
         .to_csv(output[0], sep="\t", index=False)
        )
