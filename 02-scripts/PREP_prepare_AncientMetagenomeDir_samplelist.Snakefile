################################################################################
# Project: Oral MAGs
# Part: Data preparation
# Step: Compile the sample list of dental calculus samples from
#       AncientMetagenomeDir from the release v23.06.0
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
AMDIR_SAMPLES_TSV = "https://raw.githubusercontent.com/SPAAM-community/AncientMetagenomeDir/v23.06.0/ancientmetagenome-hostassociated/samples/ancientmetagenome-hostassociated_samples.tsv"
AMDIR_LIBRARIES_TSV = "https://raw.githubusercontent.com/SPAAM-community/AncientMetagenomeDir/v23.06.0/ancientmetagenome-hostassociated/libraries/ancientmetagenome-hostassociated_libraries.tsv"
################################################################################

rule all:
    input:
        "01-resources/AncientMetagenomeDir_v23.06.0_dentalcalculus_samples.tsv"

rule amdir:
    output:
        "01-resources/AncientMetagenomeDir_v23.06.0_dentalcalculus_samples.tsv"
    message: "Extract all metagenomic dental calculus samples from AncientMetagenomeDir v23.06.0 that are sufficiently deeply sequenced for metagenome assembly"
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
            .query("material == 'dental calculus'") \
            .query("read_count >= 5000000 and library_strategy == 'WGS' and archive_project.str.startswith('PRJ')")
            .to_csv(output[0], sep="\t", index=False, float_format="%.3f")
        )
