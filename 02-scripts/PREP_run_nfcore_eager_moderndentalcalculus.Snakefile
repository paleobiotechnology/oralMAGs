################################################################################
# Project: Oral MAGs
# Part: Data preparation
# Step: Run nf-core/eager to process the modern dental calculus sequencing data
#
# nf-core/eager performs data preprocessing, e.g. adapter trimming, and
# alignment against the human reference genome hs37ds.
#
# Dependent on:
#   - PREP_download_from_ENA.Snakefile
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
samples = pd.read_csv("01-resources/modern_dental_calculus_samples.tsv", sep="\t")
################################################################################

#### Auxilliary functions ######################################################
chemistry = {
    2: ['NextSeq 500', 'Illumina NovaSeq 6000', 'HiSeq X Ten',
        'Illumina NextSeq 500', 'BGISEQ-500'],
    4: ['Illumina HiSeq 2500', 'Illumina HiSeq 4000',
        'Illumina HiSeq2000', 'Illumina HiSeq 3000',
        'Illumina HiSeq 2000', 'Illumina MiSeq'],
}
instrument_model = {sequencer: chem for chem in chemistry for sequencer in chemistry[chem]}
################################################################################

rule all:
    input:
        "04-analysis/eager/moderndentalcalc.done"

rule generate_eager_tsv:
    output:
        "04-analysis/eager/moderndentalcalc.tsv",
    message: "Generate the EAGER input TSV for the modern dental calculus samples"
    resources:
        mem_mb = 4000
    threads: 1
    run:
        eager = (samples[['run_accession', 'secondary_sample_accession',
                          'library_layout', 'instrument_model']]
            .rename({'run_accession': 'Library_ID',
                     'secondary_sample_accession': 'Sample_Name',
                     'library_layout': 'SeqType'},
                    axis=1)
        )

        eager['Lane'] = 1
        eager['Colour_Chemistry'] = [instrument_model[im]
                                     for im in eager['instrument_model'].values]
        eager['SeqType'] = ["PE" if st == "PAIRED" else "SE"
                            for st in eager['SeqType'].values]
        eager['Organism'] = "Homo sapiens"
        eager['Strandedness'] = "dsDNA"
        eager['UDG_Treatment'] = "none"
        eager['R1'] = "03-data/raw_data/" + eager['Library_ID'] + "_1.fastq.gz"
        eager['R2'] = "03-data/raw_data/" + eager['Library_ID'] + "_2.fastq.gz"
        eager['BAM'] = "NA"

        (eager[['Sample_Name', 'Library_ID', 'Lane', 'Colour_Chemistry', 'SeqType',
               'Organism', 'Strandedness', 'UDG_Treatment', 'R1', 'R2', 'BAM']]
            .sort_values(['Sample_Name'])
            .to_csv(output[0], sep="\t", index=False)
         )

rule run_eager:
    input:
        tsv = "04-analysis/eager/moderndentalcalc.tsv",
        fa = "03-data/refgenomes/hs37ds.fa"
    output:
        touch("04-analysis/eager/moderndentalcalc.done")
    message: "Run nf-core/EAGER to process the sequencing data of all modern dental calculus libraries"
    params:
        outdir = "04-analysis/eager/moderndentalcalc"
    shell:
        """
        nextflow run nf-core/eager -r 2.4.5 \
            -profile eva,archgen \
            --input "{input.tsv}" \
            --fasta "${{PWD}}/{input.fa}" \
            --skip_deduplication --skip_damage_calculation \
            --complexity_filter_poly_g \
            --skip_collapse \
            --outdir "{params.outdir}"
        """

################################################################################
