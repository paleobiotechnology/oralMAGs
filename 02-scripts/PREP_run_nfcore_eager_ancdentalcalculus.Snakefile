################################################################################
# Project: Oral MAGs
# Part: Data preparation
# Step: Run nf-core/eager to process the ancient dental calculus sequencing data
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
samples = pd.read_csv("01-resources/ancient_dental_calculus_samples.tsv", sep="\t")
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
        expand("04-analysis/eager/ancdentalcalc_{strandedness}.done", strandedness=['dsdna', 'ssdna'])

rule generate_eager_tsv:
    output:
        dsdna = "04-analysis/eager/ancdentalcalc_dsdna.tsv",
        ssdna = "04-analysis/eager/ancdentalcalc_ssdna.tsv",
    message: "Generate the EAGER input TSVs for the ancient dental calculus samples"
    resources:
        mem_mb = 4000
    threads: 1
    run:
        eager = (samples[['archive_data_accession', 'archive_sample_accession',
                          'strand_type', 'library_treatment', 'library_layout',
                          'library_strategy', 'sample_host', 'instrument_model']]
            .rename({'archive_data_accession': 'Library_ID',
                     'archive_sample_accession': 'Sample_Name',
                     'strand_type': 'Strandedness',
                     'library_treatment': 'UDG_Treatment',
                     'library_layout': 'SeqType',
                     'sample_host': 'Organism'},
                    axis=1)
        )

        eager['Lane'] = 1
        eager['Colour_Chemistry'] = [instrument_model[im]
                                     for im in eager['instrument_model'].values]
        eager['SeqType'] = ["PE" if st == "PAIRED" else "SE"
                            for st in eager['SeqType'].values]
        eager['Strandedness'] = ["dsDNA" if st == "double" else "ssDNA"
                                 for st in eager['Strandedness'].values]
        eager['R1'] = "03-data/raw_data/" + eager['Library_ID'] + "_1.fastq.gz"
        eager['R2'] = "03-data/raw_data/" + eager['Library_ID'] + "_2.fastq.gz"
        eager.loc[eager['SeqType'] == "SE", 'R2'] = "NA"
        eager['BAM'] = "NA"

        # dsDNA libraries
        (eager[['Sample_Name', 'Library_ID', 'Lane', 'Colour_Chemistry', 'SeqType',
               'Organism', 'Strandedness', 'UDG_Treatment', 'R1', 'R2', 'BAM']]
            .query('Strandedness == "dsDNA"')
            .sort_values(['Sample_Name'])
            .to_csv(output.dsdna, sep="\t", index=False)
         )
        (eager[['Sample_Name', 'Library_ID', 'Lane', 'Colour_Chemistry', 'SeqType',
               'Organism', 'Strandedness', 'UDG_Treatment', 'R1', 'R2', 'BAM']]
            .query('Strandedness == "ssDNA"')
            .sort_values(['Sample_Name'])
            .to_csv(output.ssdna, sep="\t", index=False)
         )

rule run_eager:
    input:
        tsv = "04-analysis/eager/ancdentalcalc_{strandedness}.tsv",
        fa = "03-data/refgenomes/hs37ds.fa"
    output:
        touch("04-analysis/eager/ancdentalcalc_{strandedness}.done")
    message: "Run nf-core/EAGER to process the sequencing data of all ancient dental calculus {wildcards.strandedness} libraries"
    params:
        revadapter_seq = lambda wildcards: "--clip_reverse_adaptor GGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" if wildcards.strandedness == "ssdna" else ""
        outdir = "04-analysis/eager/ancdentalcalc_{strandedness}"
    shell:
        """
        nextflow run nf-core/eager -r 2.4.5 \
            -profile eva,archgen \
            --input "{input.tsv}" \
            --fasta "${{PWD}}/{input.fa}" \
            --skip_deduplication --skip_damage_calculation \
            --complexity_filter_poly_g \
            --skip_collapse \
            {params.revadapter_seq} \
            --outdir "{params.outdir}"
        """

################################################################################
