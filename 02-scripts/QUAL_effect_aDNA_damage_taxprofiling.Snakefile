################################################################################
# Project: Oral MAGs
# Part: Quality evaluation
# Step: Quantify the effect of aDNA damage on the taxonomic profiling with
#       Kraken2 using simulated sequencing data
#
# Dependent on:
#   - PREP_run_nfcore_eager_moderndentalcalculus.Snakefile
#
# Alex Hübner
################################################################################

import os
from math import ceil

import pandas as pd

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SIMULATION PARAMETERS #####################################################
# modern dental calculus
samples = (pd.read_csv("01-resources/modern_dental_calculus_samples.tsv", sep="\t")
    .query('sample_alias.str.startswith("JAE")')
    .set_index('secondary_sample_accession')
)
# comparative ancient dental calculus
profiles = {
    'young': ['ERR6178958'], # CDM050.L0101
    'old': ['ERR3579746'], # FUM002.A0101
}
# microbial species for damage estimation
refgenomes = {
    'Doralis': 'NZ_CP021255.1',
}
################################################################################

wildcard_constraints:
    species = "[A-Za-z0-9]+",
    profile = "(old|young)"

rule all:
    input:
        expand("04-analysis/adna_damage_effect_taxprofiling/jae_samples/{sample}-{profile}-{rep}.frag.deam.fastq.gz", sample=samples.index.tolist(), profile=profiles, rep=list(range(10))),
        expand("04-analysis/adna_damage_effect_taxprofiling/jae_samples/{sample}.original.fastq.gz", sample=samples.index.tolist()),
        "05-results/QUAL_effect_aDNA_damage_taxprofiling_fragmentlength.tsv"

#### Infer DNA molecule length distribution and damage profiles ################

rule samtools_fastq:
    output:
        pe1 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/ref_data/{profile}_1.fastq.gz",
        pe2 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/ref_data/{profile}_2.fastq.gz"
    message: "Extract the read pairs that were not aligned against the human reference genome in FastQ format: {wildcards.profile}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    params:
        bam = lambda wildcards: " ".join([f"04-analysis/eager/ancdentalcalc_dsdna/mapping/bwa/{err}_PE.mapped.bam" for err in profiles[wildcards.profile]])
    threads: 1
    shell:
        """
        samtools cat {params.bam} | \
        samtools sort -n - | \
        samtools view -uh -e '(flag.paired && (flag.unmap || flag.munmap)) || (!flag.paired && flag.unmap)' - | \
        samtools fastq -1 {output.pe1} \
                       -2 {output.pe2} -
        """

rule fastp_dna_molecule_length:
    input:
        pe1 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/ref_data/{profile}_1.fastq.gz",
        pe2 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/ref_data/{profile}_2.fastq.gz"
    output:
        merged = "/mnt/local_scratch/adna_damage_effect_taxprofiling/ref_data/{profile}.merged.fastq.gz",
        not_merged = "/mnt/local_scratch/adna_damage_effect_taxprofiling/ref_data/{profile}.not_merged.fastq.gz",
        json = "04-analysis/adna_damage_effect_taxprofiling/ref_data/{profile}.fastp.json"
    message: "Infer the DNA molecule length by merging overlapping read pairs: {wildcards.profile}"
    container: "https://depot.galaxyproject.org/singularity/fastp%3A1.0.1--heae3180_0"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        """
        fastp --in1 {input.pe1} --in2 {input.pe2} \
            --out1 {output.not_merged} --out2 /dev/null --merged_out {output.merged} \
            --merge -A -G -Q -l 30 --overlap_len_require 11 \
            --json {output.json} --html /dev/null
        """

rule count_unmerged_reads:
    input:
        "/mnt/local_scratch/adna_damage_effect_taxprofiling/ref_data/{profile}.not_merged.fastq.gz"
    output:
        "04-analysis/adna_damage_effect_taxprofiling/ref_data/{profile}.not_merged.txt"
    message: "Count the number of reads that were not merged: {wildcards.profile}"
    container: "https://depot.galaxyproject.org/singularity/bioawk%3A1.0--he4a0461_9"
    resources:
        mem_mb = 4000
    threads: 1
    shell:
        """
        bioawk -c fastx 'END{{print NR}}' {input} > {output}
        """

rule download_fasta:
    output:
        "03-data/refgenomes/simulations/{species}.fna"
    message: "Download the FastA sequence from NCBI GenBank: {wildcards.species}"
    conda: "ENVS_entrez.yaml"
    resources:
        mem_mb = 4000
    params:
        accid = lambda wildcards: refgenomes[wildcards.species]
    threads: 1
    shell:
        """
        efetch -db nucleotide -id {params.accid} -mode text -format fasta > {output}
        """

rule fasta_index:
    input:
        "03-data/refgenomes/simulations/{species}.fna"
    output:
        "03-data/refgenomes/simulations/{species}.fna.fai"
    message: "Create FastA index: {wildcards.species}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 4000,
        runtime = 10,
        slurm_partition = "short"
    threads: 1
    shell:
        """
        samtools faidx {input}
        """

rule bwa_index:
    input:
        "03-data/refgenomes/simulations/{species}.fna"
    output:
        "03-data/refgenomes/simulations/{species}.amb"
    message: "Index the reference genome for the alignment with BWA: {wildcards.species}"
    container: "docker://quay.io/biocontainers/bwa:0.7.17--he4a0461_11"
    resources:
        mem_mb = 4000,
        runtime = 60,
        slurm_partition = "short,standard"
    params:
        idx = "03-data/refgenomes/simulations/{species}"
    threads: 1
    shell:
        "bwa index -p {params.idx} {input}"

rule bwa_aln:
    input:
        bwa = "03-data/refgenomes/simulations/{species}.amb",
        fastq = "/mnt/local_scratch/adna_damage_effect_taxprofiling/ref_data/{profile}.merged.fastq.gz"
    output:
        temp("/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.sai")
    message: "Align the sequencing data against the reference genome using BWA: {wildcards.profile} against the species {wildcards.species}"
    container: "docker://quay.io/biocontainers/bwa:0.7.17--he4a0461_11"
    resources:
        mem_mb = 4000,
        runtime = 60,
        slurm_partition = "short,standard"
    params:
        idx = "03-data/refgenomes/simulations/{species}"
    threads: 1
    shell:
        """
        bwa aln -t {threads} \
            {params.idx} \
            {input.fastq} \
            -n 0.01 -k 2 -l 1024 -o 2 \
            -f {output} 
        """

rule bwa_samse:
    input:
        bwa = "03-data/refgenomes/simulations/{species}.amb",
        fastq = "/mnt/local_scratch/adna_damage_effect_taxprofiling/ref_data/{profile}.merged.fastq.gz",
        sai = "/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.sai"
    output:
        temp("/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.sam")
    message: "Generate BAM file: {wildcards.profile} against species {wildcards.species}"
    container: "docker://quay.io/biocontainers/bwa:0.7.17--he4a0461_11"
    resources:
        mem_mb = 4000,
        runtime = 60,
        slurm_partition = "short,standard"
    params:
        idx = "03-data/refgenomes/simulations/{species}"
    threads: 1
    shell:
        """
        bwa samse \
            {params.idx} \
            {input.sai} {input.fastq} > {output}
        """

rule bam2sam:
    input:
        "/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.sam"
    output:
        temp("/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.bam")
    message: "Convert SAM to BAM: {wildcards.profile} against the species {wildcards.species}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 4000,
        runtime = 60,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        "samtools view -Su {input} > {output}"

rule samtools_sort_damage:
    input:
        "/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.bam"
    output:
        bam = "/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.sorted.bam",
        bai = "/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.sorted.bam.bai"
    message: "Sort BAM file by coordinate: {wildcards.profile} against the species {wildcards.species}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 4000,
        runtime = 60,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        """
        samtools sort -o {output.bam} {input} && samtools index {output.bam}
        """

rule mapdamage2:
    input:
        fasta = "03-data/refgenomes/simulations/{species}.fna",
        fai = "03-data/refgenomes/simulations/{species}.fna.fai",
        bam = "/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.sorted.bam",
        bai = "/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}.sorted.bam.bai"
    output:
        "/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}/misincorporation.txt"
    message: "Determine the misincorporation profile using mapDamage2: {wildcards.profile} against the species {wildcards.species}"
    container: "docker://quay.io/biocontainers/mapdamage2:2.2.3--py312h4711d71_0"
    resources:
        mem_mb = 16000,
        runtime = 120,
        slurm_partition = "short,standard"
    params:
        dir = "/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{profile}-{species}"
    threads: 1
    shell:
        """
        mapDamage -i {input.bam} \
            -r {input.fasta} \
            -d {params.dir} \
            --no-stats
        """

################################################################################

#### Simulate genomes using gargammel ##########################################

rule samtools_fastq_jae:
    output:
        pe1 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}_1.fastq.gz",
        pe2 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}_2.fastq.gz"
    message: "Extract the read pairs that were not aligned against the human reference genome in FastQ format: {wildcards.sample}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    params:
        bam = "04-analysis/eager/moderndentalcalc/mapping/bwa/{sample}_PE.mapped.bam"
    threads: 1
    shell:
        """
        samtools cat {params.bam} | \
        samtools sort -n - | \
        samtools view -uh -e '(flag.paired && (flag.unmap || flag.munmap)) || (!flag.paired && flag.unmap)' - | \
        samtools fastq -1 {output.pe1} \
                       -2 {output.pe2} -
        """

rule fastp_merge_overlapping_jae:
    input:
        pe1 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}_1.fastq.gz",
        pe2 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}_2.fastq.gz"
    output:
        merged = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}.merged.fastq.gz",
        pe1 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}.R1.fastq.gz",
        pe2 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}.R2.fastq.gz",
        json = "04-analysis/adna_damage_effect_taxprofiling/jae_samples/{sample}.fastp.json"
    message: "Merging overlapping read pairs: {wildcards.sample}"
    container: "https://depot.galaxyproject.org/singularity/fastp%3A1.0.1--heae3180_0"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        """
        fastp --in1 {input.pe1} --in2 {input.pe2} \
            --out1 {output.pe1} --out2 {output.pe2} --merged_out {output.merged} \
            --merge -A -G -Q -l 30 --overlap_len_require 11 \
            --json {output.json} --html /dev/null
        """

rule concat_fastqs_after_merging:
    input:
        merged = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}.merged.fastq.gz",
        pe1 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}.R1.fastq.gz",
        pe2 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}.R2.fastq.gz",
    output:
        "04-analysis/adna_damage_effect_taxprofiling/jae_samples/{sample}.original.fastq.gz"
    message: "Concatenate the FastQs into a single file: {wildcards.sample}"
    resources:
        mem_mb = 4000,
        runtime = 60,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        "cat {input.merged} {input.pe1} {input.pe2} > {output}"


rule trim_sequences:
    input:
        merged = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}.merged.fastq.gz",
        pe1 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}.R1.fastq.gz",
        pe2 = "/mnt/local_scratch/adna_damage_effect_taxprofiling/jae_samples/{sample}.R2.fastq.gz",
        short = "04-analysis/adna_damage_effect_taxprofiling/ref_data/{profile}.fastp.json",
        long = "04-analysis/adna_damage_effect_taxprofiling/jae_samples/{sample}.fastp.json",
        not_merged = "04-analysis/adna_damage_effect_taxprofiling/ref_data/{profile}.not_merged.txt"
    output:
        "/mnt/local_scratch/adna_damage_effect_taxprofiling/gargammel/{sample}-{profile}-{rep}.frag.fastq.gz",
    message: "Trim the sequences based on the fragment length distribution of an ancient samples: replicate {wildcards.rep} of sample {wildcards.sample} for profile {wildcards.profile}"
    conda: "ENVS_python.yaml"
    resources:
        mem_mb = 16000,
        runtime = 120,
        slurm_partition = "short,standard"
    threads: 1
    script:
        "pyscripts/QUAL_effect_aDNA_damage_taxprofiling-trim_sequences.py"

rule fastq_2_bam:
    input:
        "/mnt/local_scratch/adna_damage_effect_taxprofiling/gargammel/{sample}-{profile}-{rep}.frag.fastq.gz"
    output:
        "/mnt/local_scratch/adna_damage_effect_taxprofiling/gargammel/{sample}-{profile}-{rep}.frag.bam"
    message: "Convert FastQ into BAM: replicate {wildcards.rep} of sample {wildcards.sample} for profile {wildcards.profile}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        """
        samtools import -o {output} {input}
        """

rule deaminated_fragments:
    input:
        bam = "/mnt/local_scratch/adna_damage_effect_taxprofiling/gargammel/{sample}-{profile}-{rep}.frag.bam",
        misinc = lambda wildcards: f"/mnt/local_scratch/adna_damage_effect_taxprofiling/damage_estimate/{wildcards.profile}-Doralis/misincorporation.txt"
    output:
        temp("04-analysis/adna_damage_effect_taxprofiling/jae_samples/{sample}-{profile}-{rep}.frag.deam.bam")
    message: "Add aDNA damage: replicate {wildcards.rep} of sample {wildcards.sample} for profile {wildcards.profile}"
    container: "https://depot.galaxyproject.org/singularity/gargammel%3A1.1.2--he905c8f_6"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    shell:
        """
        deamSim \
            -mapdamage {input.misinc} "double" \
            -b {output} \
            {input.bam}
        """

rule bam_2_fastq:
    input:
        "04-analysis/adna_damage_effect_taxprofiling/jae_samples/{sample}-{profile}-{rep}.frag.deam.bam"
    output:
        "04-analysis/adna_damage_effect_taxprofiling/jae_samples/{sample}-{profile}-{rep}.frag.deam.fastq.gz"
    message: "Convert BAM to FastQ: replicate {wildcards.rep} of sample {wildcards.sample} for profile {wildcards.profile}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        """
        samtools fastq {input} | gzip > {output}
        """

################################################################################

#### Summarise #################################################################

rule fragment_length_dist:
    input:
        modern = expand("04-analysis/adna_damage_effect_taxprofiling/jae_samples/{sample}.fastp.json", sample=samples.index.tolist()),
        ancient = expand("04-analysis/adna_damage_effect_taxprofiling/ref_data/{profile}.fastp.json", profile=profiles)
    output:
        "05-results/QUAL_effect_aDNA_damage_taxprofiling_fragmentlength.tsv"
    message: "Summarise the fragment length distributions"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    threads: 1
    script:
        "pyscripts/QUAL_effect_aDNA_damage_taxprofiling-fragment_length_dist.py"

################################################################################
