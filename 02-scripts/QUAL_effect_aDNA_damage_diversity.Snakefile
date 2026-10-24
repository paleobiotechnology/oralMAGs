################################################################################
# Project: Oral MAGs
# Part: Quality evaluation
# Step: Quantify the effect of aDNA damage on genome diversity estimates using
#       simulated sequencing data
#
# Dependent on:
#
# Alex Hübner
################################################################################

import os
from math import ceil

import pandas as pd

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SIMULATION PARAMETERS #####################################################
REFGENOMES = {
    'Tforsythia': 'NC_016610.1',
    'Ssanguinis': 'NZ_LS483385.1',
    'Doralis': 'NZ_CP021255.1',
}
READLENGTHS = [35, 50, 75, 100]
SEQDEPTHS = [5, 10, 20, 50]
DAMAGE = [0, 5, 10, 20, 50]
BRIGGS_PARAMETERS = {
    "0": [0, 0],
    "5": [0.3, 0.07],
    "10": [0.35, 0.15],
    "20": [0.35, 0.32],
    "50": [0.20, 0.64],
}
################################################################################

#### Auxilliary functions ######################################################

def get_nreads(wildcards):
    """ Return the number of reads necessary to obtain a 1000-fold coverage at a
        given read length.
    """
    refspecies = pd.read_csv(checkpoints.summarise_genome_length.get(**wildcards).output[0],
                             sep="\t", index_col=['species'])
    genome_length = refspecies.at[wildcards.species, 'length']
    return ceil(genome_length * 1000 * 1.01 / int(wildcards.rlength))


def calculate_fraction(wildcards):
    """ Return the fraction of reads that is supposed to be sub-sampled
        starting from a 1000-fold coverage
    """
    r = int(wildcards.depth) / 1000
    return f"{r:.3f}"


def calculate_subsampling_seed(wildcards):
    """ Calculate the seed for each replicate that is used for the sub-sampling
    """
    return SEQDEPTHS.index(int(wildcards.depth)) * 10 + int(wildcards.replicate)


def pick_contigs_fasta(wildcards):
    """ Pick the correct input contigs, either the MEGAHIT reported ones or the
        corrected ones, for the minimap2 alignment.
    """
    if wildcards.correction == "raw":
        fasta_path = f"tmp/adna_damage_effect_diversity/assembly/{wildcards.species}-{wildcards.rlength}-{wildcards.damage}-{wildcards.depth}-{wildcards.replicate}/final.contigs.fa"
    else:
        fasta_path = f"tmp/adna_damage_effect_diversity/consensus_correction/{wildcards.species}-{wildcards.rlength}-{wildcards.damage}-{wildcards.depth}-{wildcards.replicate}_contigs.fa"
    return fasta_path

################################################################################

wildcard_constraints:
    species = "[A-Za-z0-9]+",
    rlength = "[0-9]+",
    damage = "[0-9]+",
    depth = "[0-9]+",
    replicate = "[0-9]+",
    correction = "(raw|corrected)"

rule all:
    input:
        "03-data/refgenomes/simulations/genome_lengths_gc.tsv",
        "03-data/adna_damage_effect_diversity/damage_estimate/damageprofiler_misincorporation.tsv",
        "05-results/QUAL_effect_aDNA_damage_diversity.tsv"

#### Prepare genomes for short-read data simulation ############################

rule download_fasta:
    output:
        "03-data/refgenomes/simulations/{species}.fna"
    message: "Download the FastA sequence from NCBI GenBank: {wildcards.species}"
    conda: "ENVS_entrez.yaml"
    resources:
        mem_mb = 4000
    params:
        accid = lambda wildcards: REFGENOMES[wildcards.species]
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

rule genome_length:
    input:
        "03-data/refgenomes/simulations/{species}.fna"
    output:
        temp("03-data/refgenomes/simulations/{species}.length.txt")
    message: "Determine the length of the reference genome: {wildcards.species}"
    container: "https://depot.galaxyproject.org/singularity/bioawk%3A1.0--he4a0461_9"
    threads: 1
    shell:
        """
        bioawk -c fastx -v species={wildcards.species} '{{
            print species "\t" length($seq) "\t" gc($seq)
        }}' {input} > {output}
        """

checkpoint summarise_genome_length:
    input:
        expand("03-data/refgenomes/simulations/{species}.length.txt", species=REFGENOMES)
    output:
        "03-data/refgenomes/simulations/genome_lengths_gc.tsv"
    message: "Summarise the genome lengths of the reference genomes"
    threads: 1
    run:
        (pd.concat([pd.read_csv(fn, sep="\t", header=None, names=['species', 'length', 'gc'])
                    for fn in input])
         .sort_values(['species'])
         .to_csv(output[0], sep="\t", index=False)
        )

################################################################################

#### Simulate genomes using gargammel ##########################################

rule fragment_genomes:
    input:
        length_gc = "03-data/refgenomes/simulations/genome_lengths_gc.tsv",
        fai = "03-data/refgenomes/simulations/{species}.fna.fai"
    output:
        "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}.bam"
    message: "Using Gargammel's FragSim to generate fragments: {wildcards.species} with read length {wildcards.rlength} bp"
    container: "https://depot.galaxyproject.org/singularity/gargammel%3A1.1.2--he905c8f_6"
    resources:
        mem_mb = 8000,
        runtime = 10,
        slurm_partition = "short,standard"
    params:
        nreads = lambda wildcards: get_nreads(wildcards),
        genome = "03-data/refgenomes/simulations/{species}.fna"
    shell:
        """
        fragSim \
            -n {params.nreads} \
            -l {wildcards.rlength} \
            -b {output} \
            {params.genome}
        """

rule deaminated_fragments:
    input:
        "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}.bam"
    output:
        "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam.bam"
    message: "Add aDNA damage: {wildcards.damage}% damage damage for {wildcards.species} with read length {wildcards.rlength} bp"
    container: "https://depot.galaxyproject.org/singularity/gargammel%3A1.1.2--he905c8f_6"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    params:
        damage = lambda wildcards: f"0.03,{BRIGGS_PARAMETERS[wildcards.damage][0]},0.01,{BRIGGS_PARAMETERS[wildcards.damage][1]}"
    shell:
        """
        if [[ "{wildcards.damage}" = "0" ]]; then
            cat {input} > {output}
        else
            deamSim \
            -damage {params.damage} \
            -b {output} \
            {input}
        fi
        """

rule bam_to_fasta:
    input:
        "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam.bam"
    output:
        temp("tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam.fasta")
    message: "Convert from BAM to FastA format: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 2000,
        runtime = 10,
        slurm_partition = "short,standard"
    shell:
        """
        samtools fasta {input} > {output}
        """

rule add_adapters:
    input:
        "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam.fasta"
    output:
        temp("tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam_adapt.fa")
    message: "Run adaptSim to add adapters: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage"
    container: "https://depot.galaxyproject.org/singularity/gargammel%3A1.1.2--he905c8f_6"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    params:
        forward_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGACTATCTCGTATGCCGTCTTCTGCTTG",
        reverse_adapter = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        chemistry = 150
    shell:
        """
        adptSim \
        -f {params.forward_adapter} \
        -s {params.reverse_adapter} \
        -l {params.chemistry} \
        -artp {output} \
        {input}
        """

rule add_illumina_error_to_adapters:
    input:
        "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam_adapt.fa"
    output:
        pe1 = temp("tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam_s_1.fq.gz"),
        pe2 = temp("tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam_s_2.fq.gz")
    message: "Add Illumina sequencing error: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage"
    container: "https://depot.galaxyproject.org/singularity/gargammel%3A1.1.2--he905c8f_6"
    resources:
        mem_mb = 16000,
        runtime = 60,
        slurm_partition = "short,standard"
    params:
        sequencer = "HS25",
        chemistry = 150,
        outputprefix = "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam_s_"
    threads: 4
    shell:
        """
        art_illumina \
            -ss {params.sequencer} \
            -amp -na \
            -p \
            -i {input} \
            -l {params.chemistry} -c 1 \
            -qs 0 \
            -qs2 0 \
            -rs 1337 \
            -o {params.outputprefix}
        gzip {params.outputprefix}*.fq
        """

rule adapterremoval:
    input:
        pe1 = "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam_s_1.fq.gz",
        pe2 = "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_deam_s_2.fq.gz"
    output:
        pe1 = temp("tmp/adna_damage_effect_diversity/reservoir/{species}-{rlength}-{damage}_1.trimmed.fq"),
        pe2 = temp("tmp/adna_damage_effect_diversity/reservoir/{species}-{rlength}-{damage}_2.trimmed.fq")
    message: "Remove adapters: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage"
    container: "https://depot.galaxyproject.org/singularity/adapterremoval%3A2.3.4--pl5321haf24da9_2"
    resources:
        mem_mb = 16000,
        runtime = 120,
        slurm_partition = "short,standard"
    params:
        prefix = "tmp/adna_damage_effect_diversity/reservoir/{species}-{rlength}-{damage}",
        forward_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGACTATCTCGTATGCCGTCTTCTGCTTG",
        reverse_adapter = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    log: "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_adapterremoval.log"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        AdapterRemoval --basename {params.prefix} \
                           --file1 {input.pe1} \
                           --file2 {input.pe2} \
                           --trimns \
                           --trimqualities \
                           --minquality 20 \
                           --minlength 30 \
                           --output1 {output.pe1} \
                           --output2 {output.pe2} \
                           --threads {threads} \
                           --qualitybase 33 \
                           --adapter1 {params.forward_adapter} \
                           --adapter2 {params.reverse_adapter} \
                           --settings {log}
        rm {params.prefix}.discarded {params.prefix}.singleton.truncated
        """

rule fastp:
    input:
        pe1 = "tmp/adna_damage_effect_diversity/reservoir/{species}-{rlength}-{damage}_1.trimmed.fq",
        pe2 = "tmp/adna_damage_effect_diversity/reservoir/{species}-{rlength}-{damage}_2.trimmed.fq"
    output:
        pe1 = "04-analysis/adna_damage_effect_diversity/reservoir/{species}-{rlength}-{damage}_1.fq.gz",
        pe2 = "04-analysis/adna_damage_effect_diversity/reservoir/{species}-{rlength}-{damage}_2.fq.gz"
    message: "Clean reads with FastP: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage"
    container: "https://depot.galaxyproject.org/singularity/fastp%3A1.0.1--heae3180_0"
    resources:
        mem_mb = 8000,
        runtime = 30,
        slurm_partition = "short,standard"
    log: "tmp/adna_damage_effect_diversity/gargammel/{species}-{rlength}-{damage}_fastp.json"
    threads: 8
    shell:
        """
        fastp --in1 {input.pe1} \
              --in2 {input.pe2} \
              --out1 {output.pe1} \
              --out2 {output.pe2} \
              -A -g -Q -L \
              -w {threads} --json {log} --html /dev/null
        """

################################################################################

#### Sub-sample to different sequencing depth 10 replicates each ###############

rule subsample:
    input:
        pe1 = "04-analysis/adna_damage_effect_diversity/reservoir/{species}-{rlength}-{damage}_1.fq.gz",
        pe2 = "04-analysis/adna_damage_effect_diversity/reservoir/{species}-{rlength}-{damage}_2.fq.gz"
    output:
        pe1 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_1.fq.gz",
        pe2 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_2.fq.gz"
    message: "Subsample FastQ files: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/seqtk:1.4--he4a0461_1"
    resources:
        mem_mb = 120000,
        runtime = 120,
        slurm_partition = "short,standard",
    params:
        frac = lambda wildcards: calculate_fraction(wildcards),
        seed = lambda wildcards: calculate_subsampling_seed(wildcards)
    threads: 1
    shell:
        """
        gunzip -c {input.pe1} | seqtk sample -s{params.seed} - {params.frac} | gzip > {output.pe1}
        gunzip -c {input.pe2} | seqtk sample -s{params.seed} - {params.frac} | gzip > {output.pe2}
        """

################################################################################

#### Damage analysis ###########################################################

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
        pe1 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_1.fq.gz",
        pe2 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_2.fq.gz"
    output:
        temp("tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}_{i}.sai")
    message: "Align the sequencing data against the reference genome using BWA: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/bwa:0.7.17--he4a0461_11"
    resources:
        mem_mb = 4000,
        runtime = 60,
        slurm_partition = "short,standard"
    params:
        idx = "03-data/refgenomes/simulations/{species}",
        fastq = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_{i}.fq.gz"
    threads: 1
    shell:
        """
        bwa aln -t {threads} \
            {params.idx} \
            {params.fastq} \
            -n 0.01 -k 2 -l 1024 -o 2 \
            -f {output} 
        """

rule bwa_sampe:
    input:
        bwa = "03-data/refgenomes/simulations/{species}.amb",
        pe1 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_1.fq.gz",
        pe2 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_2.fq.gz",
        sai1 = "tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}_1.sai",
        sai2 = "tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}_2.sai",
    output:
        temp("tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}.sam")
    message: "Generate BAM file: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
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
        bwa sampe \
            {params.idx} \
            {input.sai1} {input.sai2} \
            {input.pe1} {input.pe2} > {output}
        """

rule bam2sam:
    input:
        "tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}.sam"
    output:
        temp("tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}.bam")
    message: "Convert SAM to BAM: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
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
        "tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}.bam"
    output:
        bam = "tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}.sorted.bam",
        bai = "tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}.sorted.bam.bai"
    message: "Sort BAM file by coordinate: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
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

rule damageprofiler:
    input:
        fa = "03-data/refgenomes/simulations/{species}.fna",
        bam = "tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}.sorted.bam",
        bai = "tmp/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}.sorted.bam.bai"
    output:
        "04-analysis/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}/editDistance.txt"
    message: "Perform ancient DNA damage and read length analysis: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "https://depot.galaxyproject.org/singularity/damageprofiler%3A1.1--hdfd78af_2"
    resources:
        mem_mb = 64000,
        runtime = 60,
        slurm_partition = "short,standard"
    params:
        dir = "04-analysis/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}"
    threads: 4
    shell:
        """
        damageprofiler -Xmx48g -i {input.bam} \
            -o {params.dir} \
            -r {input.fa}
        """

rule summarise_mismatch_incorporation:
    input:
        expand("04-analysis/adna_damage_effect_diversity/damage_estimate/{species}-{rlength}-{damage}-{depth}-{replicate}/editDistance.txt",
               species=REFGENOMES,
               rlength=READLENGTHS,
               damage=DAMAGE,
               depth=[10],
               replicate=[0])
    output:
        "04-analysis/adna_damage_effect_diversity/damage_estimate/damageprofiler_misincorporation.tsv"
    message: "Summarise the misincorporation patterns calculated from Damageprofiler"
    resources:
        mem_mb = 4000,
        runtime = 60,
        slurm_partition = "short,standard"
    threads: 1
    run:
        misincorps = []
        for fn in input:
            species, rlength, damage, _, _ = os.path.basename(os.path.dirname(fn)).split("-")
            misincorp_5p = (pd.read_csv(f"{os.path.dirname(fn)}/5p_freq_misincorporations.txt",
                                        sep="\t", skiprows=3)
                .assign(species=species)
                .assign(rlength=int(rlength))
                .assign(damage=int(damage))
                .assign(end="5p")
            )
            misincorps.append(misincorp_5p.iloc[:, [15, 16, 17, 18] + list(range(15))])
            misincorp_3p = (pd.read_csv(f"{os.path.dirname(fn)}/3p_freq_misincorporations.txt",
                                        sep="\t", skiprows=3)
                .assign(species=species)
                .assign(rlength=int(rlength))
                .assign(damage=int(damage))
                .assign(end="3p")
            )
            misincorps.append(misincorp_3p.iloc[:, [15, 16, 17, 18] + list(range(15))])

        misincorps_df = (
            pd.concat(misincorps)
            .rename({'Pos': 'pos'}, axis=1)
            .sort_values(['species', 'rlength', 'damage', 'end', 'pos'])
            )
        misincorps_df.to_csv(output[0], sep="\t", index=False)

################################################################################

#### De novo assembly ##########################################################

rule megahit:
    input:
        pe1 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_1.fq.gz",
        pe2 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_2.fq.gz"
    output:
        "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa"
    message: "Assemble the metagenomic data using MEGAHIT: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/megahit:1.2.9--haf24da9_8"
    resources:
        mem_mb = 100000,
        runtime = 2880,
        slurm_partition = "standard"
    params: 
        tmpdir = "tmp/adna_damage_effect_diversity/tmp/{species}-{rlength}-{damage}-{depth}-{replicate}",
        prefix = "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}",
        minlength = 500
    threads: 4
    shell:
        """
        if [[ -d {params.prefix} ]]; then
            rm -r {params.prefix}
        fi
        if [[ ! -d {params.tmpdir} ]]; then
            mkdir -p {params.tmpdir}
        fi
        megahit \
            -t {threads} \
            -1 {input.pe1} \
            -2 {input.pe2} \
            --tmp-dir {params.tmpdir} \
            --min-contig-len {params.minlength} \
            --out-dir {params.prefix}
        """

################################################################################

#### Consensus correction ######################################################

rule build_BowTie2_index:
    input:
        "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa"
    output:
        bt1 = temp("tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.1.bt2"),
        bt2 = temp("tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.2.bt2"),
        bt3 = temp("tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.3.bt2"),
        bt4 = temp("tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.4.bt2"),
        bt_rev1 = temp("tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.rev.1.bt2"),
        bt_rev2 = temp("tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.rev.2.bt2")
    message: "Index the contigs for alignment using BowTie2: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/bowtie2:2.5.5--ha27dd3b_0"
    resources:
        mem_mb = 8000,
        runtime = 1440,
        slurm_partition = "standard",
        cores = 4
    params:
        index = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}"
    threads: 4
    shell:
        """
        if [[ -s {input} ]]; then
            bowtie2-build -f \
                {input} \
                {params.index}
        else
            touch {output.bt1} {output.bt2} {output.bt3} {output.bt4} {output.bt_rev1} {output.bt_rev2}
        fi
        """

rule BowTie2_alignment:
    input:
        fa = "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa",
        bt1 = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.1.bt2",
        bt2 = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.2.bt2",
        bt3 = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.3.bt2",
        bt4 = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.4.bt2",
        bt_rev1 = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.rev.1.bt2",
        bt_rev2 = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.rev.2.bt2",
        pe1 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_1.fq.gz",
        pe2 = "tmp/adna_damage_effect_diversity/subsamples/{species}-{rlength}-{damage}-{depth}-{replicate}_2.fq.gz"
    output:
        temp("tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.raw.sam")
    message: "Align reads back to the uncorrected contigs using BowTie2's very-sensitive setting: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/bowtie2:2.5.5--ha27dd3b_0"
    group: "ref_alignment"
    resources:
        mem_mb = 16000,
        runtime = 2880,
        slurm_partition = "standard",
        cores = 16
    params:
        index = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}"
    threads: 16
    shell:
        """
        if [[ -s {input.fa} ]]; then
            bowtie2 -p {threads} --very-sensitive -N 1 -x {params.index} \
                    -1 {input.pe1} -2 {input.pe2} -S {output}
        else
            touch {output}
        fi
        """

rule samtools_sort:
    input:
        bam = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.raw.sam",
        fasta = "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa"
    output:
        bam = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.sorted.noncorr.bam",
        bai = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.sorted.noncorr.bam.bai"
    message: "Sort the sequencing data: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    group: "ref_alignment"
    resources:
        mem_mb = 8000,
        runtime = 1440,
        slurm_partition = "standard",
        cores = 2
    threads: 2
    shell:
        """
        if [[ -s {input.fasta} ]]; then
            samtools view -Sb {input.bam} | \
            samtools calmd -u /dev/stdin {input.fasta} | \
            samtools sort -l 4 -o {output.bam} -
            samtools index {output.bam}
        else
            touch {output.bam} {output.bai}
        fi
        """

rule faidx_reffasta:
    input:
        "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa"
    output:
        "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa.fai"
    message: "Generate FastA index for de-novo assembled contigs: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 4000,
        runtime = 120,
        slurm_partition = "short"
    shell:
        """
        if [[ -s {input} ]]; then
            samtools faidx {input}
        else
            touch {output}
        fi
        """

rule freebayes_shortreads:
    input:
        fa = "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa",
        fai = "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa.fai",
        bam = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.sorted.noncorr.bam",
        bai = "tmp/adna_damage_effect_diversity/alignment/{species}-{rlength}-{damage}-{depth}-{replicate}.sorted.noncorr.bam.bai"
    output:
        temp("tmp/adna_damage_effect_diversity/consensus_correction/{species}-{rlength}-{damage}-{depth}-{replicate}.vcf")
    message: "Genotype the contigs using freeBayes using the short-read data: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/freebayes:1.3.9--hbefcdb2_3"
    group: "freebayes"
    resources:
        mem_mb = 32000,
        runtime = 5760,
        slurm_partition = "standard"
    threads: 1
    shell:
        """
        if [[ -s {input.fa} ]]; then
            freebayes -f {input.fa} -F 0.33 -p 1 -q 20 {input.bam} > {output}
        else
            touch {output}
        fi
        """

rule bcftools_filter:
    input:
        fa = "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa",
        vcf = "tmp/adna_damage_effect_diversity/consensus_correction/{species}-{rlength}-{damage}-{depth}-{replicate}.vcf"
    output:
        vcf = "tmp/adna_damage_effect_diversity/consensus_correction/{species}-{rlength}-{damage}-{depth}-{replicate}.filter.vcf.gz",
        tbi = "tmp/adna_damage_effect_diversity/consensus_correction/{species}-{rlength}-{damage}-{depth}-{replicate}.filter.vcf.gz.tbi"
    message: "Discard low-quality differences between MEGAHIT and freebayes consensus: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0"
    resources:
        mem_mb = 4000,
        runtime = 120,
        slurm_partition = "short",
        cores = 1
    shell:
        """
        if [[ -s {input.fa} ]]; then
            bcftools view \
                -v snps,mnps \
                -i 'QUAL >= 30 || (QUAL >= 20 && INFO/AO >= 3)' {input.vcf} | \
            bgzip > {output.vcf}
            bcftools index -t {output.vcf}
        else
            touch {output.vcf} {output.tbi}
        fi
        """

rule bcftools_consensus:
    input:
        fa = "tmp/adna_damage_effect_diversity/assembly/{species}-{rlength}-{damage}-{depth}-{replicate}/final.contigs.fa",
        vcf = "tmp/adna_damage_effect_diversity/consensus_correction/{species}-{rlength}-{damage}-{depth}-{replicate}.filter.vcf.gz",
        tbi = "tmp/adna_damage_effect_diversity/consensus_correction/{species}-{rlength}-{damage}-{depth}-{replicate}.filter.vcf.gz.tbi"
    output:
        "tmp/adna_damage_effect_diversity/consensus_correction/{species}-{rlength}-{damage}-{depth}-{replicate}_contigs.fa"
    message: "Correct the consensus sequence of the contigs: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate}"
    container: "docker://quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0"
    resources:
        mem_mb = 8000,
        runtime = 240,
        slurm_partition = "short",
    threads: 2
    shell:
        """
        if [[ -s {input.fa} ]]; then
            cat {input.fa} | bcftools consensus {input.vcf} > {output}
        else
            touch {output}
        fi
        """

################################################################################

#### Align with minimap2 against the reference genome and genotype #############

rule minimap2:
    input:
        fa = lambda wildcards: pick_contigs_fasta(wildcards),
        reffa = "03-data/refgenomes/simulations/{species}.fna",
        reffai = "03-data/refgenomes/simulations/{species}.fna.fai",
    output:
        temp("tmp/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.sam")
    message: "Align the de novo assembled contigs against the reference genome: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate} - {wildcards.correction}"
    container: "docker://quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
    resources:
        mem_mb = 4000,
        runtime = 120,
        slurm_partition = "short,standard"
    params:
        readgroup = lambda wildcards: rf'@RG\\tID:{wildcards.species}-{wildcards.rlength}-{wildcards.damage}-{wildcards.depth}-{wildcards.replicate}-{wildcards.correction}\\tSM:{wildcards.species}-{wildcards.rlength}-{wildcards.damage}-{wildcards.depth}-{wildcards.replicate}-{wildcards.correction}',
    threads: 1
    shell:
        """
        if [[ -s {input.fa} ]]; then
            minimap2 -a -x asm10 -o {output} -R {params.readgroup} {input.reffa} {input.fa}
        else
            touch {output}
        fi
        """

rule sam2bam:
    input:
        fa = lambda wildcards: pick_contigs_fasta(wildcards),
        sam = "tmp/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.sam"
    output:
        bam = "04-analysis/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.bam",
        bai = "04-analysis/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.bam.bai"
    message: "Convert SAM to BAM file: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate} - {wildcards.correction}"
    container: "docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"
    resources:
        mem_mb = 2000,
        runtime = 120,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        """
        if [[ -s {input.fa} ]]; then
            samtools view -bS {input.sam} | samtools sort -o {output.bam} -
            samtools index {output.bam}
        else
            touch {output}
        fi
        """

rule mpileup:
    input:
        fa = lambda wildcards: pick_contigs_fasta(wildcards),
        reffa = "03-data/refgenomes/simulations/{species}.fna",
        reffai = "03-data/refgenomes/simulations/{species}.fna.fai",
        bam = "04-analysis/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.bam",
        bai = "04-analysis/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.bam.bai"
    output:
        "tmp/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.mpileup.vcf.gz"
    message: "Generate pileup: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate} - {wildcards.correction}"
    container: "docker://quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0"
    resources:
        mem_mb = 8000,
        runtime = 120,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        """
        if [[ -s {input.fa} ]]; then
            bcftools mpileup -f {input.reffa} -O z -o {output} {input.bam}
        else
            touch {output}
        fi
        """

rule call:
    input:
        fa = lambda wildcards: pick_contigs_fasta(wildcards),
        vcf = "tmp/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.mpileup.vcf.gz"
    output:
        vcf = "04-analysis/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.vcf.gz",
        tbi = "04-analysis/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.vcf.gz.tbi"
    message: "Call genotypes: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate} - {wildcards.correction}"
    container: "docker://quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0"
    resources:
        mem_mb = 8000,
        runtime = 120,
        slurm_partition = "short,standard"
    threads: 1
    shell:
        """
        if [[ -s {input.fa} ]]; then
            bcftools call --ploidy 1 -c -O z -o {output.vcf} {input.vcf}
            bcftools index -t {output.vcf}
        else
            touch {output.vcf} {output.tbi}
        fi
        """

rule evaluate_mismatches:
    input:
        fa = lambda wildcards: pick_contigs_fasta(wildcards),
        vcf = "04-analysis/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.vcf.gz",
        tbi = "04-analysis/adna_damage_effect_diversity/minimap/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.vcf.gz.tbi"
    output:
        "04-analysis/adna_damage_effect_diversity/mismatches/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.txt"
    message: "Summarise the number of mismatches: {wildcards.species} with read length {wildcards.rlength} bp and {wildcards.damage}% damage to depth {wildcards.depth}-fold - replicate {wildcards.replicate} - {wildcards.correction}"
    container: "docker://quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0"
    resources:
        mem_mb = 8000,
        runtime = 120,
        slurm_partition = "short,standard"
    params:
        sample = lambda wildcards: f"{wildcards.species}-{wildcards.rlength}-{wildcards.damage}-{wildcards.depth}-{wildcards.replicate}-{wildcards.correction}"
    threads: 1
    shell:
        """
        if [[ -s {input.fa} ]]; then
            total_sites=$(bcftools view -H {input.vcf} | wc -l)
            mismatches=$(bcftools view -H -m 2 {input.vcf} | wc -l)
        else
            total_sites=0
            mismatches=0
        fi
        echo -e "{params.sample}\t${{total_sites}}\t${{mismatches}}" > {output}
        """

rule concat_mismatch_summary:
    input:
        expand("04-analysis/adna_damage_effect_diversity/mismatches/{species}-{rlength}-{damage}-{depth}-{replicate}-{correction}.txt",
               species=REFGENOMES,
               rlength=READLENGTHS,
               damage=DAMAGE,
               depth=SEQDEPTHS,
               replicate=list(range(10)),
               correction=["raw", "corrected"])
    output:
        "05-results/QUAL_effect_aDNA_damage_diversity.tsv"
    message: "Summarise the results of inferring the number of mismatches for each condition of the simulation experiment"
    run:
        res = []
        for fn in input:
            with open(fn, "rt") as f:
                line = f.readline().rstrip()
                sample_id, total, mismatches = line.split("\t")
                sample_info = sample_id.split("-")
                sample_info.extend([total, mismatches])
                res.append(tuple(sample_info))

        df = pd.DataFrame(res, columns=['species', 'readlength', 'aDNA_damage',
                                        'sequencing_depth', 'replicate', 'condition',
                                        'sites_reconstructed', 'total_mismatches'])
        for c in ['readlength', 'aDNA_damage', 'sequencing_depth', 'replicate',
                  'sites_reconstructed', 'total_mismatches']:
            df[c] = df[c].astype(int)
        df['mismatches_per_100kb'] = df['total_mismatches'] / (df['sites_reconstructed'] / 100000)
        (df.sort_values(['species', 'readlength', 'aDNA_damage', 'sequencing_depth', 'replicate'])
         .to_csv(output[0], sep="\t", index=False, float_format="%.2f")
        )

################################################################################
