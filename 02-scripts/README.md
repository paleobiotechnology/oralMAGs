## Overview of the folder `02-scripts`

This folder contains all the workflows and scripts that are necessary to conduct the experiments of
this project.

### `PREP`: data preparation

- `PREP_prepare_AncientMetagenomeDir_samplelist.Snakefile`: download the list of available dental calculus samples from the AncientMetagenomeDir database
- `PREP_prepare_ancient_dental_calculus_samplelist.Snakefile`: generate the sample list of all (un)published dental calculus samples used in this study
- `PREP_prepare_modern_dental_calculus_samplelist.Snakefile`: generate the sample list of modern dental calculus samples whose MAGs were previously published in Klapper, Hübner, Ibrahim _et al._ (2023)
- `PREP_prepare_saliva_samplelist.Snakefile`: generate the sample list of modern saliva samples whose MAGs were not previously published
- `PREP_download_from_ENA.Snakefile`: download the sequencing data from ENA
- `PREP_run_nfcore_eager_ancdentalcalculus.Snakefile`: pre-process the sequencing data of the ancient dental calculus samples using nf-core/eager
- `PREP_run_nfcore_eager_moderndentalcalculus.Snakefile`: pre-process the sequencing data of the modern dental calculus samples using nf-core/eager
- `PREP_run_nfcore_eager_saliva.Snakefile`: pre-process the sequencing data of the saliva samples using nf-core/eager

The files starting with the prefix `ENVS` are conda environments used in combination with Snakemake.

### `Archaea`: scripts for *Methanobrevibacter*_A and *Methanomethylophilus* analyses and figures

- `aDNA-BAMPlotter_2kb_intervals.py`: scipt modified from https://github.com/MeriamGuellil/aDNA-BAMPlotter to visualize genome mapping across across 2kb base pair intervals (extended data figure 9)
- `extended_data_8_archaea_mags_summary_figure.Rmd`: visualize methanogenic archaea MAG assembly metrics (extended data figure 8)
- `methanobrevibacter_tree_and_metadata.Rmd`: visulize *Methanobrevibacter*_A full genus tree and metadata (supplemental figure 35)
- `methanobrevibacter_cluster1187_tree_and_metadata.Rmd`: visualize Ca. *Methanobrevibacter*_A *prisca* (pr_cl_1187) single copy core gene tree and metadata (supplemental figure 39)
- `methanobrevibacter_cluster1188_tree_and_metadata.Rmd`: visualize *Methanobrevibacter*_A *oralis* (pr_cl_1188) single copy core gene tree and metadata (supplemental figure 37)
- `methanobrevibacter_cluster1189_tree_and_metadata.Rmd`: visualize Ca. *Methanobrevibacter*_A *senecta* (pr_cl_1189) single copy core gene tree and metadata (supplemental figure 38)
- `methanobrevibacter_cluster1190_tree_and_metadata.Rmd`: visualize Ca. *Methanobrevibacter*_A *cohabitans* (pr_cl_1190) single copy core gene tree and metadata (supplemental figure 40)
- `methanobrevibacter_methanogenesis_genes_tree_heatmap.Rmd`: visualize the presence/absence of genes involved in methanogenesis in *Methanobrevibacter*_A MAGs (supplemental figure 46)
- `methanobacteriaceae_species_refs_and_cluster_reps.Rmd`: visualize how *Methanobrevibacter*_A cluster representatives are related to species within the *Methanobacteriaceae* family (supplemental figure 41)
- `methanomethylophilus_tree_and_metadata.Rmd`: visulize *Methanomethylophilus* full genus tree and metadata (supplemental figure 33)
- `methanomethylophilus_cluster3192_tree.Rmd`: visualize Ca. *Methanomethylophilus fidelis* (pr_cl_3192) single copy core gene tree and metadata (supplemental figure 34)
- `archaea_cas_alignments.Rmd`: create MUSCLE and MAFFT percent identity heat maps for *Methanobrevibacter*_A cas IIIA and *Methanomethylophilus* cas IIA (supplemental figures 42 and 43)
- `mad.R`: Minimal Ancestor Deviation (MAD) rooting R script file used for single copy core gene tree rooting
