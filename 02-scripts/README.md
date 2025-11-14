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
