## Overview of the folder `02-scripts`

This folder contains all the workflows and scripts that are necessary to conduct the experiments of
this project.

### `PREP`: data preparation

- `PREP_prepare_AncientMetagenomeDir_samplelist.Snakefile`: download the list of available dental calculus samples from the AncientMetagenomeDir database
- `PREP_prepare_ancient_dental_calculus_samplelist.Snakefile`: generate the sample list of all (un)published dental calculus samples used in this study
- `PREP_prepare_saliva_samplelist.Snakefile`: generate the sample list of modern saliva samples whose MAGs were not previously published

The files starting with the prefix `ENVS` are conda environments used in combination with Snakemake.
