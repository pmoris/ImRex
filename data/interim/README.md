# Processed files

See `./data/raw/README.md` for more information on the origin of the raw data.

## vdjdb-2019-08-08

Contains datasets for the processed VDJdb files in `./data/raw/vdjdb/vdjdb-2019-08-08`. This data can be generated by running `make preprocess-vdjdb-aug-2019`, after retrieving the raw data.

All different variants are defined in the `./Makefile`. Every variant will have a normal and `_full` `.csv` file, a `.log` showing the commands with which it was created and a `.pdf` visualisation of its contents.

## immunecode-adaptive

Contains datasets for the processed Adaptive ImmuneCODE files in `./data/raw/immunecode-adaptive`. This data can be generated by running `make preprocess-adaptive`, after retrieving the raw data.

## mcpas

Contains datasets with human CDR3-epitope pairs taken from McPAS ([http://friedmanlab.weizmann.ac.il/McPAS-TCR/](http://friedmanlab.weizmann.ac.il/McPAS-TCR/)), as well as two subsets of it: sequence pairs with epitopes that occur in the VDJdb dataset and those that do not. These files can be generated by running `make preprocess-mcpas`. The raw data is included in this repository, and has underwent a number of quality checks already.