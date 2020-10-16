ImRex (Interaction Map Recognition)
==============================

TCR-epitope recognition prediction using combined sequence input represention for convolutional neural networks
------------------------------

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [Project Organization](#project-organization)
- [Usage](#usage)
  - [Training](#training)
  - [Predictions using the pre-built model](#predictions-using-the-pre-built-model)
  - [Results](#results)
- [Zenodo data repository](#zenodo-data-repository)
- [Citation](#citation)
- [Authors](#authors)
- [History of the project](#history-of-the-project)

<!-- /code_chunk_output -->



## Project Organization

    ├── LICENSE                     <- MIT license.
    ├── Makefile                    <- Makefile with commands like `make data`.
    ├── README.md                   <- The top-level README for developers using this project.
    ├── data
    │   ├── interim                 <- Intermediate data that has been transformed.
    │   └── raw                     <- The original, immutable data dump. See README for more info.
    │
    ├── models                      <- Directory wheret rained models, model predictions, summaries and
    │   │                              evaluation figures will be stored.
    │   └── models-pretrained       <- Contains a small number of pre-trained models.
    │
    ├── notebooks                   <- Jupyter notebooks. Used for generation of additional figures and
    │                                  data exploration.
    │
    ├── reports
    │   └── figures                 <- Generated graphics and figures to be used in reporting.
    │
    ├── environment.yml             <- The conda requirements file for reproducing the analysis environment, e.g.
    │                                  generated with `conda env export > environment.yml`
    │                                  Usage: `conda env create -f environment.yml`.
    ├── requirements.txt            <- The requirements file for reproducing the analysis environment, e.g.
    │                                  generated with `pip freeze > requirements.txt`.
    │                                  Usage: `pip install -r requirements.txt`.
    │
    ├── setup.py                    <- makes project pip installable (pip install -e .) so src can be imported.
    ├── src                         <- Source code for use in this project.
    │   └── scripts                 <- Scripts to process data and train or evaluate models.
    │       │
    │       ├── data_scripts        <- Bash scripts to download VDJdb data. Used by `Makefile`.
    │       │
    │       ├── preprocessing       <- Python scripts to preprocess raw data. Used by `Makefile`.
    │       │                          Can also be run from the CLI with the `--help` flag for info.
    │       ├── hpc_scripts         <- Bash scripts to train models on HPC (or via bash on local machine).
    │       │
    │       ├── evaluate            <- Interactive Python scripts to visualize and evaluate trained models.
    │       │                          Can be run from the CLI with the `--help` flag for info.
    │       └── visualization       <- Scripts to create exploratory and results oriented visualizations
    │           └── visualize.sh    <- Contains Bash functions to utilise the different Python scripts.
    │
    └── tox.ini                     <- tox file with settings for running pytest.

## Usage

### Training

The `Makefile` contains the majority of steps required to reproduce our results, but some steps will need to be performed manually to avoid accidentally starting computationally-heavy operations.

- `make requirements`: creates a conda or virtualenv environment (can be done manually).
- `make test`: test if environment and source code passes all tests.
- `make data-vdjdb-aug-2019`: download data to correct directories.
- `preprocess-vdjdb-aug-2019`: filter the data and preprocess it into the correct format.
- For training models, choose a `.pbs` (bash) file in `./src/scripts/pbs_scripts` or modify its contents to the desired architecture/data combination. Must be run from within this directory in order for the output to be stored in the correct location and the script's paths to be correct.
- For evaluating models on external datasets, choose a `.sh` script in `./src/scripts/evaluate`. These files must be run from within this directory in order for the output to be stored in the correct location and the script's paths to be correct.
- `make metrics`: computes metrics and visualisations for all models that are found in the `models` directory.
- `metrics-compare`: shows a one-liner that can be utilised to compare two (or more) model directories, after they've been placed inside the same parent directory.
- `evaluate_self`: shows a one-liner that can be utilised to perform a per-epitope evaluation on a model directory.

Further evaluation and visualisation scripts are described in `./src/scripts/evaluate/visualization.sh`. There is also `./notebooks/figures.ipybn` which can be used to create more advanced combined figures.

### Predictions using the pre-built model

We included a final model that was trained using our methodology in the `models` directory. This model was trained on the VDJdb dataset ([August 2019 release](https://github.com/antigenomics/vdjdb-db/releases/tag/2019-08-08)) that was filtered on human TRB data, no 10x data and restricted to 10-20 (CDR3) or 8-11 (epitope) amino acid residues, with negatives that were generated by shuffling (i.e. sampling an negative epitope for each positive CDR3 sequence), and downsampling the most abundant epitopes down to 400 pairs each (`./models/pretrained/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001.h5`).

To use the model, use the `./src/scripts/predict/predict.py` Python script. It requires an input dataset in .csv format with a `;` separator. The model file can be changed from the default by using the `--model` flag. The other options can be shown via the `--help` flag.

Usage:

```
python ./src/scripts/predict/predict.py --model ./models/models-pretrained/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001.h5 --input input-sequences.csv --output output-predictions.csv
```

Below is an example input file to show the required format:

    cdr3;antigen.epitope
    CAVTDDKIIF;LLWNGPMAV
    CAVDSGGYQKVTF;LLWNGPMAV
    CAGGDDKIIF;LLWNGPMAV
    CAVKDARLMF;LLWNGPMAV
    CAVGSDKIIF;LLWNGPMAV
    CAVSDARLMF;LLWNGPMAV
    CAAFDDKIIF;LLWNGPMAV
    CAASASKLIF;LLWNGPMAV
    CAVDTNAGKSTF;LLWNGPMAV

### Results

An overview of the different experiments we performed in this project is provided in `.models/README.md`. Log files for these experiments can be found in `./src/scripts/hpc_scripts`. Please refer to the publication below for more information and the Zenodo data directory for a raw copy of the results.

## Zenodo data repository

The following Zenodo repository holds all output data (raw/processed datasets, trained model files in .h5 format, CV train/test splits, etc.): [10.5281/zenodo.3973547](https://doi.org/10.5281/zenodo.3973547).

## Citation

When using our results or modelling approach in a publication, please cite our preprint: [https://doi.org/10.1101/2019.12.18.880146](https://doi.org/10.1101/2019.12.18.880146):

>Moris, Pieter, Joey De Pauw, Anna Postovskaya, Benson Ogunjimi, Kris Laukens, and Pieter Meysman. 2019. "Current challenges for epitope-agnostic TCR interaction prediction and a new perspective derived from image classification". doi:10.1101/2019.12.18.880146.

## Authors

>Pieter Moris,
>Joey De Pauw,
>Wout Bittremieux,
>Anna Postovskaya,
>Benson Ogunjimi
>Kris Laukens,
>Pieter Meysman

## History of the project

A preprint of an earlier version of this work was submitted to bioRxiv: [https://doi.org/10.1101/2019.12.18.880146](https://doi.org/10.1101/2019.12.18.880146).

This project is an extension of original work done by [Joey De Pauw](https://github.com/JoeyDP) during his master thesis at the University of Antwerp under supervision of Pieter Meysman and Pieter Moris. A fork of the original project is made available at [https://github.com/pmoris/Master-Thesis](https://github.com/pmoris/Master-Thesis).

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
