#!/bin/bash

export PROJECT_ROOT=$(readlink --canonicalize ../../../) # old method

python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" --model_type padded --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001.h5"  --name "adaptive-sars-cov" --data_path "${PROJECT_ROOT}/data/interim/immunecode-adaptive/adaptive-sars-cov.csv" --neg_shuffle --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" --per_epitope --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train.csv"
