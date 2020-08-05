#!/bin/bash

export PROJECT_ROOT=$(readlink --canonicalize ../../../../) # old method

python ${PROJECT_ROOT}/src/scripts/train/scenario_padding.py \
--batch_size 32 \
--epochs 20 \
--cv epitope_grouped \
--n_folds 5 \
--full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
--min_length_cdr3 10 \
--max_length_cdr3 20 \
--min_length_epitope 8 \
--max_length_epitope 11 \
--features "hydrophob,isoelectric,mass,hydrophil" \
--operator "absdiff" \
--data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
--model model_padded \
--depth1_1 128 \
--depth1_2 64 \
--depth2_1 128 \
--depth2_2 64 \
--activation_function_conv "relu" \
--activation_function_dense "relu" \
--dropout_conv 0.25 \
--learning_rate 0.0001 \
--regularization 0.01 \
--optimizer rmsprop \
--epitope_ratio \
--name "trbmhcidown-epitope_grouped-shuffle-padded-b32-lre4-reg001-local-epitoperatio" 2>&1 | tee -a "trbmhcidown-epitope_grouped-shuffle-padded-b32-lre4-reg001-local-epitoperatio.log"


python ${PROJECT_ROOT}/src/scripts/train/scenario_padding.py \
--batch_size 32 \
--epochs 20 \
--cv epitope_grouped \
--n_folds 5 \
--full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
--min_length_cdr3 10 \
--max_length_cdr3 20 \
--min_length_epitope 8 \
--max_length_epitope 11 \
--features "hydrophob,isoelectric,mass,hydrophil" \
--operator "absdiff" \
--data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
--model model_padded \
--depth1_1 128 \
--depth1_2 64 \
--depth2_1 128 \
--depth2_2 64 \
--activation_function_conv "relu" \
--activation_function_dense "relu" \
--dropout_conv 0.25 \
--learning_rate 0.0001 \
--regularization 0.01 \
--optimizer rmsprop \
--epitope_ratio \
--name "trbmhcidown-epitope_grouped-shuffle-padded-b32-lre4-reg001-local-epitoperatio" 2>&1 | tee -a "trbmhcidown-epitope_grouped-shuffle-padded-b32-lre4-reg001-local-epitoperatio.log"


python ${PROJECT_ROOT}/src/scripts/train/scenario_padding.py \
--batch_size 32 \
--epochs 20 \
--cv epitope_grouped \
--n_folds 5 \
--full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
--min_length_cdr3 10 \
--max_length_cdr3 20 \
--min_length_epitope 8 \
--max_length_epitope 11 \
--features "hydrophob,isoelectric,mass,hydrophil" \
--operator "absdiff" \
--data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
--model model_padded \
--depth1_1 128 \
--depth1_2 64 \
--depth2_1 128 \
--depth2_2 64 \
--activation_function_conv "relu" \
--activation_function_dense "relu" \
--dropout_conv 0.25 \
--learning_rate 0.0001 \
--regularization 0.01 \
--optimizer rmsprop \
--name "trbmhcidown-epitope_grouped-shuffle-padded-b32-lre4-reg001-local" 2>&1 | tee -a "trbmhcidown-epitope_grouped-shuffle-padded-b32-lre4-reg001-local.log"


python ${PROJECT_ROOT}/src/scripts/train/scenario_padding.py \
--batch_size 32 \
--epochs 20 \
--cv epitope_grouped \
--n_folds 5 \
--full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
--min_length_cdr3 10 \
--max_length_cdr3 20 \
--min_length_epitope 8 \
--max_length_epitope 11 \
--features "hydrophob,isoelectric,mass,hydrophil" \
--operator "absdiff" \
--data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
--model model_padded \
--depth1_1 128 \
--depth1_2 64 \
--depth2_1 128 \
--depth2_2 64 \
--activation_function_conv "relu" \
--activation_function_dense "relu" \
--dropout_conv 0.25 \
--learning_rate 0.0001 \
--regularization 0.01 \
--optimizer rmsprop \
--name "trbmhcidown-epitope_grouped-shuffle-padded-b32-lre4-reg001-local" 2>&1 | tee -a "trbmhcidown-epitope_grouped-shuffle-padded-b32-lre4-reg001-local.log"


python ${PROJECT_ROOT}/src/scripts/train/scenario_padding.py \
--batch_size 32 \
--epochs 20 \
--cv epitope_grouped_shuffle \
--n_folds 5 \
--full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
--min_length_cdr3 10 \
--max_length_cdr3 20 \
--min_length_epitope 8 \
--max_length_epitope 11 \
--features "hydrophob,isoelectric,mass,hydrophil" \
--operator "absdiff" \
--data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
--model model_padded \
--depth1_1 128 \
--depth1_2 64 \
--depth2_1 128 \
--depth2_2 64 \
--activation_function_conv "relu" \
--activation_function_dense "relu" \
--dropout_conv 0.25 \
--learning_rate 0.0001 \
--regularization 0.01 \
--optimizer rmsprop \
--name "trbmhcidown-epitope_grouped-shufflesplit-shuffle-padded-b32-lre4-reg001-local" 2>&1 | tee -a "trbmhcidown-epitope_grouped-shufflesplit-shuffle-padded-b32-lre4-reg001-local.log"


python ${PROJECT_ROOT}/src/scripts/train/scenario_padding.py \
--batch_size 32 \
--epochs 20 \
--cv epitope_grouped_shuffle \
--n_folds 5 \
--full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
--min_length_cdr3 10 \
--max_length_cdr3 20 \
--min_length_epitope 8 \
--max_length_epitope 11 \
--features "hydrophob,isoelectric,mass,hydrophil" \
--operator "absdiff" \
--data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
--model model_padded \
--depth1_1 128 \
--depth1_2 64 \
--depth2_1 128 \
--depth2_2 64 \
--activation_function_conv "relu" \
--activation_function_dense "relu" \
--dropout_conv 0.25 \
--learning_rate 0.0001 \
--regularization 0.01 \
--optimizer rmsprop \
--name "trbmhcidown-epitope_grouped-shufflesplit-shuffle-padded-b32-lre4-reg001-local" 2>&1 | tee -a "trbmhcidown-epitope_grouped-shuffle-padded-b32-lre4-reg001-local.log"


# python ${PROJECT_ROOT}/src/scripts/train/scenario_padding.py \
# --batch_size 32 \
# --epochs 20 \
# --cv epitope_grouped \
# --n_folds 5 \
# --neg_ref "${PROJECT_ROOT}/data/raw/CDR3_control_sequences.tsv" \
# --min_length_cdr3 10 \
# --max_length_cdr3 20 \
# --min_length_epitope 8 \
# --max_length_epitope 11 \
# --features "hydrophob,isoelectric,mass,hydrophil" \
# --operator "absdiff" \
# --data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
# --model model_padded \
# --depth1_1 128 \
# --depth1_2 64 \
# --depth2_1 128 \
# --depth2_2 64 \
# --activation_function_conv "relu" \
# --activation_function_dense "relu" \
# --dropout_conv 0.25 \
# --learning_rate 0.0001 \
# --regularization 0.01 \
# --optimizer rmsprop \
# --name "trbmhcidown-epitope_grouped-negref-padded-b32-lre4-reg001-local" 2>&1 | tee -a "trbmhcidown-epitope_grouped-negref-padded-b32-lre4-reg001-local.log"






# python ${PROJECT_ROOT}/src/scripts/train/scenario_separated_inputs.py \
# --batch_size 32 \
# --epochs 20 \
# --cv epitope_grouped \
# --n_folds 5 \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --min_length_cdr3 10 \
# --max_length_cdr3 20 \
# --min_length_epitope 8 \
# --max_length_epitope 11 \
# --data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
# --model nettcr \
# --learning_rate 0.0001 \
# --optimizer rmsprop \
# --name "trbmhcidown-epitope_grouped-shuffle-nettcr-b32-lre4-local" 2>&1 | tee -a "trbmhcidown-epitope_grouped-shuffle-nettcr-b32-lre4-local.log"


# python ${PROJECT_ROOT}/src/scripts/train/scenario_separated_inputs.py \
# --batch_size 32 \
# --epochs 20 \
# --cv epitope_grouped \
# --n_folds 5 \
# --neg_ref "${PROJECT_ROOT}/data/raw/CDR3_control_sequences.tsv" \
# --min_length_cdr3 10 \
# --max_length_cdr3 20 \
# --min_length_epitope 8 \
# --max_length_epitope 11 \
# --data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
# --model nettcr \
# --learning_rate 0.0001 \
# --optimizer rmsprop \
# --name "trbmhcidown-epitope_grouped-negref-nettcr-b32-lre4-local" 2>&1 | tee -a "trbmhcidown-epitope_grouped-negref-nettcr-b32-lre4-local.log"






# python ${PROJECT_ROOT}/src/scripts/train/scenario_separated_inputs.py \
# --batch_size 32 \
# --epochs 20 \
# --cv epitope_grouped \
# --n_folds 5 \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --min_length_cdr3 10 \
# --max_length_cdr3 20 \
# --min_length_epitope 8 \
# --max_length_epitope 11 \
# --data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
# --model nettcr \
# --learning_rate 0.001 \
# --optimizer rmsprop \
# --name "trbmhcidown-epitope_grouped-shuffle-nettcr-b32-lre3-local" 2>&1 | tee -a "trbmhcidown-epitope_grouped-shuffle-nettcr-b32-lre3-local.log"


# python ${PROJECT_ROOT}/src/scripts/train/scenario_separated_inputs.py \
# --batch_size 32 \
# --epochs 20 \
# --cv epitope_grouped \
# --n_folds 5 \
# --neg_ref "${PROJECT_ROOT}/data/raw/CDR3_control_sequences.tsv" \
# --min_length_cdr3 10 \
# --max_length_cdr3 20 \
# --min_length_epitope 8 \
# --max_length_epitope 11 \
# --data_path ${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down.csv \
# --model nettcr \
# --learning_rate 0.001 \
# --optimizer rmsprop \
# --name "trbmhcidown-epitope_grouped-negref-nettcr-b32-lre3-local" 2>&1 | tee -a "trbmhcidown-epitope_grouped-negref-nettcr-b32-lre3-local.log"
