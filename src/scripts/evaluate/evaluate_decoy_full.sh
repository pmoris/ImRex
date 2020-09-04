#!/bin/bash

export PROJECT_ROOT=$(readlink --canonicalize ../../../)

# # trbmhcidown - interaction map - shuffle - positives+newnegs

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type padded \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001.h5" \
# --name "train-true-fit-decoy-pos" \
# --data_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400-decoy.csv" \
# --neg_shuffle \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train.csv"

# # trbmhcidown - interaction map - shuffle - completetrain

# python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
# -i "${PROJECT_ROOT}/models/models-full/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train.csv" \
# -o "${PROJECT_ROOT}/models/models-full/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train-decoy.csv"

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type padded \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001.h5" \
# --name "train-true-fit-decoy-complete" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train-decoy.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train.csv"

# ################

# # trbmhcidown - interaction map - negref - positives+newnegs

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type padded \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001.h5" \
# --name "train-true-fit-decoy-pos" \
# --data_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400-decoy.csv" \
# --neg_ref "${PROJECT_ROOT}/data/raw/CDR3_control_sequences.tsv"  \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/iteration_no_validation/train.csv"

# # trbmhcidown - interaction map - negref - completetrain

# python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
# -i "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/iteration_no_validation/train.csv" \
# -o "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/iteration_no_validation/train-decoy.csv"

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type padded \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001.h5" \
# --name "train-true-fit-decoy-complete" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/iteration_no_validation/train-decoy.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/iteration_no_validation/train.csv"

# # trbmhcidown - interaction map - negref - traintrain

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type padded \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001.h5" \
# --name "train-true-fit-train" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/iteration_no_validation/train.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-21-20_trbmhcidown-negref-padded-b32-lre4-reg001/iteration_no_validation/train.csv"

# ################

# # trbmhcidown - dual - shuffle - positives+newnegs

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-41-31_trbmhcidown-shuffle-nettcr-b32-lre3-lrr/2020-07-24_19-41-31_trbmhcidown-shuffle-nettcr-b32-lre3-lrr.h5" \
# --name "train-true-fit-decoy-pos" \
# --data_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400-decoy.csv" \
# --neg_shuffle \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-41-31_trbmhcidown-shuffle-nettcr-b32-lre3-lrr/iteration_no_validation/train.csv"

# # trbmhcidown - dual - shuffle - completetrain

# python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
# -i "${PROJECT_ROOT}/models/models-full/2020-07-24_19-41-31_trbmhcidown-shuffle-nettcr-b32-lre3-lrr/iteration_no_validation/train.csv" \
# -o "${PROJECT_ROOT}/models/models-full/2020-07-24_19-41-31_trbmhcidown-shuffle-nettcr-b32-lre3-lrr/iteration_no_validation/train-decoy.csv"

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-41-31_trbmhcidown-shuffle-nettcr-b32-lre3-lrr/2020-07-24_19-41-31_trbmhcidown-shuffle-nettcr-b32-lre3-lrr.h5" \
# --name "train-true-fit-decoy-complete" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-24_19-41-31_trbmhcidown-shuffle-nettcr-b32-lre3-lrr/iteration_no_validation/train-decoy.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-41-31_trbmhcidown-shuffle-nettcr-b32-lre3-lrr/iteration_no_validation/train.csv"

# ################

# # trbmhcidown - dual - negref - positives+newnegs

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type padded \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-44-58_trbmhcidown-negref-nettcr-b32-lre3-lrr/2020-07-24_19-44-58_trbmhcidown-negref-nettcr-b32-lre3-lrr.h5" \
# --name "train-true-fit-decoy-pos" \
# --data_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400-decoy.csv" \
# --neg_ref "${PROJECT_ROOT}/data/raw/CDR3_control_sequences.tsv"  \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-44-58_trbmhcidown-negref-nettcr-b32-lre3-lrr/iteration_no_validation/train.csv"

# # trbmhcidown - dual - negref - completetrain

# python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
# -i "${PROJECT_ROOT}/models/models-full/2020-07-24_19-44-58_trbmhcidown-negref-nettcr-b32-lre3-lrr/iteration_no_validation/train.csv" \
# -o "${PROJECT_ROOT}/models/models-full/2020-07-24_19-44-58_trbmhcidown-negref-nettcr-b32-lre3-lrr/iteration_no_validation/train-decoy.csv"

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type padded \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-44-58_trbmhcidown-negref-nettcr-b32-lre3-lrr/2020-07-24_19-44-58_trbmhcidown-negref-nettcr-b32-lre3-lrr.h5" \
# --name "train-true-fit-decoy-complete" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-24_19-44-58_trbmhcidown-negref-nettcr-b32-lre3-lrr/iteration_no_validation/train-decoy.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-44-58_trbmhcidown-negref-nettcr-b32-lre3-lrr/iteration_no_validation/train.csv"

# ################

# # trbmhci - interaction map - shuffle - positives+newnegs

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type padded \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001/2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001.h5" \
# --name "train-true-fit-decoy-pos" \
# --data_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-decoy.csv" \
# --neg_shuffle \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train.csv"

# # trbmhci - interaction map - shuffle - completetrain

# python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
# -i "${PROJECT_ROOT}/models/models-full/2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train.csv" \
# -o "${PROJECT_ROOT}/models/models-full/2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train-decoy.csv"

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type padded \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001/2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001.h5" \
# --name "train-true-fit-decoy-complete" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train-decoy.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-30_11-30-27_trbmhci-shuffle-padded-b32-lre4-reg001/iteration_no_validation/train.csv"

# ################

# # trbmhci - dual lre4 - shuffle - positives+newnegs

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-30_11-34-13_trbmhci-shuffle-nettcr-b32-lre4-lrr/2020-07-30_11-34-13_trbmhci-shuffle-nettcr-b32-lre4-lrr.h5" \
# --name "train-true-fit-decoy-pos" \
# --data_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400-decoy.csv" \
# --neg_shuffle \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-30_11-34-13_trbmhci-shuffle-nettcr-b32-lre4-lrr/train.csv"

# # trbmhci - dual lre4 - shuffle - completetrain

# python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
# -i "${PROJECT_ROOT}/models/models-full/2020-07-30_11-34-13_trbmhci-shuffle-nettcr-b32-lre4-lrr/iteration_no_validation/train.csv" \
# -o "${PROJECT_ROOT}/models/models-full/2020-07-30_11-34-13_trbmhci-shuffle-nettcr-b32-lre4-lrr/iteration_no_validation/train-decoy.csv"

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-30_11-34-13_trbmhci-shuffle-nettcr-b32-lre4-lrr/2020-07-30_11-34-13_trbmhci-shuffle-nettcr-b32-lre4-lrr.h5" \
# --name "train-true-fit-decoy-complete" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-30_11-34-13_trbmhci-shuffle-nettcr-b32-lre4-lrr/iteration_no_validation/train-decoy.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-30_11-34-13_trbmhci-shuffle-nettcr-b32-lre4-lrr/iteration_no_validation/train.csv"

# ################

# # trbmhci - dual - shuffle - positives+newnegs

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-30_11-39-08_trbmhci-shuffle-nettcr-b32-lre3-lrr/2020-07-30_11-39-08_trbmhci-shuffle-nettcr-b32-lre3-lrr.h5" \
# --name "train-true-fit-decoy-pos" \
# --data_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400-decoy.csv" \
# --neg_shuffle \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-30_11-39-08_trbmhci-shuffle-nettcr-b32-lre3-lrr/train.csv"

# # trbmhci - dual - shuffle - completetrain

# python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
# -i "${PROJECT_ROOT}/models/models-full/2020-07-30_11-39-08_trbmhci-shuffle-nettcr-b32-lre3-lrr/iteration_no_validation/train.csv" \
# -o "${PROJECT_ROOT}/models/models-full/2020-07-30_11-39-08_trbmhci-shuffle-nettcr-b32-lre3-lrr/iteration_no_validation/train-decoy.csv"

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-30_11-39-08_trbmhci-shuffle-nettcr-b32-lre3-lrr/2020-07-30_11-39-08_trbmhci-shuffle-nettcr-b32-lre3-lrr.h5" \
# --name "train-true-fit-decoy-complete" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-30_11-39-08_trbmhci-shuffle-nettcr-b32-lre3-lrr/iteration_no_validation/train-decoy.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-30_11-39-08_trbmhci-shuffle-nettcr-b32-lre3-lrr/iteration_no_validation/train.csv"

# ################

# # trbmhcidown - dual lre4 - shuffle - positives+newnegs

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-29-13_trbmhcidown-shuffle-nettcr-b32-lre4-lrr/2020-07-24_19-29-13_trbmhcidown-shuffle-nettcr-b32-lre4-lrr.h5" \
# --name "train-true-fit-decoy-pos" \
# --data_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400-decoy.csv" \
# --neg_shuffle \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-29-13_trbmhcidown-shuffle-nettcr-b32-lre4-lrr/iteration_no_validation/train.csv"

# # trbmhcidown - dual lre4 - shuffle - completetrain

# python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
# -i "${PROJECT_ROOT}/models/models-full/2020-07-24_19-29-13_trbmhcidown-shuffle-nettcr-b32-lre4-lrr/iteration_no_validation/train.csv" \
# -o "${PROJECT_ROOT}/models/models-full/2020-07-24_19-29-13_trbmhcidown-shuffle-nettcr-b32-lre4-lrr/iteration_no_validation/train-decoy.csv"

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-29-13_trbmhcidown-shuffle-nettcr-b32-lre4-lrr/2020-07-24_19-29-13_trbmhcidown-shuffle-nettcr-b32-lre4-lrr.h5" \
# --name "train-true-fit-decoy-complete" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-24_19-29-13_trbmhcidown-shuffle-nettcr-b32-lre4-lrr/iteration_no_validation/train-decoy.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-29-13_trbmhcidown-shuffle-nettcr-b32-lre4-lrr/iteration_no_validation/train.csv"

# ################

# # trbmhcidown - dual - negref - positives+newnegs

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-32-26_trbmhcidown-negref-nettcr-b32-lre4-lrr/2020-07-24_19-32-26_trbmhcidown-negref-nettcr-b32-lre4-lrr.h5" \
# --name "train-true-fit-decoy-pos" \
# --data_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human-trb-mhci-no10x-size-down400-decoy.csv" \
# --neg_ref "${PROJECT_ROOT}/data/raw/CDR3_control_sequences.tsv"  \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-32-26_trbmhcidown-negref-nettcr-b32-lre4-lrr/iteration_no_validation/train.csv"

# # trbmhcidown - dual - negref - completetrain

# python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
# -i "${PROJECT_ROOT}/models/models-full/2020-07-24_19-32-26_trbmhcidown-negref-nettcr-b32-lre4-lrr/iteration_no_validation/train.csv" \
# -o "${PROJECT_ROOT}/models/models-full/2020-07-24_19-32-26_trbmhcidown-negref-nettcr-b32-lre4-lrr/iteration_no_validation/train-decoy.csv"

# python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
# --model_type separated \
# --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
# --model "${PROJECT_ROOT}/models/models-full/2020-07-24_19-32-26_trbmhcidown-negref-nettcr-b32-lre4-lrr/2020-07-24_19-32-26_trbmhcidown-negref-nettcr-b32-lre4-lrr.h5" \
# --name "train-true-fit-decoy-complete" \
# --data_path "${PROJECT_ROOT}/models/models-full/2020-07-24_19-32-26_trbmhcidown-negref-nettcr-b32-lre4-lrr/iteration_no_validation/train-decoy.csv" \
# --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
# --per_epitope \
# --train_dataset "${PROJECT_ROOT}/models/models-full/2020-07-24_19-32-26_trbmhcidown-negref-nettcr-b32-lre4-lrr/iteration_no_validation/train.csv"




function evaluate_decoy {
    # readlink -f $1 ;
    DIR=$(realpath -s $1)

    MODEL=$(basename ${DIR})

    DATA=$(basename "$(sed -n "s/^.*CLI argument data_path: '\(.*\)'.*$/\1/p" ${DIR}/${MODEL}.log)")
    DATA=${DATA/.csv/-decoy.csv}
    DATA="${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/${DATA}"

    echo "Creating decoy dataset based on train dataset ${DIR}/iteration_no_validation/train.csv."
    python "${PROJECT_ROOT}/src/scripts/preprocessing/decoy_epitopes.py" \
        -i "${DIR}/iteration_no_validation/train.csv" \
        -o "${DIR}/iteration_no_validation/train-decoy.csv"

    if [[ ${MODEL} == *"negref"* ]]; then
        echo "Found model built with reference negatives."
        NEG_TYPE="--neg_ref ${PROJECT_ROOT}/data/raw/CDR3_control_sequences.tsv"

    elif [[ ${MODEL} == *"shuffle"* ]]; then
        echo "Found model built with shuffled negatives."
        NEG_TYPE="--neg_shuffle"
    fi

    echo $NEG_TYPE

    if [[ ${MODEL} == *"nettcr"* ]]; then
        echo "Found model built with dual inputs."

        echo "Evaluate on decoy dataset derived from the complete train dataset (including negatives)."

        python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
        --model_type separated \
        --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
        --model "${DIR}/${MODEL}.h5" \
        --name "train-true-fit-decoy-complete" \
        --data_path "${DIR}/iteration_no_validation/train-decoy.csv" \
        --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
        --per_epitope \
        --train_dataset "${DIR}/iteration_no_validation/train.csv"

        # echo "Evaluate on decoy dataset derived from the positive train dataset (generate new negatives)."

        # python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
        # --model_type separated \
        # --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \
        # --model "${DIR}/${MODEL}.h5" \
        # --name "train-true-fit-decoy-pos" \
        # --data_path "${DATA}" \
        # ${NEG_TYPE} \
        # --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
        # --per_epitope \
        # --train_dataset "${DIR}/iteration_no_validation/train.csv"

    elif [[ ${MODEL} == *"padded"* ]]; then
        echo "Found model built with interaction map input."

        echo "Evaluate on decoy dataset derived from the complete train dataset (including negatives)."

        python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
        --model_type padded \
        --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \
        --model "${DIR}/${MODEL}.h5" \
        --name "train-true-fit-decoy-complete" \
        --data_path "${DIR}/iteration_no_validation/train-decoy.csv" \
        --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
        --per_epitope \
        --train_dataset "${DIR}/iteration_no_validation/train.csv"

        # echo "Evaluate on decoy dataset derived from the positive train dataset (generate new negatives)."

        # python "${PROJECT_ROOT}/src/scripts/evaluate/evaluate.py" \
        # --model_type padded \
        # --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \
        # --model "${DIR}/${MODEL}.h5" \
        # --name "train-true-fit-decoy-pos" \
        # --data_path "${DATA}" \
        # ${NEG_TYPE} \
        # --full_dataset_path "${PROJECT_ROOT}/data/interim/vdjdb-2019-08-08/vdjdb-human.csv" \
        # --per_epitope \
        # --train_dataset "${DIR}/iteration_no_validation/train.csv"
    fi
}

export -f evaluate_decoy

find "${PROJECT_ROOT}/models/models-full/" -maxdepth 1 -mindepth 1 -type d -exec /bin/bash -c 'evaluate_decoy "$0"' {} \;
