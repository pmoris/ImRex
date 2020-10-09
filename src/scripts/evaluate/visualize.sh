#!/bin/bash

function evaluate_self_plots {
    echo $1
    if [[ $1 == *"grouped"* ]]; then
        if [[ $1 == *"decoy"* ]]; then
            python ../../src/scripts/evaluate/visualize.py evaluate_self_plots --min_obs 30 --grouped True --decoy True "$1"
        else
            python ../../src/scripts/evaluate/visualize.py evaluate_self_plots --min_obs 30 --grouped True "$1"
        fi
    else
        if [[ $1 == *"decoy"* ]]; then
            python ../../src/scripts/evaluate/visualize.py evaluate_self_plots --min_obs 30 --min_iterations 25 --decoy True "$1"
        else
            python ../../src/scripts/evaluate/visualize.py evaluate_self_plots --min_obs 30 --min_iterations 25 "$1"
        fi
    fi
}
export -f evaluate_self_plots

function evaluate_self_comparison_plots {
    echo $1
    if [[ $1 == *"grouped"* ]]; then
        if [[ $1 == *"decoy"* ]]; then
            python ../../src/scripts/evaluate/visualize.py evaluate_self_comparison_plots --min_obs 30 --grouped True --decoy True "$1"
        else
            python ../../src/scripts/evaluate/visualize.py evaluate_self_comparison_plots --min_obs 30 --grouped True "$1"
        fi
    else
        if [[ $1 == *"decoy"* ]]; then
            python ../../src/scripts/evaluate/visualize.py evaluate_self_comparison_plots --min_obs 30 --min_iterations 25 --decoy True "$1"
        else
            python ../../src/scripts/evaluate/visualize.py evaluate_self_comparison_plots --min_obs 30 --min_iterations 25 "$1"
        fi
    fi
}
export -f evaluate_self_comparison_plots

# create metrics: requires a directory with subdirectories that contain "iteration_#" subdirectories
# find . -maxdepth 2 -mindepth 2 -type d -exec python ../../src/scripts/evaluate/visualize.py metrics --force True "{}" \;
# should already be done

# evaluate self: create predictions for test folds
# find . -maxdepth 1 -mindepth 1 -type d -name "*nettcr*" -exec python ../../src/scripts/evaluate/evaluate_self.py --input {} --model_type separated --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 \;
# find . -maxdepth 1 -mindepth 1 -type d -name "*padded*" -exec python ../../src/scripts/evaluate/evaluate_self.py --input {} --model_type padded --min_length_cdr3 10 --max_length_cdr3 20 --min_length_epitope 8 --max_length_epitope 11 --features 'hydrophob,isoelectric,mass,hydrophil' --operator absdiff \;
# should already be done

# create comparisons plots: requires a directory with subdirectories containing pre-computed iteration directories and their consolidated metrics
find . -maxdepth 1 -mindepth 1 -type d -exec python ../../src/scripts/evaluate/visualize.py compare --force True "{}" \;

# # create per-epitope plots for all "evaluate_test_folds" subdirectories
find . -maxdepth 3 -mindepth 3 -type d -name "*evaluate_test_folds" -exec /bin/bash -c 'evaluate_self_plots "$0"' {} \;

# create comparison plots for per-epitope comparisons: same requirements as above + subdirectories should contain an "evaluate_test_folds" directory with metrics_per_epitope.csv
find . -maxdepth 1 -mindepth 1 -type d -exec /bin/bash -c 'evaluate_self_comparison_plots "$0"' {} \;

# Quicker alternative one-liner that can be used inside a directory with subdirectories that need to be compared
# find . -maxdepth 1 -mindepth 1 -type d -exec python ../../../../src/scripts/evaluate/visualize.py metrics --force True "{}" \; && python ../../../../src/scripts/evaluate/visualize.py compare . --force True && python ../../../../src/scripts/evaluate/visualize.py evaluate_self_comparison_plots . --min_obs 30 --min_iterations 25 && find . -maxdepth 2 -mindepth 2 -type d -name "*evaluate_test_folds" -exec python ../../../../src/scripts/evaluate/visualize.py evaluate_self_plots --min_obs 30 --min_iterations 25 "{}" \;

# find . -maxdepth 1 -mindepth 1 -type d -exec python ../../../../src/scripts/evaluate/visualize.py metrics --force True "{}" \; && python ../../../../src/scripts/evaluate/visualize.py compare . --force True && python ../../../../src/scripts/evaluate/visualize.py evaluate_self_comparison_plots . --min_obs 30 --grouped True && find . -maxdepth 2 -mindepth 2 -type d -name "*evaluate_test_folds" -exec python ../../../../src/scripts/evaluate/visualize.py evaluate_self_plots --min_obs 30 --grouped True "{}" \;

# add --decoy True when comparing normal and decoy models
# add --grouped True for epitope grouped models

###################################################################

# diff <(sort ) <(sort )
# diff <(sort Interaction\ map\ -\ TRB\ -\ epitope-grouped\ -\ shuffled\ negatives/iteration_0/train_fold_0.csv) <(sort Dual\ input\ -\ TRB\ -\ epitope-grouped\ -\ shuffled\ negatives/iteration_0/train_fold_0.csv)
# comm -1 -2 ...
# diff <(awk -F ';' '$1 == 1 { print $1, $2, $3}' ./models/models-calcua-final/models-nettcrbn-b32-lre4-lrr/2020-06-05_23-49-12_trbmhci-repeated5fold-negref-nettcrbn-b32-lre4-lrr/iteration_0/train_fold_0.csv | sort) <(awk -F ';' '$1 == 1 { print $1, $2, $3}' ./models/models-calcua-final/nettcr/2020-06-03_00-58-41_trbmhci-repeated5fold-negref-nettcrrelu-b32-lre3-lrr/iteration_0/train_fold_0.csv | sort)

# check if two different model directories used the same postiive/negative data
# diff <(awk -F ';' '$1 == 1 { print $1, $2, $3}' ./models/models-calcua-final/models-groupedseed/2020-06-11_19-56-25_trbmhcidown-epitope_grouped-shuffle-nettcr-b32-lre3-lrr/iteration_0/train_fold_0.csv | sort) <(awk -F ';' '$1 == 1 { print $1, $2, $3}' ./models/models-calcua-final/nettcr/2020-06-03_12-54-42_trbmhcidown-epitope_grouped-shuffle-nettcrrelu-b32-lre3-lrr/iteration_0/train_fold_0.csv | sort)
# empty
