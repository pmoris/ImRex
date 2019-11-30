#!/bin/bash


for FEATURE in mutationstab
do
    for OPERATOR in prod diff absdiff layer
    do
        echo "$FEATURE - $OPERATOR"
        qsub -F "padding --features $FEATURE --operator $OPERATOR" -N "${FEATURE}_${OPERATOR}" run.pbs
    done
done


