#!/bin/bash
export NCPUS=1
export MCCE_EST_EXP="NumericVDCE"
export MCCE_MIN_RT=0.35
export MCCE_MAX_RT=10
export MCCE_CONTAM=0.02
export MCCE_MODEL="std"
export MCCE_TAG="local_profile"
export MCCE_METHOD="profile"
export MCCE_EXP_DATA="Task1_preprocessed_Accept.RDS"

if [ ! -f "src/model_estimation.R" ]
then
    echo "Command should be run from project root directory"
    echo "job_scripts/PrefEstimateAWS.sh <tag> <model> <input_data>"
    exit 1
fi

Rscript --no-save --no-restore src/model_estimation.R > estimation_$MCCE_EST_EXP.$MCCE_MODEL.$MCCE_TAG.out 2> estimation_$MCCE_EST_EXP.$MCCE_MODEL.$MCCE_TAG.err
