#!/bin/bash
export NCPUS=32
export MCCE_EST_EXP="PrefDCE"
export MCCE_MIN_RT=0.35
export MCCE_MAX_RT=10
export MCCE_CONTAM=0.02
export MCCE_METHOD="model"
# Provide tag on the command line
if [ -z "$1" ]
then
    echo "No tag argument supplied"
    echo "job_scripts/PrefEstimateAWS.sh <tag> <model> <input_data>"
    exit 1
fi
export MCCE_TAG="$1"
if [ -z "$2" ]
then
    echo "No model argument supplied"
    echo "job_scripts/PrefEstimateAWS.sh <tag> <model> <input_data>"
    exit 1
fi
export MCCE_MODEL="$2"

if [ -z "$3" ]
then
    echo "No input data file argument supplied"
    echo "job_scripts/PrefEstimateAWS.sh <tag> <model> <input_data>"
    exit 1
fi
export MCCE_EXP_DATA="$3"

if [ ! -f "src/model_estimation.R" ]
then
    echo "Command should be run from project root directory"
    echo "job_scripts/PrefEstimateAWS.sh <tag> <model> <input_data>"
    exit 1
fi

# Required for gsl (used by rtdists), due to our installation method.
export LD_LIBRARY_PATH="/usr/local/gsl/2.5/x86_64/lib64"

Rscript --no-save --no-restore src/model_estimation.R > estimation_$MCCE_EST_EXP.$MCCE_MODEL.$MCCE_TAG.out 2> estimation_$MCCE_EST_EXP.$MCCE_MODEL.$MCCE_TAG.err

outdir=$(find ~/studies -name "egress*")

cp estimation*$MCCE_TAG* $outdir
cp data/output/PrefDCE*$MCCE_TAG* $outdir

sleep 60

sudo shutdown -h now
