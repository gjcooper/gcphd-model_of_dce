#!/bin/bash
export NCPUS=32
export DCE_EST_EXP="PrefDCE"
export DCE_EXP_DATA="Pref_preprocessed.RDS"
export DCE_MIN_RT=0.35
export DCE_MAX_RT=10
export DCE_CONTAM=0.02
export VDCE_TAG="staged_burn"

# Required for gsl (used by rtdists), due to our installation method.
export LD_LIBRARY_PATH="/usr/local/gsl/2.5/x86_64/lib64"

Rscript --no-save --no-restore src/reduced_estimation.R > estimation_$DCE_EST_EXP.$VDCE_TAG.out 2> estimation_$DCE_EST_EXP.$VDCE_TAG.err

cp PrefReducedAWS.sh ~/studies/egress-store-189fa174-3bd8-40c0-a722-d93c41c19ea0/.
cp estimation* ~/studies/egress-store-189fa174-3bd8-40c0-a722-d93c41c19ea0/.
cp data/output/PrefDCE* ~/studies/egress-store-189fa174-3bd8-40c0-a722-d93c41c19ea0/.

sleep 60

sudo shutdown -h now
