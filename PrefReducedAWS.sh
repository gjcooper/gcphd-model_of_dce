#!/bin/bash
export NCPUS=32
export DCE_EST_EXP="PrefDCE"
export DCE_EXP_DATA="Pref_preprocessed.RDS"
export DCE_MIN_RT=0.35
export DCE_MAX_RT=10
export DCE_CONTAM=0.02
export VDCE_TAG="reduced"

# Required for gsl (used by rtdists), due to our installation method.
export LD_LIBRARY_PATH="/usr/local/gsl/2.5/x86_64/lib64"

function notify() {
  message="${1}"
  instance_identity="$(curl http://169.254.169.254/latest/dynamic/instance-identity/document)"
  account_id="$(echo "${instance_identity}" | jq -r '.accountId')"
  region="$(echo "${instance_identity}" | jq -r '.region')"
  sns_topic_name="uonswpoc-alerts"  # TODO: Retrieve this name from somewhere, potentially CFN exports or tags?
  topic_arn="arn:aws:sns:${region}:${account_id}:${sns_topic_name}"

  echo "${message}" && echo
  aws sns publish --message "${message}" --region "${region}" --topic-arn "${topic_arn}"
}

Rscript --no-save --no-restore src/reduced_estimation.R > estimation_$DCE_EST_EXP.$VDCE_TAG.out 2> estimation_$DCE_EST_EXP.$VDCE_TAG.err

cp PrefReducedAWS.sh ~/studies/egress-store-189fa174-3bd8-40c0-a722-d93c41c19ea0/.
cp estimation* ~/studies/egress-store-189fa174-3bd8-40c0-a722-d93c41c19ea0/.
cp data/output/PrefDCE* ~/studies/egress-store-189fa174-3bd8-40c0-a722-d93c41c19ea0/.

sleep 60

sudo shutdown -h now