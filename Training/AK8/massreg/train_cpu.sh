#!/bin/bash

#set -x

echo "args: $@"

# set the dataset dir
[[ -z $DATADIR ]] && DATADIR='/eos/cms/store/group/ml/Tagging4ScoutingHackathon/Adelina/hbb/ak8'

# set the dataset dir
[[ -z $OUTPUT ]] && OUTPUT='output'

# set a comment via `COMMENT`
suffix=${COMMENT}

# train on CPU
CMD="weaver"

data_config="data/massreg.yaml"
model_config="networks/massreg.py"

$CMD \
    --regression-mode \
    --demo \
    --data-train \
    "BulkGravitonToHHTo4Q:${DATADIR}/BulkGravitonToHH_MX960_MH82_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EE/230322_202312/0000/ntuple_1*.root" \
    "QCD:${DATADIR}/QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8/Run3Summer22EE/230328_083632/0000/ntuple_1*.root" \
    --data-config ${data_config} --network-config ${model_config} \
    --model-prefix ${OUTPUT}/ak8_massreg_{auto}${suffix}/net \
    --gpus "" \
    --num-epochs 1 \
    --optimizer ranger --log ${OUTPUT}/ak8_massreg_{auto}${suffix}.log --predict-output pred.root \
    "${@:2}"
