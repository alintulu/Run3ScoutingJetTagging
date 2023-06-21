#!/bin/bash

#set -x

echo "args: $@"

# set the dataset dir
[[ -z $DATADIR ]] && DATADIR='/eos/cms/store/group/ml/Tagging4ScoutingHackathon/Adelina/hbb/ak8'

# set the dataset dir
[[ -z $OUTPUT ]] && OUTPUT='output'

# set a comment via `COMMENT`
suffix=${COMMENT}

# set the number of gpus for DDP training via `DDP_NGPUS`
NGPUS=${DDP_NGPUS}
[[ -z $NGPUS ]] && NGPUS=1
if ((NGPUS > 1)); then
    CMD="torchrun --standalone --nnodes=1 --nproc_per_node=$NGPUS $(which weaver) --backend nccl"
else
    CMD="weaver"
fi

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
    --gpus "0" \
    --num-epochs 2 \
    --optimizer ranger --log ${OUTPUT}/ak8_massreg_{auto}${suffix}.log --predict-output pred.root \
    "${@:2}"
