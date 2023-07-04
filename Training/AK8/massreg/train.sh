#!/bin/bash

set -x

echo "args: $@"

# set the dataset dir
[[ -z $DATADIR ]] && DATADIR='/data/adlintul/input/ak8'

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

epochs=1
samples_per_epoch=$((10000 * 1024 / $NGPUS))
samples_per_epoch_val=$((10000 * 128))
dataopts="--num-workers 4 --fetch-step 0.01"
batchopts="--batch-size 256 --start-lr 5e-3"

$CMD \
    --regression-mode \
    --data-train \
    "BulkGravitonToHHTo4Q:${DATADIR}/BulkGravitonToHH_*/Run3Summer22EE/*/0000/ntuple_1.root" \
    "QCD:${DATADIR}/QCD_*/Run3Summer22EE/*/0000/ntuple_1.root" \
    --data-config ${data_config} --network-config ${model_config} \
    --model-prefix ${OUTPUT}/ak8_massreg_{auto}${suffix}/net \
    --gpus "2" \
    $dataopts $batchopts \
    --samples-per-epoch ${samples_per_epoch} --samples-per-epoch-val ${samples_per_epoch_val} \
    --num-epochs ${epochs} \
    --optimizer ranger --log ${OUTPUT}/ak8_massreg_{auto}${suffix}.log --predict-output pred.root \
    "${@:2}"
