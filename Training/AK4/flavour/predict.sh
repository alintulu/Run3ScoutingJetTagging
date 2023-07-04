#!/bin/bash

# set the dataset dir
[[ -z $DATADIR ]] && DATADIR='/data/adlintul/input/ak4'

# set the dataset dir
[[ -z $OUTPUT ]] && OUTPUT='output'

# set the model directory
[[ -z $MODEL ]] && MODEL='output'

data_config="data/flavour.yaml"
model_config="networks/flavour.py"

weaver \
	--predict \
      	--data-test "${DATADIR}/BulkGravitonToHH_MX960_MH121_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EE/230322_202038/0000/ntuple_1.root" \
	--data-config ${data_config} --network-config ${model_config} \
	--model-prefix ${MODEL}/net_best_epoch_state.pt \
	--gpus '3' --batch-size 512 \
	--predict-output ${OUTPUT}/pred.root
