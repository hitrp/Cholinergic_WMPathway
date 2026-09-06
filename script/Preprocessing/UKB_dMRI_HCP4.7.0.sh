#!/bin/bash

STUDYFOLDER=$1
SUBJID=$2

DIR_SH=tools/HCP_Pipelines/HCPpipelines-4.7.0/Examples/Scripts
outdir=${STUDYFOLDER}/${SUBJID}

# Detect --no-gpu in arguments and forward when present
NO_GPU_FLAG="--no-gpu"
for arg in "${@}"; do
	if [ "${arg}" = "--no-gpu" ]; then
		NO_GPU_FLAG="--no-gpu"
		break
	fi
done

# diffusion data
$DIR_SH/UKB_DiffusionPreprocessingBatch.sh --StudyFolder=$STUDYFOLDER --Subject=$SUBJID --runlocal ${NO_GPU_FLAG}
