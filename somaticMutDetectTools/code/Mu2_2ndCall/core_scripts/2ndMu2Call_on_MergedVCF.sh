#!/bin/bash

BamList=$(realpath $1)
TagList=$2
Spe=$3
timeID=$(timeID)

mkdir -p FilMtx_From2ndMu2Call
cp ${TagList} FilMtx_From2ndMu2Call/TagList && cd FilMtx_From2ndMu2Call
cat TagList | xargs -i BsubS -n 1 -M 6 "ForceCallonMergedVCF_SepChr_MedQue.sh {} ${BamList} UnfilteredVCF_on_1stMu2Call ${Spe} ${timeID}" Mu2ForceCall >> JobID_Mu2ForceCall_SubmitOnColonies
