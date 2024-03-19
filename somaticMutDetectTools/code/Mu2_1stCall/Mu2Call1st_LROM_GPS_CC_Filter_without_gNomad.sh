#!/bin/bash

# This is based on PIPELINE_codes_Mutect2_Call_LROM_GPS_CC_Filter_without_gNomad

module load openjdk/11.0.5-10
#module load jdk/11.0.4

RefFa=/rsrch3/scratch/reflib/iacs/reference/human/hg19/fastas/Homo_sapiens_assembly19.fasta
TimeID=$(date +%Y%m%d%H%M%S)

JavaOption='--java-options "-Xmx40G"'

# $1=10-PLATE1C1
BamList=$1
Tag=$2

BamPath=$(cat ${BamList} | grep ${Tag})

# Confirm the bamlist called is one and only.
if [[ ! $(cat ${BamList} | grep ${Tag} | wc -l) -eq 1 ]]
then
  echo "BamPath for Tag is not unique or exist. Please confirm!"
  exit
fi

CH_List="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

################
# GATK 4.2.0.0 #
################

GATKver="gatk4200"
#GnomadVCF=/rsrch3/home/leuk-rsrch/huryu/Database/hg19/GnomadDataForMutect/af-only-gnomad.raw.sites.h19.vcf
PON_PANEL_DIR="/rsrch3/scratch/leuk-rsrch/huryu/Lab_bamFile/Mutect_Line/PON_BuiltWithoutgNomad/MutectPONCreate/PON_PANEL/"
StoreDir=Mu2Call_ToFMC_NoGNOMAD_gatk4200/${Tag}/

mkdir -p ${StoreDir}

for CH in ${CH_List}
do
  CHR=chr${CH}
  PON_VCF=${PON_PANEL_DIR}/"pon_${CHR}.vcf.gz"
  mkdir -p ${StoreDir}/${CHR}/

  bsub -W 24:00 -q medium -n 2 -M 40G -R rusage[mem=40] -e ${StoreDir}/MutectCall_err.log -o ${StoreDir}/MutectCall_out.log \
       -J "MutectCall_${Tag}_${CH}_${TimeID}" \
       "${GATKver} Mutect2 -R ${RefFa} -L ${CH} -I ${BamPath} -pon ${PON_VCF} \
                     --f1r2-tar-gz ${StoreDir}/${CHR}/${Tag}.${CHR}.f1r2.tar.gz \
                     -O ${StoreDir}/${CHR}/${Tag}.${CHR}.unfiltered.vcf"
done

all_UnFilvcf=$(for CH in ${CH_List}; do   CHR=chr${CH};   printf -- "-I ${StoreDir}/${CHR}/${Tag}.${CHR}.unfiltered.vcf "; done)
all_UnFilvcfStats=$(for CH in ${CH_List}; do   CHR=chr${CH};   printf -- "-stats ${StoreDir}/${CHR}/${Tag}.${CHR}.unfiltered.vcf.stats "; done)
all_f1r2_input=$(for CH in ${CH_List}; do   CHR=chr${CH};   printf -- "-I ${StoreDir}/${CHR}/${Tag}.${CHR}.f1r2.tar.gz "; done)

MutectCallDoneSignal="done(MutectCall_${Tag}_1_${TimeID})"
for CH in ${CH_List}
do
  if [ ${CH} != "1" ]
  then
    MutectCallDoneSignal=${MutectCallDoneSignal}" && done(MutectCall_${Tag}_${CH}_${TimeID})"
  fi
done

bsub -W 24:00 -q medium -n 3 -M 20G -R rusage[mem=20] -e ${StoreDir}/MergeVCF_err.log -o ${StoreDir}/MergeVCF_out.log \
     -w "${MutectCallDoneSignal}" \
     -J "MergeVCF_${Tag}_${TimeID}" \
     "${GATKver} MergeVcfs ${all_UnFilvcf} \
                         -O ${StoreDir}/${Tag}.merged.unfiltered.vcf"

bsub -W 24:00 -q medium -n 3 -M 20G -R rusage[mem=20] -e ${StoreDir}/MMS_err.log -o ${StoreDir}/MMS_out.log \
     -w "${MutectCallDoneSignal}" \
     -J "MMS_${Tag}_${TimeID}" \
     "${GATKver} MergeMutectStats ${all_UnFilvcfStats} \
                                -O ${StoreDir}/${Tag}.merged.stats"

bsub -W 24:00 -q medium -n 3 -M 20G -R rusage[mem=20] -e ${StoreDir}/LROM_err.log -o ${StoreDir}/LROM_out.log \
     -w "${MutectCallDoneSignal}" \
     -J "LROM_${Tag}_${TimeID}" \
     "${GATKver} LearnReadOrientationModel ${all_f1r2_input} \
                                         -O ${StoreDir}/${Tag}.read-orientation-model.tar.gz"

#bsub -W 24:00 -q medium -n 4 -M 120G -R rusage[mem=120] -e ${StoreDir}/${Tag}_GPS_CC_err.log -o ${StoreDir}/${Tag}_GPS_CC_out.log \
#     -J "GPS_CC_${Tag}_${TimeID}" \
#     "${GATKver} ${JavaOption} GetPileupSummaries -I ${BamPath} -V ${GnomadVCF} -L ${GnomadVCF} -O ${StoreDir}/getpileupsummaries.table && \
#      ${GATKver} ${JavaOption} CalculateContamination -I ${StoreDir}/getpileupsummaries.table -tumor-segmentation ${StoreDir}/segments.table -O ${StoreDir}/calculatecontamination.table"

FMC_CMD="${GATKver} ${JavaOption} FilterMutectCalls -R ${RefFa} \
                                  -V ${StoreDir}/${Tag}.merged.unfiltered.vcf \
                                  --ob-priors ${StoreDir}/${Tag}.read-orientation-model.tar.gz \
                                  --stats ${StoreDir}/${Tag}.merged.stats \
                                  -O ${StoreDir}/${Tag}.merged.filtered.vcf"

bsub -W 24:00 -q medium -n 2 -M 40G -R rusage[mem=40] -e ${StoreDir}/${Tag}_FMC_err.log -o ${StoreDir}/${Tag}_FMC_out.log \
     -w "done(MergeVCF_${Tag}_${TimeID}) && done(MMS_${Tag}_${TimeID}) && done(LROM_${Tag}_${TimeID})" \
     -J "FMC_${Tag}_${TimeID}" \
     "${FMC_CMD}"
#     "${GATKver} ${JavaOption} FilterMutectCalls -R ${RefFa} \
#                                               -V ${StoreDir}/${Tag}.merged.unfiltered.vcf \
#                                               --tumor-segmentation  ${StoreDir}/segments.table \
#                                               --contamination-table ${StoreDir}/calculatecontamination.table \
#                                               --ob-priors ${StoreDir}/${Tag}.read-orientation-model.tar.gz \
#                                               --stats ${StoreDir}/${Tag}.merged.stats \
#                                               -O ${StoreDir}/${Tag}.merged.filtered.vcf"

exit
