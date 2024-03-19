#!/bin/bash

TagList=$1

mkdir -p 7FilTagged_VCF_Dir/TagTrimmedVCF/
mkdir -p 7FilTagged_VCF_Dir/WholeLocus
mv 7FilTagged_*.vcf 7FilTagged_VCF_Dir/WholeLocus/.
mv 7Fil_passed_7FilTagged_*.vcf 7FilTagged_VCF_Dir/.
rm SNV_count_?Fil.tsv MNV_count_?Fil.tsv INDEL_count_?Fil.tsv SNV_SampleNameCol.tsv MNV_SampleNameCol.tsv INDEL_SampleNameCol.tsv JobID_*

cd 7FilTagged_VCF_Dir/TagTrimmedVCF/
ls ../*.vcf > VCFList

# This part may change the name of the TagList. This works well.
cp ../../${TagList} ./TagList
mkdir -p Annot_Annovar Phylo SigPro TaggedVCF TagTrimmed_SNV_VCF

for i in $(cat TagList)
do
  TargetVCF=$(cat VCFList | grep ${i})
  cat ${TargetVCF} | cut -f 1-10 > 7Fil_TagTrimmed.${i}.vcf
done

mv ../*.vcf TaggedVCF/.

for i in $(cat TagList)
do
  VCF=$(ls TaggedVCF/*vcf | grep "${i}")
  cat ${VCF} | grep ^# | cut -f 1-10 > TagTrimmed_SNV_VCF/${i}_TagTrimmed_SNV.vcf
  cat ${VCF} | grep "SNV" | cut -f 1-10 >> TagTrimmed_SNV_VCF/${i}_TagTrimmed_SNV.vcf
done

mkdir -p TagTrimmed_INDEL_VCF/
for i in $(cat TagList)
do
 VCF=$(ls TaggedVCF/*vcf | grep "${i}")
 cat ${VCF} | grep ^# | cut -f 1-10 > TagTrimmed_INDEL_VCF/${i}_TagTrimmed_INDEL.vcf
 cat ${VCF} | grep "INDEL" | cut -f 1-10 >> TagTrimmed_INDEL_VCF/${i}_TagTrimmed_INDEL.vcf
done

rm TagList VCFList
mv Annot_Annovar/ ../.
mv Phylo/ ../.
mv SigPro/ ../.
mv TaggedVCF/ ../.
mv TagTrimmed_INDEL_VCF/ ../.
mv TagTrimmed_SNV_VCF/ ../.
