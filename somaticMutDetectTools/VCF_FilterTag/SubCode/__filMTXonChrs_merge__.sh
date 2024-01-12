#!/bin/bash

# child_process_of_FilTagMatrixOnMergeVCFs.sh
# post processing of FilTagMatrixOnMergeVCFs_onChr.sh

CH_List="1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n21\n22\nX\nY\nMT"
#echo -e ${CH_List}
#exit 0

cat chr1/chr1_MTX_for_VCF_filter.tsv | grep "^#" > MTX_for_VCF_filter.tsv

for i in $(echo -e ${CH_List})
do
  cat chr${i}/chr${i}_MTX_for_VCF_filter.tsv | grep -v "^#" >> MTX_for_VCF_filter.tsv
done
exit 0
