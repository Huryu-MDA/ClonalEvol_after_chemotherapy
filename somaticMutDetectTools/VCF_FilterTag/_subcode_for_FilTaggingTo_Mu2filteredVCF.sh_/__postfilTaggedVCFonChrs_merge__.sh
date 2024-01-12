#!/bin/bash

CH_List="1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n21\n22\nX\nY"
prefix=$1

# VCF_header
cat chr1/7FilTagged_${prefix}_chr1.merged.filtered.vcf | grep "^#" > ../7FilTagged_${prefix}.merged.filtered.vcf

# VCF_body
for i in $(echo -e ${CH_List})
do
  cat chr${i}/7FilTagged_${prefix}_chr${i}.merged.filtered.vcf | grep -v "^#" >> ../7FilTagged_${prefix}.merged.filtered.vcf
done

exit 0
