module load bcftools

ref_fa="/rsrch3/home/leuk-rsrch/huryu/Database/hg19/RefFasta/Homo_sapiens_assembly19.fasta"
vcf_raw=$1
vcf_raw_base=$(basename ${vcf_raw} .vcf)

bcftools norm -m- ${vcf_raw} | bcftools norm -f ${ref_fa} > ${vcf_raw_base}.biallelic.vcf
# bcftools norm -m- ${vcf_raw} > ${vcf_raw_base}.biallelic.vcf

bcftools view -v snps ${vcf_raw_base}.biallelic.vcf -o ${vcf_raw_base}.snvs_only.vcf
# bcftools view -v snps PID0002-S11015246.merged.filtered.vcf.mu2_passed.vcf -o PID0002-S11015246.merged.filtered.vcf.mu2_passed.snvs_only.vcf

rm ${vcf_raw_base}.biallelic.vcf
