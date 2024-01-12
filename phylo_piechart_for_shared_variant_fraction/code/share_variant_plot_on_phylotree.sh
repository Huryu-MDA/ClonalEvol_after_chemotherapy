#!/bin/bash

export LD_LIBRARY_PATH=/rsrch3/home/leuk-rsrch/huryu/miniconda3/lib:${LD_LIBRARY_PATH}

prefix=$(timeID)
FLG_LOH="FALSE"
while getopts c:p:P:v:t:q: OPT
do
  case $OPT in
    "c" ) FLG_LOH="TRUE" ; cnv_file_path="$OPTARG" ;;
    "p" ) prefix="$OPTARG" ;;
    "P" ) input_phylo="$OPTARG" ;;
    "t" ) input_nwk="$OPTARG" ;;
    "v" ) tMN_vcf_path="$OPTARG" ;;
    "q" ) queryMergeSamples="$OPTARG" ;; 
  esac
done
shift "$((OPTIND-1))"

tMN_vcf_base=$(basename ${tMN_vcf_path} .vcf)
tMN_snv_vcf=${tMN_vcf_base}.snvs_only.vcf

# argv[1] and argv[2] are used as input. argv[3] is the name of output file.
echo "tree_and_TipPhylo_To_PhydatWithNode.r ${input_nwk} ${input_phylo} out_phydat_${prefix}.fasta"
tree_and_TipPhylo_To_PhydatWithNode.r ${input_nwk} ${input_phylo} out_phydat_${prefix}.fasta

# Input file: out_phydat_${prefix}.fasta, Output: phydat_reshape_matrix_${prefix}.tsv
echo "phydat_reshape_to_mtx.sh out_phydat_${prefix}.fasta ${prefix}"
phydat_reshape_to_mtx.sh out_phydat_${prefix}.fasta ${prefix}

# convert tMN_unfiltered_merged_vcf to snv_based_vcf ---> output is snv_limited_vcf (${tMN_snv_vcf})
echo "snvs_vcf_convert.sh ${tMN_vcf_path}"
snvs_vcf_convert.sh ${tMN_vcf_path}

# Output: tMN_seq_${prefix}.txt (dummy_fasta with lower character)
echo "tMN_vcf_to_seq.py ${queryMergeSamples} ${tMN_snv_vcf} ${prefix}"
tMN_vcf_to_seq.py ${queryMergeSamples} ${tMN_snv_vcf} ${prefix}


# The process to create share_variant_fraction_matrix
if [ ${FLG_LOH} = "TRUE" ] 
then
  # Output: shared_var_frac_in_order_${prefix}.tsv ---> This part will ignore the unshared variant in nodes/tips of colonies that are located in the LOH region of tMN (bulk).
  phydat_matrix_to_shared_var_frac_with_CNV.py tMN_seq_${prefix}.txt phydat_reshape_matrix_${prefix}.tsv ${prefix} ${cnv_file_path} ${queryMergeSamples}
else
  # Output: shared_var_frac_in_order_${prefix}.tsv
  echo "phydat_matrix_to_shared_var_frac.py tMN_seq_${prefix}.txt phydat_reshape_matrix_${prefix}.tsv ${prefix}"
  phydat_matrix_to_shared_var_frac.py tMN_seq_${prefix}.txt phydat_reshape_matrix_${prefix}.tsv ${prefix}
fi

# plot_tree with pie_chart (or thermometer)
echo "_ans_seq_recon_.r ${input_nwk} shared_var_frac_in_order_${prefix}.tsv ${prefix}"
_ans_seq_recon_.r ${input_nwk} shared_var_frac_in_order_${prefix}.tsv ${prefix}

exit 0
