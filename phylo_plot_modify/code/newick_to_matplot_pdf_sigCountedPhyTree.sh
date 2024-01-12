#!/bin/bash

module -q load R/4.1.0

linewidth_horiz=15
linewidth_vert=3
figsize_horiz=35
figsize_vert=30
figfmt="svg"

while getopts h:v:H:V:f: OPT
do
  case $OPT in
    h) linewidth_horiz="$OPTARG" ;;
    v) linewidth_vert="$OPTARG" ;;
    H) figsize_horiz="$OPTARG" ;;
    V) figsize_vert="$OPTARG" ;;
    f) figfmt="$OPTARG" ;;
  esac
done

shift "$((OPTIND-1))"

PID_prefix=$1
nwk_path=$2
ML_vcf_dir_path=$3
ML_mmSig_rel_matrix_path=$4
# e.g. PID0007 -> PHYML/ML_mmSig/mmsig_relative_20230327140128.tsv

# echo $PID_prefix $nwk_path $ML_vcf_dir_path $ML_mmSig_rel_matrix_path $linewidth_horiz $linewidth_vert $figsize_horiz $figsize_vert $figfmt
# exit 0

echo -e "\n"
echo -e ${PID_prefix}"\n"

read_nwk_info.r ${nwk_path}

# this part process PHYML/ML/*.vcf to count the tip / node variant. ---> output file: node_x.tsv, the matrix of ["node label", "SNV counts"] 
# this output would be processed by  phy_tree_plot_on_VarCounts_mmSig.py to plot the signature-based tree.
tip_node_var_counter.sh ${ML_vcf_dir_path}

phy_tree_plot.py ${PID_prefix}
#phy_tree_plot_on_VarCounts_mmSig.py --fg_f ${figfmt} --lw_v ${linewidth_vert} --lw_h ${linewidth_horiz} --fg_v ${figsize_vert} --fg_h ${figsize_horiz} ${PID_prefix} ${ML_mmSig_rel_matrix_path} node_x.tsv
phy_tree_plot_on_VarCounts_mmSig.py --fg_f ${figfmt} --lw_v ${linewidth_vert} --lw_h ${linewidth_horiz} --fg_v ${figsize_vert} --fg_h ${figsize_horiz} ${PID_prefix} ${ML_mmSig_rel_matrix_path}

# code confirmation ---> then this will record the last code in the code.log file.
#echo "phy_tree_plot_on_VarCounts_mmSig.py --fg_f ${figfmt} --lw_v ${linewidth_vert} --lw_h ${linewidth_horiz} --fg_v ${figsize_vert} --fg_h ${figsize_horiz} ${PID_prefix} ${ML_mmSig_rel_matrix_path} node_x.tsv" | tee code.log
echo "phy_tree_plot_on_VarCounts_mmSig.py --fg_f ${figfmt} --lw_v ${linewidth_vert} --lw_h ${linewidth_horiz} --fg_v ${figsize_vert} --fg_h ${figsize_horiz} ${PID_prefix} ${ML_mmSig_rel_matrix_path}" | tee code.log

exit 0
