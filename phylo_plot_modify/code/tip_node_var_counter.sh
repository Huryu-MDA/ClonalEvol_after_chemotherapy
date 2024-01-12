#!/bin/bash

ML_vcf_dir_path=$1

# Part to decide the X_position from ML?*.vcf for each phyml
# The base tree in PID0007: PHYML/ROOTED_QueryMTX_SampleMerged_MTXs_phylo_seq58591_REFSEQaddedAsRoot.newick
# Should be arraged for the file path.
wc -l ${ML_vcf_dir_path}/*.vcf | rev | cut -f 2 -d " " | rev > node_x_X.txt
wc -l ${ML_vcf_dir_path}/*.vcf | rev | cut -f 1 -d " " | cut -f 1 -d "/" | cut -f 2 -d "." | rev > node_x_NodeLabel.txt
paste node_x_NodeLabel.txt node_x_X.txt | grep -v total > node_x.tsv
echo -e "REFSEQ""\t"0 >> node_x.tsv
rm node_x_X.txt node_x_NodeLabel.txt
