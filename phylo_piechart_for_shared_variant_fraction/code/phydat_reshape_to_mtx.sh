#!/bin/bash

input_phydat=$1
#prefix shoud be added!
prefix=$2

cat ${input_phydat} | grep "^>" | cut -f 2 -d ">" > tip_node_labels_${prefix}.txt
cat ${input_phydat} | grep -v "^>" > tip_node_fasta_${prefix}.txt
paste tip_node_labels_${prefix}.txt tip_node_fasta_${prefix}.txt -d "\t" > phydat_reshape_matrix_${prefix}.tsv

rm tip_node_labels_${prefix}.txt tip_node_fasta_${prefix}.txt
