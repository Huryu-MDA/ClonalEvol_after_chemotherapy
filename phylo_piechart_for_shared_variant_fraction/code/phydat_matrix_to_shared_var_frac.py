#!/usr/bin/env python
# coding: utf-8

import sys, os, glob
import numpy as np
import pandas as pd

tMN_seq_file = sys.argv[1]
#tMN_seq_file = "tMN_seq"
tMN_seq = next(open(tMN_seq_file, "r")).split("\n")[0]

tip_node_snv_mtx_input = sys.argv[2]
#tip_node_snv_mtx_input = "phydat_reshape_matrix.tsv"
tip_node_snv_mtx_df = pd.read_csv(tip_node_snv_mtx_input, sep = "\t", names=["labels", "seq"])

prefix = sys.argv[3]
#prefix = "XXXXXX"
output_df_label_sharedVarFrac = "shared_var_frac_in_order_" + prefix + ".tsv"

refseq_seq_ser = tip_node_snv_mtx_df[tip_node_snv_mtx_df["labels"] == "REFSEQ"]["seq"]
refseq_seq = refseq_seq_ser.iloc[0]

def somatic_mut_set(refseq_seq, alt_seq):
    if len(refseq_seq) != len(alt_seq):
        print("Use same length sequence as args.")
    
    somatic_list = [[ind, alt_seq[ind]] for ind in range(len(refseq_seq)) if (refseq_seq[ind] != alt_seq[ind])]
    #somatic_df = pd.DataFrame(somatic_list, columns = ["INDEX", "ALT"])
    somatic_set = set([str(elm[0]) + "_" + elm[1] for elm in somatic_list])
    #return somatic_set, somatic_df
    return somatic_set


tMN_mut_set = somatic_mut_set(refseq_seq, tMN_seq)

somat_set_on_labels = {}
#somat_df_on_labels = {}
shared_var_frac_dict = {}

for ind in tip_node_snv_mtx_df.index:
    label = tip_node_snv_mtx_df.loc[ind, "labels"]
    alt_seq = tip_node_snv_mtx_df.loc[ind, "seq"]
    #somat_set, somat_df = somatic_mut_set_or_df(refseq_seq, alt_seq)
    somat_set = somatic_mut_set(refseq_seq, alt_seq)
    somat_set_on_labels[label] = somat_set
    #somat_df_on_labels[label] = somat_df
    
    shared_var = tMN_mut_set & somat_set_on_labels[label]
    if len(somat_set_on_labels[label]) == 0:
        shared_var_frac = 1.0
    else:
        shared_var_frac = len(shared_var) / len(somat_set_on_labels[label])
    #print(label, shared_var_frac)
    shared_var_frac_dict[label] = shared_var_frac
    #break

shared_var_df = pd.DataFrame([[key, shared_var_frac_dict[key]] for key in shared_var_frac_dict.keys()], columns=["labels", "shared_var_frac"])

label_seq_sharedVar_merged_df = pd.merge(left = tip_node_snv_mtx_df, right=shared_var_df, on="labels", how = "outer")

sharedVar_merged_df = label_seq_sharedVar_merged_df[["labels", "shared_var_frac"]].copy()
sharedVar_merged_df.to_csv(output_df_label_sharedVarFrac, sep = "\t", index=False, header = False)


exit()


"""
somat_set_on_labels["118"] - tMN_mut_set
len(somat_set_on_labels["198"])
len(somat_set_on_labels["197"])

len(somat_set_on_labels["198"]
len(somat_set_on_labels["198"] - tMN_mut_set)

# 197 will be the true node to connect
len(somat_set_on_labels["197"] - tMN_mut_set)

somat_set_on_labels["198"] - tMN_mut_set
somat_set_on_labels["197"] - tMN_mut_set

test_set = (somat_set_on_labels["197"] - tMN_mut_set) - (somat_set_on_labels["198"] - tMN_mut_set)


(somat_set_on_labels["198"] - tMN_mut_set) - (somat_set_on_labels["197"] - tMN_mut_set)
"""
exit()
