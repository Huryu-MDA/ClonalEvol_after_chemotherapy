#!/usr/bin/env python
# coding: utf-8

import sys, os, glob
import numpy as np
import pandas as pd

def dtype_set(df, col_name):
    if col_name.upper() in ["CH", "CHR", "#CHROM"]:
        df[col_name] = df[col_name].astype("str")
    elif col_name.upper() in ["POS", "POSITION", "START", "END"]:
        df[col_name] = df[col_name].astype("int")
    elif col_name.upper() in ["REF", "ALT"]:
        df[col_name] = df[col_name].astype("str")
    else:
        pass

def somatic_mut_set(refseq_seq, alt_seq):
    if len(refseq_seq) != len(alt_seq):
        print("Use same length sequence as args.")
    
    somatic_list = [[ind, alt_seq[ind]] for ind in range(len(refseq_seq)) if (refseq_seq[ind] != alt_seq[ind])]
    #somatic_df = pd.DataFrame(somatic_list, columns = ["INDEX", "ALT"])
    somatic_set = set([str(elm[0]) + "_" + elm[1] for elm in somatic_list])
    #return somatic_set, somatic_df
    return somatic_set

# LOH_df should be processed with dtype_set before being used this function
def chrom_pos_in_LOH(chrom, pos, LOH_df):
    chrom = str(chrom).split(" ")[0]
    pos = int(pos)
    LOH_chrom_df = LOH_df[LOH_df["#CHROM"] == chrom].copy()
    overlap_LOH_region = LOH_chrom_df[(LOH_chrom_df["START"] <= pos) & (LOH_chrom_df["END"] >= pos)]
    
    if len(overlap_LOH_region) == 1:
        return True
    elif len(overlap_LOH_region) == 0:
        return False
    else:
        print("#CHROM:", chrom, ", POS:", pos, ", CNV_file:", cnv_region)
        raise Exception("The variant exist in more than two region in CNV_LOH_region.")

# id_mut example: "1034_t", "10462_t", "10463_g", ...
# query_df: colony_merged_df (merge with bcftools)
# LOH_df: LOH region that was derived from ASCAT
def id_mut_on_LOH_of_bulk(id_mut, query_df, LOH_df):
    mut_index = int(id_mut.split("_")[0])
    chrom_of_mut = query_df["#CHROM"][mut_index]
    pos_of_mut = query_df["POS"][mut_index]
    on_LOH = chrom_pos_in_LOH(chrom_of_mut, pos_of_mut, LOH_df)
    return on_LOH

tMN_seq_file = sys.argv[1]
#tMN_seq_file = "tMN_seq_PID0002-S11015246.txt"
#tMN_seq_file = "tMN_seq"
tMN_seq = next(open(tMN_seq_file, "r")).split("\n")[0]

tip_node_snv_mtx_input = sys.argv[2]
#tip_node_snv_mtx_input = "phydat_reshape_matrix_PID0002-S11015246.tsv"
#tip_node_snv_mtx_input = "phydat_reshape_matrix.tsv"
tip_node_snv_mtx_df = pd.read_csv(tip_node_snv_mtx_input, sep = "\t", names=["labels", "seq"])

prefix = sys.argv[3]
#prefix = "PID0002-S11015246"
output_df_label_sharedVarFrac = "shared_var_frac_in_order_" + prefix + ".tsv"

cnv_region = sys.argv[4]
#cnv_region = "PID0002-S11015246-50xWGS.copynumber.caveman.sorted.bed"
cnv_df = pd.read_csv(cnv_region, sep = "\t")
LOH_region_df = cnv_df[cnv_df["Hetero_Or_LOH"] == "LOH"][["#CHROM", "START", "END"]].copy()

queryMergeSamples = sys.argv[5]
#queryMergeSamples = "QueryMTX_SampleMerged_MTXs.tsv"
queryColony_merge_df = pd.read_csv(queryMergeSamples, sep = "\t", header = None, usecols=[0,1,2,3], names=["#CHROM", "POS", "REF", "ALT"])
#queryColony_merge_df.head()

for col in LOH_region_df.columns:
    dtype_set(LOH_region_df, col)

for col in queryColony_merge_df.columns:
    dtype_set(queryColony_merge_df, col)

# This line is to remove the extra space after 12, 21, 22, X, Y. This bug may be derived from bcftools (older version?) or pandas read function
queryColony_merge_df["#CHROM"] = queryColony_merge_df["#CHROM"].str.split(" ", expand = True)[0]

refseq_seq_ser = tip_node_snv_mtx_df[tip_node_snv_mtx_df["labels"] == "REFSEQ"]["seq"]
refseq_seq = refseq_seq_ser.iloc[0]

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
    unshared_var = somat_set_on_labels[label] - tMN_mut_set
    unshared_var_sorted = sorted(unshared_var)
    # unshared_var_of_colony/node on LOH_of_bulk should be removed from fraction counts because of its impossibility.
    unshared_var_on_LOH = [var for var in unshared_var_sorted 
                           if (id_mut_on_LOH_of_bulk(var, queryColony_merge_df, LOH_region_df) == True)]
    unshared_var_off_LOH = [var for var in unshared_var_sorted 
                           if (id_mut_on_LOH_of_bulk(var, queryColony_merge_df, LOH_region_df) == False)]
    
    if len(somat_set_on_labels[label]) == 0:
        shared_var_frac = 1.0
    else:
        #shared_var_frac = len(shared_var) / len(somat_set_on_labels[label])
        shared_var_frac = len(shared_var) / (len(shared_var) + len(unshared_var_off_LOH))
    #print(label, shared_var_frac)
    shared_var_frac_dict[label] = shared_var_frac
    print(len(shared_var), len(unshared_var_off_LOH), len(unshared_var_on_LOH), len(somat_set))
    #break

shared_var_df = pd.DataFrame([[key, shared_var_frac_dict[key]] for key in shared_var_frac_dict.keys()], columns=["labels", "shared_var_frac"])

label_seq_sharedVar_merged_df = pd.merge(left = tip_node_snv_mtx_df, right=shared_var_df, on="labels", how = "outer")

sharedVar_merged_df = label_seq_sharedVar_merged_df[["labels", "shared_var_frac"]].copy()
sharedVar_merged_df.to_csv(output_df_label_sharedVarFrac, sep = "\t", index=False, header = False)

exit()
