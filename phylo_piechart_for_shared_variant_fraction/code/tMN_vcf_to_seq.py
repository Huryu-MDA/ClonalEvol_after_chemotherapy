#!/usr/bin/env python
# coding: utf-8

import sys, os, glob
import numpy as np
import pandas as pd

def vcf_file_to_df(vcf_path):
    Read = open(vcf_path)
    for line in Read:
        if line[:2] == "##":
            pass
            #print(line.split("\n")[0], file = open(SaveFile, "a"))
            #print(line.split("\n")[0])
        if line[:2] != "##":
            ColumnHead = line.split()
            #Print out this part just before the print out of df part.
            break

    VCFbody_List = [line.split() for line in Read]
    VCF_df = pd.DataFrame(VCFbody_List, columns = ColumnHead)
    return VCF_df

def dtype_set(df, col_name):
    if col_name.upper() in ["CH", "CHR", "#CHROM"]:
        df[col_name] = df[col_name].astype("str")
    elif col_name.upper() in ["POS", "POSITION"]:
        df[col_name] = df[col_name].astype("int")
    elif col_name.upper() in ["REF", "ALT"]:
        df[col_name] = df[col_name].astype("str")
    else:
        pass

sample_merge_base = sys.argv[1]
#sample_merge_base = "PID0002_2010_phyml_trsf_box/QueryMTX_SampleMerged_MTXs.tsv"

bulk_vcf = sys.argv[2]
#bulk_vcf = "bulk_vcf/PID0002-S11015246.merged.filtered.vcf.mu2_passed.snvs_only.vcf"
#bulk_vcf = "bulk_vcf/PID0002-S11015246.merged.filtered.snvs_only.vcf"

prefix = sys.argv[3]
#prefix = "XXXXXX"

sample_merge_base_df = pd.read_csv(sample_merge_base, usecols=[0, 1, 2, 3], sep = "\t", header=None, 
                                   names=["#CHROM", "POS", "REF", "ALT"])

for col in sample_merge_base_df.columns:
    dtype_set(sample_merge_base_df, col)

# Remove bugs from previous bugs in bcftools
sample_merge_base_df["REF"] = sample_merge_base_df["REF"].str.split(" ", expand = True)[0]
sample_merge_base_df["ALT"] = sample_merge_base_df["ALT"].str.split(" ", expand = True)[0]
sample_merge_base_df["#CHROM"] = sample_merge_base_df["#CHROM"].str.split(" ", expand = True)[0]

vcf_df_whole = vcf_file_to_df(bulk_vcf)
vcf_df = vcf_df_whole[["#CHROM", "POS", "REF", "ALT"]].copy()

for col in vcf_df.columns:
    dtype_set(vcf_df, col)

vcf_df["ref_or_alt"] = "ALT"
vcf_df["alt_on_tMN"] = vcf_df["ALT"]

merged_df = pd.merge(left = sample_merge_base_df, right = vcf_df, on = ["#CHROM", "POS", "REF", "ALT"], how = "left")
merged_df["ref_or_alt"] = merged_df["ref_or_alt"].fillna("REF")
merged_df["alt_on_tMN"] = merged_df["alt_on_tMN"].fillna(merged_df["REF"])

bulk_seq = "".join(list(merged_df["alt_on_tMN"]))
bulk_seq_lower = bulk_seq.lower()

print(bulk_seq_lower, file=open("tMN_seq_" + prefix + ".txt", "w"))

exit()
