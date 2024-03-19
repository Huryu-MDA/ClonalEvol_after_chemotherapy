#!/usr/bin/env python
# coding: utf-8


# This code is for the biallelic VCF (not for multiallelic)
import sys, os, glob, gzip
import numpy as np
import pandas as pd

# TargetVCF = "ChrToMu2F_INFO_VarTypeTagged__20211020225104.vcf"
# TargetVCF = sys.argv[1]
# Read = open(TargetVCF, "r")
Read = sys.stdin

Ref_Mtx = sys.argv[1]
# Ref_Mtx = "VariantDepthMatrix_Germl_Filter__20211020225104.tsv"
Ref_Mtx_Read = open(Ref_Mtx, "r")

# SaveFile = Filter + "FilterTagged_" + os.path.basename(TargetVCF)
# Filter = "LowMAP"
# TargetVCF = "GermBinomRemPlusHeader_9-PLATE1B5.vcf"
# LowMAPQ_Filter = "lowVAF_matrix.txt"


# Cut header part

# Print VCF_Header to save file and define column header of VCF_body
for line in Read:
    if line[:2] == "##":
        #print(line.split("\n")[0], file = open(SaveFile, "a"))
        print(line.split("\n")[0])
    if line[:2] != "##":
        ColumnNameOfBody = line.split()
        #Print out this part just before the print out of df part.
        break


filterTag_List = [line.split()[0] + "_" + line.split()[1] + "_" + line.split()[2] + "_" + line.split()[3] 
                  for line in Ref_Mtx_Read]
filterTag_Ser = pd.Series(filterTag_List)


VCFbody_List = [line.split() for line in Read]
df = pd.DataFrame(VCFbody_List, columns = ColumnNameOfBody)

# This part was changed from "df.iloc" to "df.loc".
# CHR_POS_REF_ALT_as_dfTab = [str(df.iloc[ind, 0]) + "_" + 
#                             str(df.iloc[ind, 1]) + "_" + 
#                             str(df.iloc[ind, 3]) + "_" + 
#                             str(df.iloc[ind, 4]) for ind in df.index]

CHR_POS_REF_ALT_as_dfTab = [str(df.loc[ind, "#CHROM"]) + "_" + 
                            str(df.loc[ind, "POS"]) + "_" + 
                            str(df.loc[ind, "REF"]) + "_" + 
                            str(df.loc[ind, "ALT"]) for ind in df.index]


filterName = "GermL_filter"
Filter = "GermLineBinom"
df["filter_Tag"] = pd.Series(CHR_POS_REF_ALT_as_dfTab)
df[filterName] = "None"
df[filterName][df["filter_Tag"].isin(filterTag_Ser)] = Filter
del df["filter_Tag"]


#df.head()

df.to_csv(sys.stdout,index=False, sep="\t")
