#!/usr/bin/env python
# coding: utf-8

import sys, os
import numpy as np
import pandas as pd


#targetVCF = sys.argv[1]
targetVCF = sys.stdin

#filter_matrix = sys.argv[2]
filter_matrix = sys.argv[1]

#SaveFile = "germlF_LowVafSampleNumF_LowVAFmeanF_" + targetVCF

FilMtx = pd.read_csv(filter_matrix, sep = "\t")

# Make Tag column to concat VCF df

Chr_ser = FilMtx["#CHROM"].astype(str).copy()
POS_ser = FilMtx["POS"].astype(str).copy()
REF_ser = FilMtx["REF"].copy()
ALT_ser = FilMtx["ALT"].copy()

ConcatTag_ser = Chr_ser.copy()
ConcatTag_ser = ConcatTag_ser.str.cat(POS_ser, sep = "_")
ConcatTag_ser = ConcatTag_ser.str.cat(REF_ser, sep = "_")
ConcatTag_ser = ConcatTag_ser.str.cat(ALT_ser, sep = "_")

FilMtx["ConcatTag"] = ConcatTag_ser

#VCF_read = open(targetVCF)
VCF_read = targetVCF

# Header part is written into SaveFile
for line in VCF_read:
    if line[:2] == "##":
        #print(line.split("\n")[0], file = open(SaveFile, "a"))
        print(line.split("\n")[0])
        #pass
    elif line[:1] == "#":
        ColumnHead = line.split()
        #print(line.split("\n")[0], file = open(SaveFile, "a"))
        #print(ColumnHead)
        break

# Body part is converted into dataframe
VCF_df = pd.DataFrame([line.split() for line in VCF_read], columns=ColumnHead)

# Make Tag column to concat FilMtx df
Chr_ser_VCF = VCF_df["#CHROM"].astype(str).copy()
POS_ser_VCF = VCF_df["POS"].astype(str).copy()
REF_ser_VCF = VCF_df["REF"].copy()
ALT_ser_VCF = VCF_df["ALT"].copy()

ConcatTag_ser_VCF = Chr_ser_VCF.copy()
ConcatTag_ser_VCF = ConcatTag_ser_VCF.str.cat(POS_ser_VCF, sep = "_")
ConcatTag_ser_VCF = ConcatTag_ser_VCF.str.cat(REF_ser_VCF, sep = "_")
ConcatTag_ser_VCF = ConcatTag_ser_VCF.str.cat(ALT_ser_VCF, sep = "_")

VCF_df["ConcatTag"] = ConcatTag_ser_VCF

# FilTagNum + ConcatTag = 5. If the numver of filter in matrix increase, change here!

tag_df = FilMtx[FilMtx.columns[-5:]]

# merge VCF_body and tag using concat_tag
# tagged_VCF_df = pd.merge(left = VCF_df, right = tag_df, on = "ConcatTag", how = "outer")
tagged_VCF_df = pd.merge(left = VCF_df, right = tag_df, on = "ConcatTag", how = "inner")

del tagged_VCF_df["ConcatTag"]
tagged_VCF_df.fillna("None", inplace = True)


#tagged_VCF_df.to_csv(SaveFile,sep ="\t", index = False, mode = "a")
tagged_VCF_df.to_csv(sys.stdout,sep ="\t", index = False)

