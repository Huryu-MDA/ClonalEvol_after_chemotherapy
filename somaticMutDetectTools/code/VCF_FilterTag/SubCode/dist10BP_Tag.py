#!/usr/bin/env python
# coding: utf-8

# This code is origitated from FilterTagFor_SNV_INDEL_within10BP.py

import sys, os, glob, functools
import numpy as np
import pandas as pd

# targetVCF = "test100K_mono.vcf"
# Read = open(targetVCF, "r")
Read = sys.stdin


# In case of using SaveFile as data storing, arrange this part.
# SaveFile = "Processed_" + os.path.basename(targetVCF)

FilterName = "DistWithin10BP"

# To use functools.partial, the argument to use with map function as iter should be positioned in the 1st place.
def Within10BP_tag(ind, df):
    CH = df.loc[ind, "#CHROM"]
    POS = df.loc[ind, "POS"]
    # The condition of dist10BP (if there are equal or more than one loci that was not itself.)
    cond_dist10BP =  (CH == df["#CHROM"]) & (abs(df["POS"] - POS) <= 10) & (df["POS"] != POS)
    
    Within10BP_df = df[cond_dist10BP]
    if len(Within10BP_df) >= 1:
        return "WithIn10BP"
    else:
        return "None"

# Print VCF_Header to save file and define column header of VCF_body
for line in Read:
    if line[:2] == "##":
        #print(line.split("\n")[0], file = open(SaveFile, "a")) -> arrange this part when using SaveFile as data storing.
        print(line.split("\n")[0])
    if line[:2] != "##":
        ColumnNameOfBody = line.split()
        break

VCFbody_List = [line.split() for line in Read]
df = pd.DataFrame(VCFbody_List, columns = ColumnNameOfBody)

df["POS"] = df["POS"].astype("int")

ResultList = list(map(functools.partial(Within10BP_tag, df = df), df.index))
df[FilterName] = pd.Series(ResultList)

# This will be adjusted when using SaveFile as the data storing.
#df.to_csv(SaveFile,index=False, sep="\t")
df.to_csv(sys.stdout,index=False, sep="\t")
