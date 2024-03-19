#!/usr/bin/env python
# coding: utf-8

import sys, os, glob, functools
import numpy as np
import pandas as pd
#from concurrent.futures import ProcessPoolExecutor

#targetVCF = sys.argv[1]
#targetVCF = "VCF_body_ChrToMu2F_INFO__20211020225104.vcf"

# In case to set the SaveFile, not stdout, arrange this part.
#SaveFile = "Processed_" + os.path.basename(targetVCF)

#Read = open(targetVCF, "r")
Read = sys.stdin

# Print VCF_Header to save file and define column header of VCF_body
for line in Read:
    if line[:2] == "##":
        #print(line.split("\n")[0], file = open(SaveFile, "a"))
        print(line.split("\n")[0])
    if line[:2] != "##":
        ColumnNameOfBody = line.split()
        break


#ColumnNameOfBody

VCFbody_List = [line.split() for line in Read]
df = pd.DataFrame(VCFbody_List, columns = ColumnNameOfBody)

df["POS"] = df["POS"].astype("int")
df["REF_len"] = pd.Series([len(elm) for elm in df["REF"]])
df["ALT_len"] = pd.Series([len(elm) for elm in df["ALT"]])
df["VarType"] = None

df.loc[(df["REF_len"] == 1) & (df["ALT_len"] == 1), "VarType"] = "SNV"
df.loc[(df["REF_len"] != 1) & (df["REF_len"] == df["ALT_len"]), "VarType"] = "MNV"
df["VarType"].fillna("INDEL", inplace = True)

# These codes below are commented off because of the SettingWithCopyWarning
# Also see "https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy"
#df["type"][(df["REF_len"] == 1) & (df["ALT_len"] == 1) ] = "SNV"
#df["type"][(df["REF_len"] != 1) & (df["REF_len"] == df["ALT_len"]) ] = "MNV"
#df["type"].fillna("INDEL", inplace = True)

# column pruning
del df["REF_len"]
del df["ALT_len"]

# This will be adjusted when using SaveFile as the data storing.
#df.to_csv(SaveFile,index=False, sep="\t")

df.to_csv(sys.stdout,index=False, sep="\t")
