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
#for line in Read:
#    if line[:2] == "##":
#        #print(line.split("\n")[0], file = open(SaveFile, "a"))
#        print(line.split("\n")[0])
#    if line[:2] != "##":
#        ColumnNameOfBody = line.split()
#        break

def func_vartype(line):
    Line = line.split()
    REF = Line[3]
    ALT = Line[4]
    REF_len = len(REF)
    ALT_len = len(ALT)
    if (REF_len == 1) and (ALT_len == 1):
        return "SNV"
    elif (REF_len != 1) and (REF_len == ALT_len):
        return "MNV"
    else:
        return "INDEL"

#ColumnNameFile=sys.argv[1]
#ColumnNameOfBody=next(open(ColumnNameFile)).split("\n")[0].split()

#for line in Read:
def VarType_AddToLine(line):
    Line = line.split()
    CHR = Line[0]
    POS = Line[1]
    ID  = Line[2]
    REF = Line[3]
    ALT = Line[4]
    QUAL = Line[5]
    FILTER = Line[6]
    INFO = Line[7]
    VarType = func_vartype(line)
    Line.append(VarType)
    output_line = "\t".join(Line)
    print(output_line)

list(map(VarType_AddToLine, Read))
