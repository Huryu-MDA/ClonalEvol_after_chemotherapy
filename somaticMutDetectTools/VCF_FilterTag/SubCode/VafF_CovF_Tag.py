#!/usr/bin/env python
# coding: utf-8

import sys, os, glob, argparse, functools
import numpy as np
import pandas as pd

from decimal import *


parser = argparse.ArgumentParser(description='This code will recieve the stdin and yield the matrix attached tag_filter of VAF_cutoff and COV_cutoff.')

parser.add_argument("Variant_Read_Mtx", help="Args for Variant Read coverage matrix made with bcftools")
parser.add_argument("-mc", "--MinCov")
parser.add_argument("-lc", "--LessCovAltSample")
parser.add_argument("-g", "--Gender")

args = parser.parse_args()

CovMtx = args.Variant_Read_Mtx
CovCutoff = int(args.MinCov)
SampleNumCutOfLowCov = int(args.LessCovAltSample)
Gender = args.Gender

# setting the cutoff of Alt_Coverage on SexChr.
if CovCutoff < 5:
    CovCutoff_forSexChr = 2
else:
    CovCutoff_forSexChr = 3


#print(CovCutoff, SampleNumCutOfLowCov, Gender)
#sys.exit()

#CovCutoff = 5
#SampleNumCutOfLowCov = 5

def AF_compute(ALT_cov,DP):
    if DP == 0:
        AF = 0
    else:
        AF_deci = Decimal(ALT_cov) / Decimal(DP)
        AF = float(round(AF_deci, 3))
    return AF

# list or array-like object can be used.
# ALT_cov_list: Coverage of ALT from mutect for all samples
# DP_list: All depth (quality-filtered) on the locus of samples

# compute AF for each sample and return if the ALT on samples met criteria of sequential error or not

def LowVafonSample(ALT_cov_list, DP_list):
    
    Filter_LowVaf_Alert = "LowVAFinTooManySamples"
    
    AF_List = [AF_compute(Alt,DP) for Alt, DP in zip(Alt_CoV_List, DP_List)]
    AF_nonZero = [af for af in AF_List if af != 0]
    
    SampleNum = len(AF_List)
    # Warning!!! AF_nonZero should be used in this part.
    SampleNum_withVafUnder10 = len([elm for elm in AF_nonZero if elm < 0.10])
    
    if SampleNum_withVafUnder10 > SampleNum / 10:
        return Filter_LowVaf_Alert
    else:
        return "None"

def LowVafMean(ALT_cov_list, DP_list):
    
    Filter_LowVafMean_Alert = "LowVAFMeanInMutSamples"
    
    AF_List = [AF_compute(Alt,DP) for Alt, DP in zip(ALT_cov_list, DP_List)]
    AF_nonZero = [af for af in AF_List if af != 0]
    AF_nonZero_median = np.median(AF_nonZero)    
    
    # Checking if the max AF >= 0.4 and the the number of samples whose ALT coverage (<= 2) is less than 2 
    if (max(AF_List) >= 0.35) and (len([elm for elm in Alt_CoV_List if (elm > 0) and (elm <= 2)]) <= 3):
        AF_List = AFList_Reviser(ALT_cov_list, DP_List)    
        AF_nonZero = [af for af in AF_List if af != 0]
        AF_nonZero_median = np.median(AF_nonZero)

    if AF_nonZero_median <= 0.3:
        return Filter_LowVafMean_Alert
    else:
        return "None"
    
# Trimmer of small number of artifact (only less than two samples which less than 2 artifact read)
def AFList_Reviser(ALT_cov_list, DP_List):
    # Part of small artifact removing
    Index_artifact = [ind for ind, cov in enumerate(ALT_cov_list) if (cov > 0) and (cov <= 2)]
    for ind in Index_artifact:
        ALT_cov_list[ind] = 0
    # remake AF_List
    revised_AF_List = [AF_compute(Alt,DP) for Alt, DP in zip(ALT_cov_list, DP_List)]
    return revised_AF_List


# Function for LowCov On Samples

def CovOnSamples_cutoff(ALT_cov_list, Gender, CH, CovCutoff, CovCutoff_forSexChr = 3):
    
    Non_Zero_ALT_cov_list = [cov for cov in ALT_cov_list if cov != 0]
    
    if Gender == "F":
        if (CH != "Y"):
            SampleNum_LowCov = len([cov for cov in Non_Zero_ALT_cov_list if cov < CovCutoff])
            if SampleNum_LowCov >= SampleNumCutOfLowCov:
                return "TooManySampleOfLowCov"
            else:
                return "None"
            
        if (CH == "Y"):
            return "Yread_on_Female"
               
    elif Gender == "M":
        if (CH != "Y") and (CH != "X"):
            SampleNum_LowCov = len([cov for cov in Non_Zero_ALT_cov_list if cov < CovCutoff])
            if SampleNum_LowCov >= SampleNumCutOfLowCov:
                return "TooManySampleOfLowCov"
            else:
                return "None"
            
        if (CH == "Y") or (CH == "X"):
            SampleNum_LowCov = len([cov for cov in Non_Zero_ALT_cov_list if cov < CovCutoff_forSexChr])
            if SampleNum_LowCov >= SampleNumCutOfLowCov:
                return "TooManySampleOfLowCov"
            else:
                return "None"

if (Gender != "M") and (Gender != "F"):
    raise ValueError('sys.argv[i] should be "F" or "M"!')
    sys.exit(1)

# The concept of VAF in Mutect2 is slightly different what we thought. They set the value with the hypothesis that what value of AF is if the ALT read exist in that locus on the sample.
# VafMtx = "VarAfMtx_BiAlMergeVCF__20211021202117.tsv"


# 2 input files are required for this code. #1: variant read matrix for compute AF list of Coverage list #2: FilTagged_Chr_Pos_Ref_Alt_Mtx (The former body part of VCF).
#CovMtx = "Passed_VarCovMtx_BiAlMergeVCF__20211021202117.tsv"
#CovMtx = sys.argv[1]
Read_Cov = open(CovMtx, "r")

#InputFile = "Passed_VCF_body_ChrToMu2F_INFO_VarType_GermLfil_dist10BPfil__20211021202117.vcf"
#Read_Std = open(InputFile, "r")
Read_Std = sys.stdin


for line in Read_Std:
    if line[:2] == "##":
        print(line)
    elif line[:1] == "#":
        ColumnHead = line.split()
        break
    print("File format has no header part of VCF.")
    sys.exit()

ColumnHead.append("LowVafSampleNum_filter")
ColumnHead.append("LowVafMean_filter")
ColumnHead.append("SamplesOfLowCov_filter")

# The code to remove the LowCov on each sample should be written as other code (final process of filtering VCF)
#print("\t".join(ColumnHead), file = open("TempCheck.txt", "a"))
print("\t".join(ColumnHead))

for line in Read_Cov:
    StdIn = next(Read_Std)
    StdInLine = StdIn.split()
    CovLine = line.split()
    
    # Check if the position of StdIn and the VAF_matrix is same!.
    if (StdInLine[0] != CovLine[0]) or (StdInLine[1] != CovLine[1]) or (StdInLine[3] != CovLine[2]) or (StdInLine[4] != CovLine[3]):
        print("The position of stdin and VAF matrix is misaligned. Check both of their format.")
        sys.exit()
    
    CH, POS, REF, ALT = CovLine[0:4]
    POS = int(POS)

    Ref_Alt_Cov_List = CovLine[4:]
    Ref_CoV_List = [int(elm.split(",")[0]) for elm in Ref_Alt_Cov_List]
    Alt_CoV_List = [int(elm.split(",")[1]) for elm in Ref_Alt_Cov_List]
    DP_List = [Ref + Alt for Ref,Alt in zip(Ref_CoV_List, Alt_CoV_List)]

    FilTag_LowVaf = LowVafonSample(Alt_CoV_List, DP_List)
    FilTag_LowVafMean = LowVafMean(Alt_CoV_List, DP_List)
    FilTag_LowCovSampleNum = CovOnSamples_cutoff(Alt_CoV_List, Gender, CH, CovCutoff, CovCutoff_forSexChr)
    
    #print(StdIn + "\t" + FilTag_LowVaf + "\t" + FilTag_LowVafMean + "\t" + FilTag_LowCovSampleNum)
    #Alt_CoV_List = [int(elm.split(",")[1]) for elm in Ref_Alt_Cov_List]
    
    print(StdIn.split("\n")[0] + "\t" + FilTag_LowVaf + "\t" + FilTag_LowVafMean + "\t" + FilTag_LowCovSampleNum)
    #print(" ".join(map(str,DP_List)), file = open("TempCheck.txt", "a"))
    #print(" ".join(map(str,Alt_CoV_List)), file = open("TempCheck.txt", "a"))
    #print("", file = open("TempCheck.txt", "a"))
