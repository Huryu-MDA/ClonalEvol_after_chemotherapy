#!/usr/bin/env python
# coding: utf-8

# This is a code to define germline from the stdin from bcftools query info for AD?

import sys, os
import numpy as np
from scipy.stats import binom_test

VarMtx=sys.argv[1]
# VarMtx = "Variant_Depth_Matrix.tsv"

Read = open(file = VarMtx)

for line in Read:
    Line = line.split()
    CH  = Line[0]
    POS = Line[1]
    Ref = Line[2]
    Alt = Line[3]
    
    DepthList_Str = Line[4:]
    
    #print(CH, POS, Ref, Alt)
    
    ReadList = list(Ref.split(","))
    ReadList.extend(Alt.split(","))
    DepthList = np.array([list(map(int, elm.split(","))) for elm in DepthList_Str]).sum(axis = 0)
    
    #print(DepthList)
    RefDep = DepthList[0]
    
    # Checkpoint if the Alt has alt2/alt3/ or ....
    #if len(AltDepthListDescending) >= 2:
    #    print(CH, POS, Ref, Alt)
    #    print(AltDepthListDescending)
    #    break
    #else:
    #    continue
    
    # This is when there are second allele, third allele....
    AltDepthListDescending = DepthList[1:].argsort()[::-1]
    MaxAltIndex = AltDepthListDescending[0]
    MaxAlt = ReadList[1:][MaxAltIndex]
    MaxAltDepth = DepthList[1:][MaxAltIndex]
    
    TotalDepth = DepthList.sum()
    
    # minimum Depth of ALT check. if this is not "0", the ALT should be considered as germline.
    ALT_depth = [int(elm.split(",")[1]) for elm in DepthList_Str]
    minALT_depth = min(ALT_depth)
    if minALT_depth >= 1:
        print(str(CH) + "\t" + str(POS) + "\t" + Ref + "\t" + MaxAlt + "\t" + str(RefDep) + "|" + str(MaxAltDepth) + "|" + str(TotalDepth) + "\t" + str(minALT_depth))
        continue
        # break
    
    # BinomialDistributionTest
    BinomVal = binom_test(MaxAltDepth, TotalDepth, p=0.5, alternative="less")
    
    if (BinomVal > 10 ** (-10)) or (MaxAltDepth >= RefDep):
        #Alt_GermLine_P = True
        #print(line, file = open("GermLineCand.txt", mode = "a"))
        #print([BinomVal, Ref, RefDep, MaxAlt, MaxAltDepth, TotalDepth], file = open("GermLineCand_Stat.txt", mode = "a"))
        print(str(CH) + "\t" + str(POS) + "\t" + Ref + "\t" + MaxAlt + "\t" + str(RefDep) + "|" + str(MaxAltDepth) + "|" + str(TotalDepth) + "\t" + str(BinomVal))
    
