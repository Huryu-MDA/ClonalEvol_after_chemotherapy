#!/usr/bin/env python
# coding: utf-8

# This is a code to define germline from the stdin from bcftools query info for AD?
# this can be recived only biallelic tsv input. 

import sys, os
import numpy as np
from scipy.stats import binom_test

# stdin = sys.stdin
# line = next(stdin)
tsvfile = sys.argv[1]
Read=open(tsvfile)

# [depth list in bcftools matrix (on column 4-)] into REF_int, ALT_int, Total_int depth list.
def depth_list_conversion(depth_list_str):
    ref_cov_arr = np.array([int(elm.split(",")[0]) for elm in depth_list_str])
    alt_cov_arr = np.array([int(elm.split(",")[1]) for elm in depth_list_str])
    total_cov_arr = ref_cov_arr + alt_cov_arr
    return ref_cov_arr, alt_cov_arr, total_cov_arr

# EUD: enough depth on total depth
def cov_arr_EUD(ref_cov_arr, alt_cov_arr, cutoff):
    total_cov_arr = ref_cov_arr + alt_cov_arr
    eu_depth_ind = np.where(total_cov_arr >= cutoff)[0]
    ref_cov_arr_EUD = ref_cov_arr[eu_depth_ind]
    alt_cov_arr_EUD = alt_cov_arr[eu_depth_ind]
    total_cov_arr_EUD = ref_cov_arr_EUD + alt_cov_arr_EUD
    return ref_cov_arr_EUD, alt_cov_arr_EUD, total_cov_arr_EUD


# for line in Read:
def GermlF_output(line):
    Line = line.split()
    CH, POS, REF, ALT  = Line[0:4]
    POS = int(POS)
    
    depth_list_str = Line[4:]
    ref_cov_arr, alt_cov_arr, total_cov_arr = depth_list_conversion(depth_list_str)
    ref_cov_sum = sum(ref_cov_arr)
    alt_cov_sum = sum(alt_cov_arr)
    total_cov_sum = ref_cov_sum + alt_cov_sum
    
    # alternate allele coverage for total depth >= 5 ( as cut off for 5)
    total_depth_cutoff = min(max(total_cov_arr) ,5)
    # EUD: Index of enough depth of total depth (REF_cov + ALT_cov)
    ref_cov_arr_EUD, alt_cov_arr_EUD, total_cov_arr_EUD = cov_arr_EUD(ref_cov_arr, alt_cov_arr, total_depth_cutoff)
    ref_cov_sum_EUD, alt_cov_sum_EUD, total_cov_sum_EUD = np.sum([ref_cov_arr_EUD, alt_cov_arr_EUD, total_cov_arr_EUD], axis = 1)
    #print(ref_cov_sum_EUD, alt_cov_sum_EUD, total_cov_sum_EUD)
    #break
    
    if  min(alt_cov_arr_EUD) >= 1:
        print(str(CH) + "\t" + str(POS) + "\t" + REF + "\t" + ALT + "\t" + str(ref_cov_sum) + "|" + str(alt_cov_sum) + "|" + str(total_cov_sum) + "\t" + "Min Alt coverage on enough depth" + " " + str(min(alt_cov_arr_EUD)))
    elif alt_cov_sum >= ref_cov_sum:
        print(str(CH) + "\t" + str(POS) + "\t" + REF + "\t" + ALT + "\t" + str(ref_cov_sum) + "|" + str(alt_cov_sum) + "|" + str(total_cov_sum) + "\t" + "RefCov <= AltCov")

    #BinomVal = binom_test(AltCov, TotalCov, p=0.5, alternative="less")
    # Null hypothesis: Alt is heterozygote germline 
    # Alternative hypothesis: Alt is somatic / hymozygoto germline (defined with p_value under 10*(-10))
    
    #BinomVal = binom_test(AltCovSum, TotalCovSum, p=0.5)
    #if BinomVal > 10 ** (-10):
    elif binom_test(alt_cov_sum_EUD, total_cov_sum_EUD, p=0.5) > 10 ** (-10):
        BinomVal = binom_test(alt_cov_sum_EUD, total_cov_sum_EUD, p=0.5)
        #heterozygote germline or homozygoto germline
        print(str(CH) + "\t" + str(POS) + "\t" + REF + "\t" + ALT + "\t" + str(ref_cov_sum_EUD) + "|" + str(alt_cov_sum_EUD) + "|" + str(total_cov_sum_EUD) + "\t" + "Binom of Alt to TotalCov on enough depth" + " " + str(BinomVal))
    else:
        # Somatic mutation candidate
        print(str(CH) + "\t" + str(POS) + "\t" + REF + "\t" + ALT + "\t" + str(ref_cov_sum) + "|" + str(alt_cov_sum) + "|" + str(total_cov_sum) + "\t" + "Somatic mutation candidate")
        #continue

GermlF_output_Run = [GermlF_output(line) for line in Read]

