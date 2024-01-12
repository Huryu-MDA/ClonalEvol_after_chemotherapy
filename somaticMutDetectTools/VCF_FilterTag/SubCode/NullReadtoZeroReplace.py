#!/usr/bin/env python
# coding: utf-8

import sys, os

#file = "BcfQuery_Biallelic_SampleMerged_VCFs.vcf"
#Read = open(file)
Read = sys.stdin

for line in Read:
    line = line.replace(",.\t", ",0\t")
    line = line.replace(".\t", "0,0\t")
    print(line, end = "")
    #break
