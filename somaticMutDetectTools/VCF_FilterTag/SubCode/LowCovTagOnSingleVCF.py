#!/usr/bin/env python

import sys, argparse

parser = argparse.ArgumentParser()

parser.add_argument('-g', '--gender') 
parser.add_argument('-m', '--min_cov', default="6", help="ALT whose read is more than this value can be considered as putative somatic mutation.")

args = parser.parse_args()

# Gender will be "M or F"
Gender = args.gender
Min_Cov = int(args.min_cov)

Min_Cov_onSexChr = 3
if Min_Cov <= 4:
    Min_Cov_onSexChr = 2



FilTagName = "LowCovOnThisVCF"

Read = sys.stdin

for line in Read:
    if line[:2] == "##":
        print(line, end = "")
        continue
    elif line[:1] == "#":
        VCF_body_header = line.split("\n")[0] + "\t" + FilTagName + "\n"
        print(VCF_body_header, end = "")
        break

for line in Read:
    Line = line.split()
    CH = Line[0]
    AltCov = int(Line[9].split(":")[1].split(",")[1])

    if Gender == "M" and ((CH == "X") or (CH == "Y")):
        if AltCov < Min_Cov_onSexChr:
            LowCovTag = "LowCov_VCF"
        else:
            LowCovTag = "None"
    else:
        if AltCov < Min_Cov:
            LowCovTag = "LowCov_VCF"
        else:
            LowCovTag = "None"

    Line.append(LowCovTag)
    print("\t".join(Line))
