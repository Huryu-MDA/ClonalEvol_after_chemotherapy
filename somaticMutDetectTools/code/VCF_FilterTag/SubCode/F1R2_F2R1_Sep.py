#!/usr/bin/env python
# coding: utf-8

import sys, os

file = os.path.realpath(sys.argv[1])

def Return_RefAlt(file):
    Ref_Save = os.path.dirname(file) + "/" + "Ref_" + os.path.basename(file)
    Alt_Save = os.path.dirname(file) + "/" + "Alt_" + os.path.basename(file)
    
    Read = open(file)
    
    for line in Read:
        Line = line.split()
        #print(Line)
        Sample_Ref = [elm.split(",")[0] for elm in Line[4:]]
        Sample_Alt = [elm.split(",")[1] for elm in Line[4:]]
        print("\t".join(Sample_Ref), file = open(Ref_Save, "a"))
        print("\t".join(Sample_Alt), file = open(Alt_Save, "a"))

Return_RefAlt(file)
