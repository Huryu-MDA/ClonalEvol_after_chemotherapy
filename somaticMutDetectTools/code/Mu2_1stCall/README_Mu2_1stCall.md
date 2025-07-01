# Mu2_1stCall — Initial Somatic Mutation Calling

This directory contains bash scripts for the **first step of somatic mutation detection** using GATK Mutect2. Each script performs per-chromosome Mutect2 calling with Panel of Normals (PON) support.

## Overview

The scripts in this folder launch parallel Mutect2 jobs across chromosomes for each sample tag listed in a BAM list. They support flexible configurations, including variant calling with or without gnomAD filtering and with support for multiple GATK versions (e.g., 4.2.0.0 or 4.3.0.0).

## Available Scripts

- `Mu2Call1st_LROM_GPS_CC_Filter_without_gNomad.sh`  
  Default script used for most patient samples on our HPC system.  
- `Mu2Call1st_LROM_GPS_CC_Filter_without_gNomad_BulkSample.sh`  
  Variant tuned for bulk sequencing samples; includes higher memory allocations.  
- `Mu2Call1st_LROM_GPS_CC_Filter_without_gNomad_gatk4300.sh`  
  Version-specific script for GATK 4.3.0.0.

## Typical Usage

This is the most commonly used batch submission pattern:

```bash
cat TagList_409632_2000_2003 | grep -v "^#" | \
xargs -i BsubS -n 2 -M 4 "Mu2Call1st_LROM_GPS_CC_Filter_without_gNomad.sh BamList {}" Mu2Call_to_FMC \
>> JobID_Mu2CalltoFMC
```

Where:
- `TagList_409632_2000_2003`: Contains sample tags
- `BamList`: Maps each tag to a BAM file path
- `BsubS`: HPC wrapper for `bsub` job submission
- `Mu2Call_to_FMC`: Queue or label used in job logging

Each job will output unfiltered VCFs per chromosome and `.f1r2.tar.gz` files under:
```
Mu2Call_ToFMC_NoGNOMAD_gatk4200/SAMPLE_TAG/chr*/SAMPLE_TAG.chr*.unfiltered.vcf
```

## Notes

- All scripts assume LSF job scheduling (`bsub`)
- You can modify memory and thread settings depending on sample size or environment
- GATK version must match the script used (`gatk4200`, `gatk4300`, etc.)
