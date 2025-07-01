# Mu2_2ndCall — Forced Calling and Tag Matrix Generation

This directory includes the scripts used in the **second stage** of somatic mutation detection. It performs forced Mutect2 calls using merged VCF loci and generates per-sample tag matrices suitable for downstream filtering and visualization.

---

## Directory Structure

```
Mu2_2ndCall/
├── core_scripts/                # High-level orchestration for each stage
│   ├── VCF_merge_from_1stMu2call.sh
│   ├── 2ndMu2Call_on_MergedVCF.sh
│   ├── VCF_concat_merge_to_FilMTX.sh
│   ├── ForceCallonMergedVCF_SepChr_MedQue.sh
│   ├── VCF_SampleMerger.sh
│   └── VCF_SampleMerger_toChrSep.sh
├── filter_engine/              # Scripts to apply per-chromosome variant filters
│   ├── FilTagMatrixOnMergeVCFs.sh
│   ├── FilTagMatrixOnMergeVCFs_onChr.sh
│   └── __filMTXonChrs_merge__.sh
├── subcode_for_filtering/      # Helper Python scripts
│   ├── GermlF_map.py
│   └── NullReadtoZeroReplace.py
└── README.md                   # This file
```

---

## Overview of Workflow

1. **Merge per-sample VCFs** from `Mu2_1stCall`
2. **Split merged VCF by chromosome**
3. **Run Mutect2 with --alleles** to force-call variants at specified positions
4. **Merge resulting VCFs** per sample
5. **Generate a matrix of variants** and apply filtering

---

## Typical Execution

```bash
VCF_merge_from_1stMu2call.sh BamList TagList_on_PID02
2ndMu2Call_on_MergedVCF.sh BamList TagList_on_PID02 hg19
VCF_concat_merge_to_FilMTX.sh TagList_on_PID02 hg19 F
```

---

## Script Roles

### 🔧 `core_scripts/`
- Orchestrates the key phases:
  - Merging VCFs
  - Running Mutect2 forced calling
  - Concatenating VCFs
  - Launching matrix filtering

### 🧪 `filter_engine/`
- Applies variant-type and read-depth filters
- Produces final filtered matrix for downstream analysis

### 🧬 `subcode_for_filtering/`
- Python helpers used internally:
  - `GermlF_map.py`: Germline filtering
  - `NullReadtoZeroReplace.py`: Handles missing depth entries

---

## Output

Final output matrix:
```
MTX_for_VCF_filter.tsv
```

Other outputs:
- Per-chromosome filtered matrices
- Intermediate unfiltered/filtered VCFs
- Merge logs and job status summaries

---

## Dependencies

- **GATK** 4.2+ (Mutect2, IndexFeatureFile)
- **bcftools**
- **Python ≥ 3.6**
- LSF-based job system with:
  - `BsubS`, `BsubM`, `BsubL`
  - `WaitSignalMaker.sh`
  - `File_pickByTag.sh`

> These will be managed under a `LSF_utils/` module (to be added).

---

## Notes

- All scripts assume an LSF job scheduler (`bsub`)
- Output directories are automatically created under `FilMtx_From2ndMu2Call/`
- You can adapt scripts for SLURM/local use by replacing job submission logic

---

