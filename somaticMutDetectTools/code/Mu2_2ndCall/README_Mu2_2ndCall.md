# Mu2_2ndCall — Forced Calling and Tag Matrix Generation

This directory includes the scripts used in the **second stage** of somatic mutation detection. It performs forced Mutect2 calls using merged VCF loci and generates per-sample tag matrices suitable for downstream filtering and visualization.

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

## Script Breakdown

### 🧩 `VCF_merge_from_1stMu2call.sh`
- Entry point for merging `.merged.unfiltered.vcf` files from the first Mutect2 call.
- Prepares and calls `VCF_SampleMerger_toChrSep.sh`.

### 🧩 `VCF_SampleMerger_toChrSep.sh`
- Merges VCFs using `bcftools` and splits them into `Chr*.vcf` files.
- Indexes each chromosome VCF for later use in Mutect2.

### 🧩 `2ndMu2Call_on_MergedVCF.sh`
- Dispatches sample-wise jobs for Mutect2 forced calling on merged VCF loci.
- Internally uses `ForceCallonMergedVCF_SepChr_MedQue.sh`.

### 🧩 `ForceCallonMergedVCF_SepChr_MedQue.sh`
- Runs Mutect2 for each chromosome with `--alleles` and `-L` options.
- Produces `.Chr*.unfiltered.vcf` files and prepares sample-level merge commands.

### 🧩 `VCF_concat_merge_to_FilMTX.sh`
- Waits until all forced-called VCFs are complete.
- Merges per-sample `.ChrMerged.unfiltered.vcf` files using `VCF_SampleMerger.sh`.
- Initiates matrix generation using `FilTagMatrixOnMergeVCFs.sh`.

### 🧩 `FilTagMatrixOnMergeVCFs.sh`
- Splits merged VCF into chromosomes and applies filtering on each via:
  - `FilTagMatrixOnMergeVCFs_onChr.sh` (per-chromosome)
  - `__filMTXonChrs_merge__.sh` (combines into single TSV matrix)

---

## Output

Final matrix file:  
```
MTX_for_VCF_filter.tsv
```

Other outputs:
- Chromosome-separated intermediate VCFs
- VCF merge logs
- AD-based variant matrices

---

## Dependencies

- GATK 4.2+ (Mutect2, IndexFeatureFile)
- bcftools
- Python scripts (e.g., `NullReadtoZeroReplace.py`, `GermlF_map.py`, etc.)
- LSF job submission environment with `BsubS`, `BsubM`, `WaitSignalMaker.sh`, etc.

> These utilities will be grouped under `LSF_utils/` and added soon.

---

## Notes

- All steps rely on `bsub` job scheduling.
- Be sure to check memory (`-M`) and CPU thread settings (`-n`) for each stage.
- Output paths and intermediate directories are auto-generated under `FilMtx_From2ndMu2Call/`.

