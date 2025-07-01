# Clonal Evolution of Hematopoietic Stem Cells After Cancer Chemotherapy

This repository contains the full pipeline and scripts used in the analysis presented in:

**Uryu, H. et al.**  
*Clonal evolution of hematopoietic stem cells after autologous stem cell transplantation.*  
*Nature Genetics* (2025) [https://doi.org/10.1038/s41588-025-02235-w](https://doi.org/10.1038/s41588-025-02235-w)

---

## Repository Overview

This project is divided into modular, checkpointed components that align with the somatic mutation detection and phylogenetic analysis process:

```
somaticMutDetectTools/
├── code/
│   ├── LSF_utils/         # LSF job wrappers and helper utilities
│   ├── Mu2_1stCall/       # Initial GATK Mutect2 per-chromosome calling
│   ├── Mu2_2ndCall/       # Forced calling and tag matrix construction
│   └── VCF_FilterTag/     # Final filtering and annotation of VCFs
├── test_data/             # Panel of Normals and TP53 BAMs for testing
├── test_run/              # Output from running test data
└── README.md              # This file
```

Additional analysis components:
```
phylo_plot_modify/                           # Convert tree branches to mutational signature profiles
phylo_piechart_for_shared_variant_fraction/  # Shared variant pie chart visualization
```

---

## Requirements and Setup

This project relies primarily on shell scripting and the GATK toolkit for somatic mutation detection.

### Software Dependencies

- **Bash** (tested on Linux environment)
- **Java** 11 (e.g., `openjdk/11.0.5`)
- **GATK** 4.2.0.0
- **LSF Job Scheduler** (uses `bsub` for job submission)

> Adapt `bsub` and `module load` commands if you're not using an LSF-based HPC system.

### Required Reference Files

| File                          | Description                                                  |
|-------------------------------|--------------------------------------------------------------|
| `Homo_sapiens_assembly19.fasta` | Reference genome (hg19)                                      |
| `pon_chr*.vcf.gz`               | Panel of Normals for each chromosome                         |
| `bamlist.txt`                   | File containing full paths to BAM files with associated tags |
| *(optional)* `gnomAD.vcf`       | For contamination filtering (currently commented out)        |

---

## How to Use

Each major step is designed to be run **independently and modularly**.

### 🧪 1. `Mu2_1stCall/`
Run GATK Mutect2 for initial per-chromosome variant calling.

```bash
bash Mu2_1stCall/Mu2Call1st_LROM_GPS_CC_Filter_without_gNomad.sh BamList SampleTag
```

### 🔁 2. `Mu2_2ndCall/`
Force call known variant loci and create a tag matrix for filtering.

```bash
bash Mu2_2ndCall/core_scripts/VCF_merge_from_1stMu2call.sh BamList TagList
bash Mu2_2ndCall/core_scripts/2ndMu2Call_on_MergedVCF.sh BamList TagList hg19
bash Mu2_2ndCall/core_scripts/VCF_concat_merge_to_FilMTX.sh TagList hg19 F
```

### 🧬 3. `VCF_FilterTag/`
Final filtering and annotation using the matrix from 2nd call.

```bash
bash VCF_FilterTag/FilTaggingToVCF.sh TagList F
```

---

## Job Management (LSF)

This pipeline depends on the LSF job scheduler and includes helper tools in:

```
somaticMutDetectTools/code/LSF_utils/
```

These include:
- `BsubS`, `BsubM`, `BsubL`: Job wrappers for different types
- `WaitSignalMaker.sh`: Job dependency tracker
- `File_pickByTag*.sh`: Extract paths by sample tag

---

## Reproducibility Strategy

Due to common instability in HPC environments:
- Each stage is checkpointed
- Logs (`*.log`, `*.err`) are preserved
- Only rerun failed steps as needed

See: `Recommended_Usage_Recovery.md` for more details.

---

## Citation

> Uryu, H. et al. *Clonal evolution of hematopoietic stem cells after autologous stem cell transplantation.*  
> *Nature Genetics* (2025). [https://doi.org/10.1038/s41588-025-02235-w](https://doi.org/10.1038/s41588-025-02235-w)

---

## Contact

**Hidetaka Uryu**  
📧 huryu@mdanderson.org

Please feel free to reach out for questions or collaborations.
