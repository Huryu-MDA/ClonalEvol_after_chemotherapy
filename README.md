# Clonal Evolution of Hematopoietic Stem Cells  After Autologous Stem Cell Transplantation

Welcome to the official repository accompanying our research article,\
**"Clonal Evolution of Hematopoietic Stem Cells  After Autologous Stem Cell Transplantation."**

This repository contains code, scripts, and test data supporting our somatic mutation analysis, variant filtering, and phylogenetic reconstruction using patient-derived hematopoietic stem cell samples.

---

## Repository Structure

```
├── dir_list/                                 # Utility or directory listings
├── phylo_piechart_for_shared_variant_fraction/
│   ├── code/                                 # Scripts to express the shared variant frequencies in each node of phylogenetic trees
│   ├── data/                                 # Input data
│   └── test/                                 # Example/test data
├── phylo_plot_modify/
│   ├── code/                                 # Custom plot scripts for phylogenies; tree branches can be converted to mutational signature compositions
│   ├── data/                                 # Input data
│   ├── README.md                             # Module-specific documentation
│   └── test/                                 # Example/test outputs
├── somaticMutDetectTools/
│   ├── code/                                 # Core somatic mutation pipeline (stepwise execution)
│   │   ├── Mu2_1stCall/                      # Initial Mutect2 per-chromosome calls
│   │   ├── Mu2_2ndCall/                      # Merge and filter VCFs; orientation model
│   │   └── VCF_FilterTag/                    # Custom postprocessing and tagging
│   ├── README.md                             # In-depth guide for using this module
│   ├── test_data/                            # BAM test inputs and reference resources
│   │   ├── PON_PANEL_on_hg19/                # Panel of Normals (chr-split VCFs)
│   │   └── test_data_on_tp53/                # BAM/metadata for TP53 test cases
│   └── test_run/                             # Results from running on test data
│       ├── README.md                         # Explains test procedure and expected outputs
│       ├── test_run_backup/                  # Archived or previous runs
│       └── test_run_results/                 # Output of test pipeline (VCFs, logs, etc.)
└── README.md                                 # This file
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

| File                            | Description                                                  |
| ------------------------------- | ------------------------------------------------------------ |
| `Homo_sapiens_assembly19.fasta` | Reference genome (hg19)                                      |
| `pon_chr*.vcf.gz`               | Panel of Normals for each chromosome                         |
| `bamlist.txt`                   | File containing full paths to BAM files with associated tags |
| *(optional)* `gnomAD.vcf`       | For contamination filtering (currently commented out)        |

---

## Running the Somatic Mutation Detection Pipeline

Navigate to `somaticMutDetectTools/code/` and execute the steps in order:

1. **Mu2\_1stCall**: Run initial per-chromosome somatic variant calls using GATK Mutect2
2. **Mu2\_2ndCall**: Merge outputs, learn read orientation model, and apply filtering
3. **VCF\_FilterTag**: Postprocess and annotate filtered VCFs as needed

Each subfolder contains scripts that can be executed independently or integrated in a workflow.

Example:

```bash
bash Mu2_1stCall/run_Mutect2_by_chr.sh bamlist.txt SAMPLE_TAG
bash Mu2_2ndCall/run_merge_filter.sh SAMPLE_TAG
bash VCF_FilterTag/apply_custom_filter.sh SAMPLE_TAG
```

Final results are stored under a subdirectory like:

```
Mu2Call_ToFMC_NoGNOMAD_gatk4200/SAMPLE_TAG/
```

---

## Additional Tools and Analysis

Other subdirectories contain code for:

- **Phylogenetic plotting** — Custom scripts for converting phylogenetic tree branches into mutational signature compositions
- **Variant sharing visualization with pie charts** — Scripts to express the shared variant frequencies in each node of phylogenetic trees
- **Post-processing and figure generation for publication**

Please refer to each subfolder's `README.md` for module-specific documentation and test data.

---

## Citation

If you use this code or data in your research, please cite:

> Uryu, H. et al. *Clonal evolution of hematopoietic stem cells after autologous stem cell transplantation*. *Nature Genetics*. [https://doi.org/10.1038/s41588-025-02235-w](https://doi.org/10.1038/s41588-025-02235-w) (2025).

---

## Contact

For any questions, bug reports, or collaboration inquiries, please contact:

**Hidetaka Uryu**\
📧 [huryu@mdanderson.org](mailto\:huryu@mdanderson.org)

We welcome contributions and feedback to improve this research resource.

