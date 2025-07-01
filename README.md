# Clonal Evolution of Hematopoietic Stem Cells After Cancer Chemotherapy

Welcome to the official repository accompanying our research article,\
**"Clonal evolution of hematopoietic stem cells after autologous stem cell transplantation."**

This repository contains code, scripts, and test data supporting our somatic mutation analysis, variant filtering, and phylogenetic reconstruction using patient-derived hematopoietic stem cell samples.

---

## Repository Structure

```
├── dir_list/                                 # Utility or directory listings
├── phylo_piechart_for_shared_variant_fraction/
│   ├── code/                                 # Scripts for pie chart generation
│   ├── data/                                 # Input data
│   └── test/                                 # Example/test data
├── phylo_plot_modify/
│   ├── code/                                 # Custom plot scripts for phylogenies
│   ├── data/                                 # Input data
│   ├── README.md                             # Module-specific documentation
│   └── test/                                 # Example/test outputs
├── somaticMutDetectTools/
│   ├── code/                                 # Main GATK-based variant calling pipeline
│   ├── README.md                             # Detailed guide for mutation detection
│   ├── test_data/                            # Sample input BAMs and metadata
│   └── test_run/                             # Expected output from test run
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

Navigate to `somaticMutDetectTools/code/` and run the main script:

```bash
bash Mutect2_pipeline.sh bamlist.txt SAMPLE_TAG
```

This will launch chromosome-wise Mutect2 jobs for the sample tag, and produce:

- Unfiltered VCFs per chromosome
- Merged unfiltered VCF
- Read orientation model
- Filtered final VCF: `${Tag}.merged.filtered.vcf`

Output is stored under:

```
Mu2Call_ToFMC_NoGNOMAD_gatk4200/SAMPLE_TAG/
```

---

## Additional Tools and Analysis

Other subdirectories contain code for:

- **Phylogenetic plotting**
- **Variant sharing visualization with pie charts**
- **Post-processing and figure generation for publication**

Please refer to each subfolder's `README.md` for module-specific documentation and test data.

---

## Citation

If you use this code or data in your research, please cite:

> Uryu, H., et al. "Clonal evolution of hematopoietic stem cells after autologous stem cell transplantation." *[https://doi.org/10.1038/s41588-025-02235-w]*

---

## Contact

For any questions, bug reports, or collaboration inquiries, please contact:

**Hidetaka Uryu**\
📧 [huryu@mdanderson.org](mailto\:huryu@mdanderson.org)

We welcome contributions and feedback to improve this research resource.

