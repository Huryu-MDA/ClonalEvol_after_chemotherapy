# VCF_FilterTag — Final Annotation and Filtering of VCFs

This module applies the final layer of annotations and biological/technical filters to the output VCFs from Mutect2. These filters include low variant allele frequency (VAF), proximity to nearby mutations, low coverage, and more.

---

## 🔁 Typical Execution

```bash
FilTaggingToVCF.sh TagList_on_PID02 F
```

- `TagList_on_PID02`: A list of sample tags
- `F`: Gender (use `"M"` or `"F"`)

---

## 📂 Directory Contents

### Core Scripts

| Script | Description |
|--------|-------------|
| **FilTaggingToVCF.sh** | Launches the tagging pipeline per sample |
| **FilTaggingTo_Mu2filteredVCF.sh** | Splits VCF by chromosome, annotates, and filters each part |

### SubCode — Python Filtering Modules

Located in `SubCode/`:

| Script | Function |
|--------|----------|
| **VafF_CovF_Tag.py** | Adds tags based on variant allele frequency and sample-level low coverage |
| **dist10BP_Tag.py** | Identifies variants within 10bp of another SNV/INDEL |
| **FilTag_ConcatToVCFBody.py** | Merges filtered matrix with VCF body by chromosome-position-ref-alt |
| **LowCovTagOnSingleVCF.py** | Flags variants below coverage thresholds, customized for sex chromosomes |

---

## 🧬 Process Summary

1. **Filtered VCFs** from `Mu2_1stCall` are split by chromosome
2. **Filter matrices** from `Mu2_2ndCall` are read for corresponding chromosomes
3. Each line in the VCF is annotated with tags:
   - VAF filters
   - Distance-to-other-variant filter
   - Low coverage alert
4. Filtered and annotated VCFs are merged back into a per-sample VCF

---

## 📥 Input Requirements

- `TagList`: Sample tags that match filenames of `.merged.filtered.vcf`
- `FilMTX_onChr_list.txt`: List of filter matrix files from `Mu2_2ndCall`
- VCFs must reside under: `../Mu2Call_ToFMC_NoGNOMAD_gatk4200/*/*.merged.filtered.vcf`

---

## 📤 Output

- `TagVCF_withFilMTX/{Sample}/7FilTagged_{Sample}_chr*.merged.filtered.vcf`
- Chromosome-wise VCFs with added INFO tags
- Merged and compressed versions submitted to downstream modules

---

## Notes

- Uses `bsub`-based wrappers like `BsubM` and `BsubS`
- Dependent on `WaitSignalMaker.sh` for job coordination
- Gender-aware tagging is applied for sex chromosomes
