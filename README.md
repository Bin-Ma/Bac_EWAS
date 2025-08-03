# 🧬 EWAS Analysis Pipeline

A Python script for **Epigenome-Wide Association Study (EWAS)** using bacterial pan-methylome data and gene annotations. It generates genotype and phenotype data in [PLINK](https://zzz.bwh.harvard.edu/plink/) format to facilitate downstream association analyses of methylation and strain adaptation.

---

##  Project Structure

```
.
├── ewas_analysis.py      # Main analysis script
├── README.md             # This help documentation
├── data/                 # Example input files (optional)
└── output/               # Output directory for PLINK files
```

---

## Requirements

- Python 3.6+
- Dependencies:

```bash
pip install pandas numpy
```

---

## Quick Start

```bash
python ewas_analysis.py <index> <gene_index> <out> <methy_info> <gene_info> <classA> \
    [--methy_name 6mA|5mC] [--methy_type complete|pos|neg] [--gene_type G|N]
```

---

## Positional Arguments

| Argument       | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `index`        | Path to strain list (FASTA or identifier list for methylation presence).    |
| `gene_index`   | Gene index file (used to remove `.fa` suffixes).                            |
| `out`          | Output prefix or directory for results.                                     |
| `methy_info`   | Methylation site data (CSV) with per-strain methylation status in column 8. |
| `gene_info`    | Gene annotation table with columns: `Gene`, `start`, `length`.              |
| `classA`       | List of strain names adapted to the environment (treated as fitness = 1).   |

---

## Optional Arguments

| Flag             | Default   | Description                                                              |
|------------------|-----------|---------------------------------------------------------------------------|
| `--methy_name`   | `6mA`     | Type of methylation mark: `6mA` or `5mC`.                                 |
| `--methy_type`   | `complete`| Site filtering strategy: `complete`, `pos`, or `neg`.                    |
| `--gene_type`    | `G`       | Prefix for gene type: `G` for gene, `B` for non-coding.                  |

---

##  Output Files

| File              | Description                                                                          |
|-------------------|--------------------------------------------------------------------------------------|
| `plink.ped`       | PED file containing methylation genotype matrix and binary phenotype column.         |
| `plink.map`       | MAP file with methylation site locations (based on gene start + site position).      |

---

##  Example

```bash
python ewas_analysis.py \
    strain_index.csv \
    gene_index.csv \
    output/ \
    methylation_data.csv \
    gene_annotation.csv \
    fit_strains.csv \
    --methy_name 6mA \
    --methy_type complete \
    --gene_type G
```

---

##  Notes

- The script encodes methylation status using artificial "alleles" such as `AMMA`, `CNNC`, `AAAA`, etc., for compatibility with PLINK.
- Strains in `classA` are assigned **phenotype = 1 (adapted)**, all others are assigned **phenotype = 2 (non-adapted)**.
- Site matching occurs only if the methylation site is within a valid gene region based on annotation.

---

##  Input File Format Overview

### `methy_info` CSV (example columns)

| Gene_ID | Site | ... | Methylated_Strain_List |
|---------|------|-----|-------------------------|
| 00001   | 125  | ... | strainA,strainC,...     |

### `gene_info` CSV

| Gene     | start | end |  length  |
|----------|-------|--------|----|
| G_00001  | 1000  | 1500    |    500|

### `classA` (TXT)

```
strainA
strainB
strainF
```

---

##  License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

---




## ✉️ Contact

For questions or feedback, please contact:  
📧 `your_email@example.com`
