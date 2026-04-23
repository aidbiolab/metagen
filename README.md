# Metagen: Shotgun Metagenomics Pipeline for Oral Samples


**Author Pipeline:** Fabian Tobar-Tosse  
**Institution:** Pontificia Universidad Javeriana Cali, and Aidbio (Bioinformatics and Data Science Solutions)  
**ORCID:** 0000-0001-5334-4286  
**Contact:** ftobar@javerianacali.edu.co, https://www.aidbio.com.co/

**Author related Paper:** Alveiro Erira, et al.,  
**Institution:** Universidad Cooperativa de Colombia, Bogotá D.C   
**ORCID:** 0000-0003-1509-5631  
**Contact:** alveiro.erira@campusucc.edu.co

---

## Overview

**Metagen** is a modular bioinformatics pipeline designed for the analysis of **shotgun metagenomic sequencing data** (unpaired FASTQ files) from oral samples (dental plaque, saliva, and tumor tissue).

This pipeline was developed as part of the study:

> **"Bacterial metagenome in plaque, saliva, and tumor samples from individuals with and without oral squamous cell carcinoma by Next Generation Sequencing"**  
> (Manuscript in preparation for *GigaByte*, 2026)

**Data Availability:**  
NCBI BioProject: [PRJNA1165312](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1165312)

---

## Pipeline Features

- Quality control and trimming (`fastp` + `FastQC` + `MultiQC`)
- Taxonomic classification (`Kraken2`)
- Abundance estimation (`Bracken`)
- Interactive visualization (`Krona`)
- Alpha and Beta diversity analysis (`scikit-bio`)
- Generation of publication-ready figures (heatmaps, PCoA, Venn diagrams)

---

## Folder Structure

```bash
metagen_puj/
├── data/
│   ├── MetsNCBI/                 # Raw .fastq.gz files
│   └── kraken/                   # Kraken2 database (k2_standard_20260226)
├── outputs/
│   ├── fastp/
│   ├── fastqc/
│   ├── multiqc/
│   ├── kraken2/
│   ├── bracken/
│   ├── krona/
│   └── diversity/
├── pipeline_meta.py              # Main pipeline script
├── README.md
└── environment.yml               # (optional) Conda environment
```
## Requirements
- Linux server (tested on 32 cores, 120 GB RAM)
- Conda / Mamba

### Create the environment:

```bash
mamba env create -f environment.yml
conda activate metagen-qc
```
### Or manually:
```bash
mamba create -n metagen-qc -c conda-forge -c bioconda \
  fastp fastqc multiqc pigz kraken2 bracken krona pandas scikit-bio matplotlib seaborn -y
```
### Usage
1. Place your raw FASTQ files in ./data/MetsNCBI/
2. Download and extract the Kraken2 database into ./data/kraken/k2_standard_20260226/
3. Edit pipeline_meta.py to uncomment the steps you want to run
4. Execute:
```bash
conda activate metagen-qc
python3 pipeline_meta.py
```
### Main Outputs
- outputs/bracken/merged_rel_abund.csv → Merged abundance table
- outputs/diversity/alpha_diversity_counts.csv
- outputs/diversity/pcoa_braycurtis.csv
- outputs/krona/*.krona.html → Interactive Krona charts
### Citation
If you use this pipeline in your research, please cite:
- Erira A, et al. (2026). Bacterial metagenome in plaque, saliva, and tumor samples from individuals with and without oral squamous cell carcinoma. GigaByte (in preparation).
- GitHub repository: https://github.com/aidbiolab/metagen
### License
This project is licensed under the MIT License — feel free to use and modify it.
