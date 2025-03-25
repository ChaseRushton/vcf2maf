# vcf2maf (Python Implementation)

A Python implementation of the [mskcc/vcf2maf](https://github.com/mskcc/vcf2maf) tool, which converts VCF (Variant Call Format) files to MAF (Mutation Annotation Format) files.

## Overview

This tool converts a VCF file into a MAF file, where each variant is annotated to only one of all possible gene isoforms. The implementation follows the functionality of the original Perl-based vcf2maf tool but is written entirely in Python.

## Features

- Convert VCF files to MAF format
- Annotate variants using Ensembl VEP (Variant Effect Predictor)
- Support for tumor/normal sample pairs
- Proper handling of different VCF formats from various variant callers
- Option to convert MAF files back to VCF format

## Requirements

- Python 3.6+
- Ensembl VEP (optional, for variant annotation)
- Python packages:
  - pysam
  - pandas
  - biopython

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```bash
python vcf2maf.py --input-vcf input.vcf --output-maf output.maf
```

### With Tumor/Normal Sample IDs

```bash
python vcf2maf.py --input-vcf input.vcf --output-maf output.maf --tumor-id TUMOR_SAMPLE --normal-id NORMAL_SAMPLE
```

### Using VEP for Annotation

```bash
python vcf2maf.py --input-vcf input.vcf --output-maf output.maf --vep-path /path/to/vep --vep-data /path/to/vep/data
```

### Converting MAF to VCF

```bash
python maf2vcf.py --input-maf input.maf --output-vcf output.vcf
```


