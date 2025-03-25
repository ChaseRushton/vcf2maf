#!/usr/bin/env python3
"""
maf2vcf.py - Convert a MAF file to VCF format

This script converts a MAF (Mutation Annotation Format) file to VCF (Variant Call Format),
preserving as much information as possible.
Python implementation of the original Perl-based maf2vcf tool from MSKCC.

Original: https://github.com/mskcc/vcf2maf
"""

import os
import sys
import argparse
import logging
import csv
import re
from collections import defaultdict
import pandas as pd
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('maf2vcf')

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Convert MAF files to VCF format')
    
    # Required arguments
    parser.add_argument('--input-maf', required=True, help='Path to input MAF file')
    parser.add_argument('--output-vcf', required=True, help='Path to output VCF file')
    
    # Optional arguments
    parser.add_argument('--ref-fasta', help='Path to reference FASTA file')
    parser.add_argument('--tumor-id', help='Tumor sample ID to use in the VCF')
    parser.add_argument('--normal-id', help='Normal sample ID to use in the VCF')
    parser.add_argument('--vcf-tumor-id', help='Tumor sample ID to use in VCF genotype columns')
    parser.add_argument('--vcf-normal-id', help='Normal sample ID to use in VCF genotype columns')
    parser.add_argument('--per-tn-vcfs', action='store_true', help='Create a VCF for each tumor-normal pair')
    parser.add_argument('--tum-depth-col', default='t_depth', help='Column name for tumor depth in MAF')
    parser.add_argument('--tum-rad-col', default='t_ref_count', help='Column name for tumor reference allele depth in MAF')
    parser.add_argument('--tum-vad-col', default='t_alt_count', help='Column name for tumor variant allele depth in MAF')
    parser.add_argument('--nrm-depth-col', default='n_depth', help='Column name for normal depth in MAF')
    parser.add_argument('--nrm-rad-col', default='n_ref_count', help='Column name for normal reference allele depth in MAF')
    parser.add_argument('--nrm-vad-col', default='n_alt_count', help='Column name for normal variant allele depth in MAF')
    parser.add_argument('--verbose', action='store_true', help='Display more detailed progress')
    
    return parser.parse_args()

def check_prerequisites(args):
    """Check if required files and tools are available"""
    # Check if input MAF exists
    if not os.path.exists(args.input_maf):
        logger.error(f"Input MAF file not found: {args.input_maf}")
        sys.exit(1)
    
    # Check if reference FASTA is needed and available
    if args.ref_fasta and not os.path.exists(args.ref_fasta):
        logger.error(f"Reference FASTA file not found: {args.ref_fasta}")
        sys.exit(1)

def normalize_chromosome(chrom):
    """Normalize chromosome names to ensure consistency"""
    # Remove 'chr' prefix if present
    if chrom.lower().startswith('chr'):
        chrom = chrom[3:]
    
    # Convert to string and handle special cases
    if chrom == '23':
        return 'X'
    elif chrom == '24':
        return 'Y'
    elif chrom == '25':
        return 'MT'
    else:
        return str(chrom)

def get_vcf_genotype_format():
    """Define the FORMAT field for the VCF file"""
    return ['GT', 'AD', 'DP', 'AF']

def get_vcf_header(args, samples):
    """Generate VCF header"""
    header = [
        '##fileformat=VCFv4.2',
        f'##fileDate={pd.Timestamp.now().strftime("%Y%m%d")}',
        '##source=maf2vcf.py',
        '##reference=GRCh38',  # Default, can be modified based on MAF
        '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">',
        '##INFO=<ID=VT,Number=1,Type=String,Description="Variant type">',
        '##INFO=<ID=Gene,Number=1,Type=String,Description="Gene symbol">',
        '##INFO=<ID=Transcript,Number=1,Type=String,Description="Transcript ID">',
        '##INFO=<ID=HGVSc,Number=1,Type=String,Description="HGVS coding sequence name">',
        '##INFO=<ID=HGVSp,Number=1,Type=String,Description="HGVS protein sequence name">',
        '##INFO=<ID=VarClass,Number=1,Type=String,Description="Variant classification">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">'
    ]
    
    # Add the column header line
    header.append(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}')
    
    return header

def parse_maf_to_vcf(args):
    """Parse MAF file and convert to VCF format"""
    logger.info(f"Converting MAF file {args.input_maf} to VCF format...")
    
    # Read MAF file
    try:
        maf_df = pd.read_csv(args.input_maf, sep='\t', comment='#', low_memory=False)
    except Exception as e:
        logger.error(f"Error reading MAF file: {e}")
        sys.exit(1)
    
    # Check required columns
    required_columns = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']
    for col in required_columns:
        if col not in maf_df.columns:
            logger.error(f"Required column '{col}' not found in MAF file")
            sys.exit(1)
    
    # Determine tumor and normal sample IDs
    tumor_id = args.tumor_id
    normal_id = args.normal_id
    vcf_tumor_id = args.vcf_tumor_id or tumor_id or 'TUMOR'
    vcf_normal_id = args.vcf_normal_id or normal_id or 'NORMAL'
    
    # If tumor_id not specified, try to get it from the MAF
    if not tumor_id and 'Tumor_Sample_Barcode' in maf_df.columns:
        tumor_id = maf_df['Tumor_Sample_Barcode'].iloc[0]
        if not args.vcf_tumor_id:
            vcf_tumor_id = tumor_id
        logger.info(f"No tumor ID specified, using first sample in MAF: {tumor_id}")
    
    # If normal_id not specified, try to get it from the MAF
    if not normal_id and 'Matched_Norm_Sample_Barcode' in maf_df.columns:
        normal_id = maf_df['Matched_Norm_Sample_Barcode'].iloc[0]
        if not args.vcf_normal_id:
            vcf_normal_id = normal_id
        logger.info(f"No normal ID specified, using first sample in MAF: {normal_id}")
    
    # Define sample string for VCF header
    samples = vcf_tumor_id
    if normal_id:
        samples = f"{vcf_normal_id}\t{vcf_tumor_id}"
    
    # Generate VCF header
    vcf_header = get_vcf_header(args, samples)
    
    # Prepare VCF records
    vcf_records = []
    
    # Process each variant in the MAF
    for _, row in tqdm(maf_df.iterrows(), total=len(maf_df), desc="Processing variants"):
        chrom = normalize_chromosome(str(row['Chromosome']))
        pos = int(row['Start_Position'])
        ref_allele = row['Reference_Allele']
        alt_allele = row['Tumor_Seq_Allele2']
        
        # Skip if ref and alt are the same
        if ref_allele == alt_allele:
            continue
        
        # Get variant ID (dbSNP if available)
        variant_id = '.'
        if 'dbSNP_RS' in row and pd.notna(row['dbSNP_RS']) and row['dbSNP_RS'] != '':
            variant_id = row['dbSNP_RS']
        
        # Prepare INFO field
        info_fields = []
        
        # Add somatic flag if mutation status is available
        if 'Mutation_Status' in row and row['Mutation_Status'] == 'Somatic':
            info_fields.append('SOMATIC')
        
        # Add variant type
        if 'Variant_Type' in row:
            info_fields.append(f'VT={row["Variant_Type"]}')
        
        # Add gene symbol
        if 'Hugo_Symbol' in row and pd.notna(row['Hugo_Symbol']):
            info_fields.append(f'Gene={row["Hugo_Symbol"]}')
        
        # Add transcript ID
        if 'Transcript_ID' in row and pd.notna(row['Transcript_ID']):
            info_fields.append(f'Transcript={row["Transcript_ID"]}')
        
        # Add HGVS notation
        if 'HGVSc' in row and pd.notna(row['HGVSc']):
            info_fields.append(f'HGVSc={row["HGVSc"]}')
        if 'HGVSp' in row and pd.notna(row['HGVSp']):
            info_fields.append(f'HGVSp={row["HGVSp"]}')
        
        # Add variant classification
        if 'Variant_Classification' in row:
            info_fields.append(f'VarClass={row["Variant_Classification"]}')
        
        # Combine INFO fields
        info = ';'.join(info_fields) if info_fields else '.'
        
        # Prepare FORMAT field
        format_fields = get_vcf_genotype_format()
        format_str = ':'.join(format_fields)
        
        # Prepare sample genotypes
        tumor_gt = {}
        normal_gt = {}
        
        # Extract tumor depth and allele counts
        t_depth = int(row[args.tum_depth_col]) if args.tum_depth_col in row and pd.notna(row[args.tum_depth_col]) else 0
        t_ref_count = int(row[args.tum_rad_col]) if args.tum_rad_col in row and pd.notna(row[args.tum_rad_col]) else 0
        t_alt_count = int(row[args.tum_vad_col]) if args.tum_vad_col in row and pd.notna(row[args.tum_vad_col]) else 0
        
        # Calculate tumor genotype
        if t_depth == 0:
            tumor_gt['GT'] = './.'
            tumor_gt['AD'] = '0,0'
            tumor_gt['DP'] = '0'
            tumor_gt['AF'] = '0'
        else:
            # Determine genotype based on allele counts
            if t_alt_count == 0:
                tumor_gt['GT'] = '0/0'  # Homozygous reference
            elif t_ref_count == 0:
                tumor_gt['GT'] = '1/1'  # Homozygous alternate
            else:
                tumor_gt['GT'] = '0/1'  # Heterozygous
            
            tumor_gt['AD'] = f'{t_ref_count},{t_alt_count}'
            tumor_gt['DP'] = str(t_depth)
            tumor_gt['AF'] = f'{t_alt_count/t_depth:.4f}' if t_depth > 0 else '0'
        
        # Extract normal depth and allele counts if available
        if normal_id:
            n_depth = int(row[args.nrm_depth_col]) if args.nrm_depth_col in row and pd.notna(row[args.nrm_depth_col]) else 0
            n_ref_count = int(row[args.nrm_rad_col]) if args.nrm_rad_col in row and pd.notna(row[args.nrm_rad_col]) else 0
            n_alt_count = int(row[args.nrm_vad_col]) if args.nrm_vad_col in row and pd.notna(row[args.nrm_vad_col]) else 0
            
            # Calculate normal genotype
            if n_depth == 0:
                normal_gt['GT'] = './.'
                normal_gt['AD'] = '0,0'
                normal_gt['DP'] = '0'
                normal_gt['AF'] = '0'
            else:
                # Determine genotype based on allele counts
                if n_alt_count == 0:
                    normal_gt['GT'] = '0/0'  # Homozygous reference
                elif n_ref_count == 0:
                    normal_gt['GT'] = '1/1'  # Homozygous alternate
                else:
                    normal_gt['GT'] = '0/1'  # Heterozygous
                
                normal_gt['AD'] = f'{n_ref_count},{n_alt_count}'
                normal_gt['DP'] = str(n_depth)
                normal_gt['AF'] = f'{n_alt_count/n_depth:.4f}' if n_depth > 0 else '0'
        
        # Format genotype strings
        tumor_gt_str = ':'.join([tumor_gt.get(f, '.') for f in format_fields])
        
        # Prepare the full VCF record
        if normal_id:
            normal_gt_str = ':'.join([normal_gt.get(f, '.') for f in format_fields])
            sample_gt = f"{normal_gt_str}\t{tumor_gt_str}"
        else:
            sample_gt = tumor_gt_str
        
        vcf_record = f"{chrom}\t{pos}\t{variant_id}\t{ref_allele}\t{alt_allele}\t.\tPASS\t{info}\t{format_str}\t{sample_gt}"
        vcf_records.append(vcf_record)
    
    # Write VCF file
    with open(args.output_vcf, 'w') as vcf_file:
        # Write header
        for header_line in vcf_header:
            vcf_file.write(f"{header_line}\n")
        
        # Write records
        for record in vcf_records:
            vcf_file.write(f"{record}\n")
    
    logger.info(f"Conversion complete. VCF file written to {args.output_vcf}")
    logger.info(f"Total variants processed: {len(vcf_records)}")

def main():
    """Main function"""
    # Parse command line arguments
    args = parse_arguments()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Check prerequisites
    check_prerequisites(args)
    
    # Convert MAF to VCF
    parse_maf_to_vcf(args)

if __name__ == '__main__':
    main()
