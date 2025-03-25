#!/usr/bin/env python3
"""
vcf2maf.py - Convert a VCF file to MAF format

This script converts a VCF (Variant Call Format) file to MAF (Mutation Annotation Format),
with each variant annotated to only one of all possible gene isoforms.
Python implementation of the original Perl-based vcf2maf tool from MSKCC.

Original: https://github.com/mskcc/vcf2maf
"""

import os
import sys
import argparse
import subprocess
import tempfile
import logging
import re
import csv
import gzip
from collections import defaultdict
from datetime import datetime
import vcf  # PyVCF3
import pandas as pd
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('vcf2maf')

# MAF column headers based on the specification
MAF_COLUMNS = [
    'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome',
    'Start_Position', 'End_Position', 'Strand', 'Variant_Classification',
    'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
    'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
    'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1',
    'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2',
    'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase',
    'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer',
    'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short',
    'Transcript_ID', 'Exon_Number', 'n_depth', 't_depth', 'n_ref_count',
    'n_alt_count', 't_ref_count', 't_alt_count', 'all_effects'
]

# Variant classifications mapping
VARIANT_CLASSIFICATIONS = {
    'Frame_Shift_Del': 'frameshift deletion',
    'Frame_Shift_Ins': 'frameshift insertion',
    'In_Frame_Del': 'nonframeshift deletion',
    'In_Frame_Ins': 'nonframeshift insertion',
    'Missense_Mutation': 'missense',
    'Nonsense_Mutation': 'nonsense',
    'Silent': 'silent',
    'Splice_Site': 'splice site',
    'Translation_Start_Site': 'translation start site',
    'Nonstop_Mutation': 'nonstop',
    'Three_Prime_UTR': "3'UTR",
    'Five_Prime_UTR': "5'UTR",
    'Intron': 'intron',
    'RNA': 'RNA',
    'Targeted_Region': 'targeted region',
    'IGR': 'IGR',
    'Splice_Region': 'splice region',
    'Splice_Donor': 'splice donor',
    'Splice_Acceptor': 'splice acceptor',
    'Synonymous': 'synonymous',
    'Stop_Codon_Del': 'stop codon deletion',
    'Stop_Codon_Ins': 'stop codon insertion'
}

# Variant type mapping
VARIANT_TYPES = {
    'SNP': 'SNP',
    'DNP': 'DNP',
    'TNP': 'TNP',
    'ONP': 'ONP',
    'INS': 'INS',
    'DEL': 'DEL',
    'COMPLEX': 'COMPLEX'
}

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Convert VCF files to MAF format')
    
    # Required arguments
    parser.add_argument('--input-vcf', required=True, help='Path to input VCF file')
    parser.add_argument('--output-maf', required=True, help='Path to output MAF file')
    
    # Optional arguments
    parser.add_argument('--tumor-id', help='Tumor sample ID in the VCF')
    parser.add_argument('--normal-id', help='Normal sample ID in the VCF')
    parser.add_argument('--vcf-tumor-id', help='Tumor sample ID used in VCF genotype columns')
    parser.add_argument('--vcf-normal-id', help='Normal sample ID used in VCF genotype columns')
    parser.add_argument('--vep-path', help='Path to VEP installation directory')
    parser.add_argument('--vep-data', help='Path to VEP data directory')
    parser.add_argument('--vep-forks', type=int, default=4, help='Number of parallel processes to run VEP')
    parser.add_argument('--ref-fasta', help='Path to reference FASTA file')
    parser.add_argument('--species', default='homo_sapiens', help='Species name for VEP')
    parser.add_argument('--ncbi-build', default='GRCh38', help='NCBI build version')
    parser.add_argument('--cache-version', help='Version of VEP cache to use')
    parser.add_argument('--maf-center', default='Unknown', help='Center name to report in MAF')
    parser.add_argument('--retain-info', help='Comma-separated list of INFO fields to retain')
    parser.add_argument('--retain-fmt', help='Comma-separated list of FORMAT fields to retain')
    parser.add_argument('--min-hom-vaf', type=float, default=0.7, help='Minimum VAF to call a variant homozygous')
    parser.add_argument('--inhibit-vep', action='store_true', help='Skip running VEP')
    parser.add_argument('--verbose', action='store_true', help='Display more detailed progress')
    
    return parser.parse_args()

def check_prerequisites(args):
    """Check if required files and tools are available"""
    # Check if input VCF exists
    if not os.path.exists(args.input_vcf):
        logger.error(f"Input VCF file not found: {args.input_vcf}")
        sys.exit(1)
    
    # Check if VEP is needed and available
    if not args.inhibit_vep:
        if args.vep_path:
            vep_script = os.path.join(args.vep_path, 'vep')
            if not os.path.exists(vep_script):
                logger.error(f"VEP script not found at: {vep_script}")
                sys.exit(1)
        else:
            # Try to find VEP in PATH
            try:
                subprocess.run(['vep', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
            except FileNotFoundError:
                logger.error("VEP not found in PATH. Please install VEP or specify --vep-path")
                sys.exit(1)
    
    # Check if reference FASTA is needed and available
    if args.ref_fasta and not os.path.exists(args.ref_fasta):
        logger.error(f"Reference FASTA file not found: {args.ref_fasta}")
        sys.exit(1)

def run_vep(args, vcf_file):
    """Run VEP on the input VCF file"""
    if args.inhibit_vep:
        logger.info("Skipping VEP annotation as requested")
        return vcf_file
    
    logger.info("Running VEP annotation...")
    
    # Create a temporary file for VEP output
    vep_output = tempfile.NamedTemporaryFile(delete=False, suffix='.vcf')
    vep_output.close()
    
    # Build VEP command
    vep_cmd = ['vep']
    if args.vep_path:
        vep_cmd = [os.path.join(args.vep_path, 'vep')]
    
    vep_cmd.extend([
        '--input_file', vcf_file,
        '--output_file', vep_output.name,
        '--format', 'vcf',
        '--vcf',
        '--symbol',
        '--terms', 'SO',
        '--tsl',
        '--hgvs',
        '--fasta', args.ref_fasta if args.ref_fasta else '.',
        '--offline',
        '--cache',
        '--pick',
        '--fork', str(args.vep_forks),
        '--species', args.species
    ])
    
    if args.vep_data:
        vep_cmd.extend(['--dir', args.vep_data])
    
    if args.cache_version:
        vep_cmd.extend(['--cache_version', args.cache_version])
    
    # Run VEP
    try:
        logger.info(f"Running VEP command: {' '.join(vep_cmd)}")
        subprocess.run(vep_cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running VEP: {e}")
        sys.exit(1)
    
    return vep_output.name

def determine_variant_type(ref, alt):
    """Determine the variant type based on reference and alternate alleles"""
    if len(ref) == 1 and len(alt) == 1:
        return 'SNP'
    elif len(ref) == 2 and len(alt) == 2:
        return 'DNP'
    elif len(ref) == 3 and len(alt) == 3:
        return 'TNP'
    elif len(ref) > 3 and len(alt) > 3 and len(ref) == len(alt):
        return 'ONP'
    elif len(ref) > len(alt):
        return 'DEL'
    elif len(ref) < len(alt):
        return 'INS'
    else:
        return 'COMPLEX'

def extract_vep_data(info_field):
    """Extract VEP annotation data from the INFO field"""
    vep_data = {}
    
    # Look for the CSQ field which contains VEP annotations
    csq_match = re.search(r'CSQ=([^;]+)', info_field)
    if not csq_match:
        return vep_data
    
    csq_data = csq_match.group(1)
    
    # Parse the CSQ format from the VCF header
    # For simplicity, we'll use a default format if not available
    default_csq_format = [
        'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type',
        'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position',
        'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation',
        'DISTANCE', 'STRAND', 'FLAGS', 'SYMBOL_SOURCE', 'HGNC_ID'
    ]
    
    # Get the first CSQ entry (the picked one if --pick was used with VEP)
    csq_values = csq_data.split(',')[0].split('|')
    
    # Map values to their field names
    for i, field in enumerate(default_csq_format):
        if i < len(csq_values):
            vep_data[field] = csq_values[i]
    
    return vep_data

def map_variant_classification(consequence):
    """Map VEP consequence terms to MAF variant classifications"""
    # Default mapping for common terms
    consequence_map = {
        'frameshift_variant': 'Frame_Shift_Del',  # Will be adjusted based on variant type
        'stop_gained': 'Nonsense_Mutation',
        'stop_lost': 'Nonstop_Mutation',
        'start_lost': 'Translation_Start_Site',
        'initiator_codon_variant': 'Translation_Start_Site',
        'missense_variant': 'Missense_Mutation',
        'protein_altering_variant': 'Missense_Mutation',
        'coding_sequence_variant': 'Missense_Mutation',
        'conservative_missense_variant': 'Missense_Mutation',
        'rare_amino_acid_variant': 'Missense_Mutation',
        'transcript_ablation': 'Splice_Site',
        'splice_acceptor_variant': 'Splice_Site',
        'splice_donor_variant': 'Splice_Site',
        'splice_region_variant': 'Splice_Region',
        'synonymous_variant': 'Silent',
        'stop_retained_variant': 'Silent',
        'start_retained_variant': 'Silent',
        '3_prime_UTR_variant': 'Three_Prime_UTR',
        '5_prime_UTR_variant': 'Five_Prime_UTR',
        'intron_variant': 'Intron',
        'intergenic_variant': 'IGR',
        'upstream_gene_variant': 'IGR',
        'downstream_gene_variant': 'IGR',
        'TF_binding_site_variant': 'IGR',
        'regulatory_region_variant': 'IGR'
    }
    
    # Handle multiple consequences (comma-separated)
    if ',' in consequence:
        consequences = consequence.split(',')
        # Use the first consequence that has a mapping
        for cons in consequences:
            if cons in consequence_map:
                return consequence_map[cons]
    elif consequence in consequence_map:
        return consequence_map[consequence]
    
    # Default to Missense_Mutation if no mapping found
    return 'Missense_Mutation'

def adjust_variant_classification(var_class, var_type):
    """Adjust variant classification based on variant type"""
    if var_class == 'Frame_Shift_Del' and var_type in ['INS', 'COMPLEX']:
        return 'Frame_Shift_Ins'
    elif var_class == 'Frame_Shift_Ins' and var_type in ['DEL', 'COMPLEX']:
        return 'Frame_Shift_Del'
    return var_class

def parse_vcf_to_maf(args):
    """Parse VCF file and convert to MAF format"""
    logger.info(f"Converting VCF file {args.input_vcf} to MAF format...")
    
    # Run VEP if needed
    vcf_file = args.input_vcf
    if not args.inhibit_vep:
        vcf_file = run_vep(args, vcf_file)
    
    # Open VCF file
    vcf_reader = vcf.Reader(filename=vcf_file)
    
    # Determine tumor and normal sample IDs
    tumor_id = args.tumor_id
    normal_id = args.normal_id
    vcf_tumor_id = args.vcf_tumor_id or tumor_id
    vcf_normal_id = args.vcf_normal_id or normal_id
    
    # If tumor_id not specified, use the first sample in the VCF
    if not tumor_id and vcf_reader.samples:
        tumor_id = vcf_reader.samples[0]
        vcf_tumor_id = tumor_id
        logger.info(f"No tumor ID specified, using first sample in VCF: {tumor_id}")
    
    # If normal_id not specified but there are multiple samples, use the second one
    if not normal_id and len(vcf_reader.samples) > 1:
        normal_id = vcf_reader.samples[1]
        vcf_normal_id = normal_id
        logger.info(f"No normal ID specified, using second sample in VCF: {normal_id}")
    
    # Prepare MAF data
    maf_data = []
    
    # Process each variant in the VCF
    for record in tqdm(vcf_reader, desc="Processing variants"):
        # For each alternate allele
        for alt_idx, alt in enumerate(record.ALT):
            # Skip if alt is None or not a proper allele
            if alt is None:
                continue
            
            alt_allele = str(alt)
            ref_allele = str(record.REF)
            
            # Determine variant type
            var_type = determine_variant_type(ref_allele, alt_allele)
            
            # Extract genotype information
            tumor_gt = None
            normal_gt = None
            t_depth = 0
            n_depth = 0
            t_ref_count = 0
            t_alt_count = 0
            n_ref_count = 0
            n_alt_count = 0
            
            # Get tumor genotype
            if vcf_tumor_id in record.samples:
                tumor_sample = record.samples[vcf_tumor_id]
                tumor_gt = tumor_sample.gt_type  # 0=hom_ref, 1=het, 2=hom_alt, None=missing
                
                # Extract depth and allele counts if available
                if hasattr(tumor_sample.data, 'DP'):
                    t_depth = tumor_sample.data.DP
                
                # Try to get ref and alt counts from various FORMAT fields
                if hasattr(tumor_sample.data, 'AD') and tumor_sample.data.AD:
                    t_ref_count = tumor_sample.data.AD[0]
                    if len(tumor_sample.data.AD) > alt_idx + 1:
                        t_alt_count = tumor_sample.data.AD[alt_idx + 1]
            
            # Get normal genotype
            if vcf_normal_id in record.samples:
                normal_sample = record.samples[vcf_normal_id]
                normal_gt = normal_sample.gt_type
                
                # Extract depth and allele counts if available
                if hasattr(normal_sample.data, 'DP'):
                    n_depth = normal_sample.data.DP
                
                # Try to get ref and alt counts from various FORMAT fields
                if hasattr(normal_sample.data, 'AD') and normal_sample.data.AD:
                    n_ref_count = normal_sample.data.AD[0]
                    if len(normal_sample.data.AD) > alt_idx + 1:
                        n_alt_count = normal_sample.data.AD[alt_idx + 1]
            
            # Extract VEP data if available
            vep_data = extract_vep_data(str(record.INFO))
            
            # Determine variant classification
            var_class = 'Missense_Mutation'  # Default
            if 'Consequence' in vep_data:
                var_class = map_variant_classification(vep_data['Consequence'])
                var_class = adjust_variant_classification(var_class, var_type)
            
            # Prepare MAF record
            maf_record = {
                'Hugo_Symbol': vep_data.get('SYMBOL', 'Unknown'),
                'Entrez_Gene_Id': vep_data.get('Gene', '0'),
                'Center': args.maf_center,
                'NCBI_Build': args.ncbi_build,
                'Chromosome': record.CHROM,
                'Start_Position': record.POS,
                'End_Position': record.POS + len(ref_allele) - 1,
                'Strand': '+',
                'Variant_Classification': var_class,
                'Variant_Type': var_type,
                'Reference_Allele': ref_allele,
                'Tumor_Seq_Allele1': ref_allele if tumor_gt in [0, 1] else alt_allele,
                'Tumor_Seq_Allele2': alt_allele if tumor_gt in [1, 2] else ref_allele,
                'dbSNP_RS': record.ID if record.ID else '',
                'dbSNP_Val_Status': '',
                'Tumor_Sample_Barcode': tumor_id,
                'Matched_Norm_Sample_Barcode': normal_id if normal_id else '',
                'Match_Norm_Seq_Allele1': ref_allele if normal_gt in [0, 1] else alt_allele,
                'Match_Norm_Seq_Allele2': alt_allele if normal_gt in [1, 2] else ref_allele,
                'Tumor_Validation_Allele1': '',
                'Tumor_Validation_Allele2': '',
                'Match_Norm_Validation_Allele1': '',
                'Match_Norm_Validation_Allele2': '',
                'Verification_Status': '',
                'Validation_Status': '',
                'Mutation_Status': 'Somatic' if normal_gt == 0 and tumor_gt in [1, 2] else 'Unknown',
                'Sequencing_Phase': '',
                'Sequence_Source': '',
                'Validation_Method': '',
                'Score': '',
                'BAM_File': '',
                'Sequencer': '',
                'Tumor_Sample_UUID': '',
                'Matched_Norm_Sample_UUID': '',
                'HGVSc': vep_data.get('HGVSc', ''),
                'HGVSp': vep_data.get('HGVSp', ''),
                'HGVSp_Short': vep_data.get('HGVSp', '').replace('%3D', '=').split(':')[-1] if 'HGVSp' in vep_data else '',
                'Transcript_ID': vep_data.get('Feature', ''),
                'Exon_Number': vep_data.get('EXON', ''),
                'n_depth': n_depth,
                't_depth': t_depth,
                'n_ref_count': n_ref_count,
                'n_alt_count': n_alt_count,
                't_ref_count': t_ref_count,
                't_alt_count': t_alt_count,
                'all_effects': ''
            }
            
            maf_data.append(maf_record)
    
    # Write MAF file
    with open(args.output_maf, 'w', newline='') as maf_file:
        writer = csv.DictWriter(maf_file, fieldnames=MAF_COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(maf_data)
    
    logger.info(f"Conversion complete. MAF file written to {args.output_maf}")
    logger.info(f"Total variants processed: {len(maf_data)}")
    
    # Clean up temporary files
    if not args.inhibit_vep and vcf_file != args.input_vcf:
        os.unlink(vcf_file)

def main():
    """Main function"""
    # Parse command line arguments
    args = parse_arguments()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Check prerequisites
    check_prerequisites(args)
    
    # Convert VCF to MAF
    parse_vcf_to_maf(args)

if __name__ == '__main__':
    main()
