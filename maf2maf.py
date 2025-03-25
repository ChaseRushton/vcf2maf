#!/usr/bin/env python3
"""
maf2maf.py - Reannotate a MAF file

This script reannotates a MAF file by converting it to VCF format and then back to MAF format.
Python implementation of the original Perl-based maf2maf tool from MSKCC.

Original: https://github.com/mskcc/vcf2maf
"""

import os
import sys
import argparse
import logging
import tempfile
import subprocess
from pathlib import Path
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('maf2maf')

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Reannotate MAF files')
    
    # Required arguments
    parser.add_argument('--input-maf', required=True, help='Path to input MAF file')
    parser.add_argument('--output-maf', required=True, help='Path to output MAF file')
    
    # Optional arguments
    parser.add_argument('--tumor-id', help='Tumor sample ID to use')
    parser.add_argument('--normal-id', help='Normal sample ID to use')
    parser.add_argument('--vcf-tumor-id', help='Tumor sample ID to use in VCF genotype columns')
    parser.add_argument('--vcf-normal-id', help='Normal sample ID to use in VCF genotype columns')
    parser.add_argument('--ref-fasta', help='Path to reference FASTA file')
    parser.add_argument('--vep-path', help='Path to VEP installation directory')
    parser.add_argument('--vep-data', help='Path to VEP data directory')
    parser.add_argument('--vep-forks', type=int, default=4, help='Number of parallel processes to run VEP')
    parser.add_argument('--species', default='homo_sapiens', help='Species name for VEP')
    parser.add_argument('--ncbi-build', default='GRCh38', help='NCBI build version')
    parser.add_argument('--cache-version', help='Version of VEP cache to use')
    parser.add_argument('--maf-center', default='Unknown', help='Center name to report in MAF')
    parser.add_argument('--retain-info', help='Comma-separated list of INFO fields to retain')
    parser.add_argument('--retain-fmt', help='Comma-separated list of FORMAT fields to retain')
    parser.add_argument('--min-hom-vaf', type=float, default=0.7, help='Minimum VAF to call a variant homozygous')
    parser.add_argument('--tum-depth-col', default='t_depth', help='Column name for tumor depth in MAF')
    parser.add_argument('--tum-rad-col', default='t_ref_count', help='Column name for tumor reference allele depth in MAF')
    parser.add_argument('--tum-vad-col', default='t_alt_count', help='Column name for tumor variant allele depth in MAF')
    parser.add_argument('--nrm-depth-col', default='n_depth', help='Column name for normal depth in MAF')
    parser.add_argument('--nrm-rad-col', default='n_ref_count', help='Column name for normal reference allele depth in MAF')
    parser.add_argument('--nrm-vad-col', default='n_alt_count', help='Column name for normal variant allele depth in MAF')
    parser.add_argument('--inhibit-vep', action='store_true', help='Skip running VEP')
    parser.add_argument('--verbose', action='store_true', help='Display more detailed progress')
    
    return parser.parse_args()

def check_prerequisites(args):
    """Check if required files and tools are available"""
    # Check if input MAF exists
    if not os.path.exists(args.input_maf):
        logger.error(f"Input MAF file not found: {args.input_maf}")
        sys.exit(1)
    
    # Check if maf2vcf.py and vcf2maf.py are available
    script_dir = os.path.dirname(os.path.abspath(__file__))
    maf2vcf_path = os.path.join(script_dir, 'maf2vcf.py')
    vcf2maf_path = os.path.join(script_dir, 'vcf2maf.py')
    
    if not os.path.exists(maf2vcf_path):
        logger.error(f"maf2vcf.py not found at: {maf2vcf_path}")
        sys.exit(1)
    
    if not os.path.exists(vcf2maf_path):
        logger.error(f"vcf2maf.py not found at: {vcf2maf_path}")
        sys.exit(1)
    
    # Check if reference FASTA is needed and available
    if args.ref_fasta and not os.path.exists(args.ref_fasta):
        logger.error(f"Reference FASTA file not found: {args.ref_fasta}")
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

def run_maf2vcf(args, temp_vcf):
    """Run maf2vcf.py to convert MAF to VCF"""
    logger.info("Converting MAF to VCF...")
    
    # Build command
    script_dir = os.path.dirname(os.path.abspath(__file__))
    maf2vcf_path = os.path.join(script_dir, 'maf2vcf.py')
    
    cmd = [
        sys.executable,
        maf2vcf_path,
        '--input-maf', args.input_maf,
        '--output-vcf', temp_vcf
    ]
    
    # Add optional arguments
    if args.tumor_id:
        cmd.extend(['--tumor-id', args.tumor_id])
    if args.normal_id:
        cmd.extend(['--normal-id', args.normal_id])
    if args.vcf_tumor_id:
        cmd.extend(['--vcf-tumor-id', args.vcf_tumor_id])
    if args.vcf_normal_id:
        cmd.extend(['--vcf-normal-id', args.vcf_normal_id])
    if args.ref_fasta:
        cmd.extend(['--ref-fasta', args.ref_fasta])
    if args.tum_depth_col:
        cmd.extend(['--tum-depth-col', args.tum_depth_col])
    if args.tum_rad_col:
        cmd.extend(['--tum-rad-col', args.tum_rad_col])
    if args.tum_vad_col:
        cmd.extend(['--tum-vad-col', args.tum_vad_col])
    if args.nrm_depth_col:
        cmd.extend(['--nrm-depth-col', args.nrm_depth_col])
    if args.nrm_rad_col:
        cmd.extend(['--nrm-rad-col', args.nrm_rad_col])
    if args.nrm_vad_col:
        cmd.extend(['--nrm-vad-col', args.nrm_vad_col])
    if args.verbose:
        cmd.append('--verbose')
    
    # Run maf2vcf
    try:
        logger.info(f"Running command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running maf2vcf: {e}")
        sys.exit(1)

def run_vcf2maf(args, temp_vcf):
    """Run vcf2maf.py to convert VCF to MAF"""
    logger.info("Converting VCF to MAF...")
    
    # Build command
    script_dir = os.path.dirname(os.path.abspath(__file__))
    vcf2maf_path = os.path.join(script_dir, 'vcf2maf.py')
    
    cmd = [
        sys.executable,
        vcf2maf_path,
        '--input-vcf', temp_vcf,
        '--output-maf', args.output_maf
    ]
    
    # Add optional arguments
    if args.tumor_id:
        cmd.extend(['--tumor-id', args.tumor_id])
    if args.normal_id:
        cmd.extend(['--normal-id', args.normal_id])
    if args.vcf_tumor_id:
        cmd.extend(['--vcf-tumor-id', args.vcf_tumor_id])
    if args.vcf_normal_id:
        cmd.extend(['--vcf-normal-id', args.vcf_normal_id])
    if args.vep_path:
        cmd.extend(['--vep-path', args.vep_path])
    if args.vep_data:
        cmd.extend(['--vep-data', args.vep_data])
    if args.vep_forks:
        cmd.extend(['--vep-forks', str(args.vep_forks)])
    if args.ref_fasta:
        cmd.extend(['--ref-fasta', args.ref_fasta])
    if args.species:
        cmd.extend(['--species', args.species])
    if args.ncbi_build:
        cmd.extend(['--ncbi-build', args.ncbi_build])
    if args.cache_version:
        cmd.extend(['--cache-version', args.cache_version])
    if args.maf_center:
        cmd.extend(['--maf-center', args.maf_center])
    if args.retain_info:
        cmd.extend(['--retain-info', args.retain_info])
    if args.retain_fmt:
        cmd.extend(['--retain-fmt', args.retain_fmt])
    if args.min_hom_vaf:
        cmd.extend(['--min-hom-vaf', str(args.min_hom_vaf)])
    if args.inhibit_vep:
        cmd.append('--inhibit-vep')
    if args.verbose:
        cmd.append('--verbose')
    
    # Run vcf2maf
    try:
        logger.info(f"Running command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running vcf2maf: {e}")
        sys.exit(1)

def main():
    """Main function"""
    # Parse command line arguments
    args = parse_arguments()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Check prerequisites
    check_prerequisites(args)
    
    # Create a temporary file for the intermediate VCF
    with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as temp_vcf_file:
        temp_vcf = temp_vcf_file.name
    
    try:
        # Convert MAF to VCF
        run_maf2vcf(args, temp_vcf)
        
        # Convert VCF back to MAF
        run_vcf2maf(args, temp_vcf)
        
        logger.info(f"Reannotation complete. MAF file written to {args.output_maf}")
    
    finally:
        # Clean up temporary files
        if os.path.exists(temp_vcf):
            os.unlink(temp_vcf)

if __name__ == '__main__':
    main()
