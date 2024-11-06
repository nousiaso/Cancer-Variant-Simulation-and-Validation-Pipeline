#!/usr/bin/env python
# File: scripts/convert_truth_to_maf.py

import os
import pandas as pd
import numpy as np
from collections import OrderedDict

REQUIRED_COLUMNS = [
    'Hugo_Symbol',
    'Entrez_Gene_Id',
    'Center',
    'NCBI_Build',
    'Chromosome',
    'Start_Position',
    'End_Position',
    'Strand',
    'Variant_Classification',
    'Variant_Type',
    'Reference_Allele',
    'Tumor_Seq_Allele1',
    'Tumor_Seq_Allele2',
    'Tumor_Sample_Barcode',
    'Matched_Norm_Sample_Barcode'
]

def convert_truth_to_maf(truth_file, output_maf, ref_gene_file=None):
    """Convert truth files to MAF format with proper specifications."""
    
    # Read simulation output file
    df = pd.read_csv(truth_file, sep='\s+')  # Using \s+ to handle variable whitespace
    
    # Extract sample name from filename
    sample_name = os.path.basename(truth_file).replace('.simulation_output.txt', '')
    
    # Calculate read counts based on VAF and coverage
    df['alt_count'] = (df['input_VAF'] * df['input_coverage']).round().astype(int)
    df['ref_count'] = df['input_coverage'] - df['alt_count']
    
    # Create MAF data with all required fields
    maf_data = OrderedDict({
        'Hugo_Symbol': 'Unknown',
        'Entrez_Gene_Id': 0,
        'Center': 'Simulation_Center',
        'NCBI_Build': 'GRCh38',
        'Chromosome': df['chromosome'].astype(str),
        'Start_Position': df.iloc[:, 1],  # First position column
        'End_Position': df.iloc[:, 2],    # Second position column
        'Strand': '+',
        'Variant_Classification': 'SNP',
        'Variant_Type': 'SNP',
        'Reference_Allele': df['ref_allele'].str.upper(),
        'Tumor_Seq_Allele1': df['ref_allele'].str.upper(),
        'Tumor_Seq_Allele2': df['alt_allele'].str.upper(),
        'Tumor_Sample_Barcode': sample_name,
        'Matched_Norm_Sample_Barcode': 'normal',
        't_depth': df['input_coverage'],
        't_ref_count': df['ref_count'],
        't_alt_count': df['alt_count'],
        'n_depth': df['input_coverage'],
        'n_ref_count': df['input_coverage'],
        'n_alt_count': 0,
        'Validation_Status': 'Unknown',
        'Mutation_Status': 'Somatic',
        'HGVSp_Short': '',
        'Transcript_ID': '',
        'Exon_Number': '',
        'dbSNP_RS': '',
        'COSMIC_ID': '',
        't_var_freq': df['input_VAF'],
        'n_var_freq': 0.0,
        'IMPACT': 'MODERATE',
        'FILTER': 'PASS',
        'actual_vaf': df['output_VAF']  # Store the actual VAF for validation
    })
    
    # Convert to DataFrame and ensure proper types
    maf_df = pd.DataFrame(maf_data)
    
    # Ensure integer types for relevant columns
    int_columns = ['Start_Position', 'End_Position', 'Entrez_Gene_Id',
                  't_depth', 't_ref_count', 't_alt_count',
                  'n_depth', 'n_ref_count', 'n_alt_count']
    
    for col in int_columns:
        maf_df[col] = maf_df[col].astype(int)
    
    # Add variant_id column
    maf_df['variant_id'] = maf_df.apply(
        lambda x: f"{x['Chromosome']}_{x['Start_Position']}_{x['Reference_Allele']}_{x['Tumor_Seq_Allele2']}", 
        axis=1
    )
    
    # Categorize chromosomes
    chrom_categories = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                      '11', '12', '13', '14', '15', '16', '17', '18', '19', 
                      '20', '21', '22', 'X', 'Y', 'MT']
    maf_df['Chromosome'] = pd.Categorical(
        maf_df['Chromosome'],
        categories=chrom_categories,
        ordered=True
    )
    
    # Sort by chromosome and position
    maf_df = maf_df.sort_values(['Chromosome', 'Start_Position'])
    
    # Ensure all required columns are present
    for col in REQUIRED_COLUMNS:
        if col not in maf_df.columns:
            maf_df[col] = ''
    
    # Write to MAF file with correct formatting
    maf_df.to_csv(output_maf, sep='\t', index=False)
    
    return maf_df

def main():
    # Get input and output files from snakemake
    sim_output = snakemake.input.sim_output
    output_maf = snakemake.output.maf
    
    # Convert file
    convert_truth_to_maf(sim_output, output_maf)

if __name__ == '__main__':
    main()
