#!/usr/bin/env python3
import os
import pandas as pd
import gzip
import logging
import requests
from typing import Dict, Optional
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SimpleGeneAnnotator:
    def __init__(self, gtf_path: str = None):
        """Initialize with GTF file"""
        if gtf_path is None:
            gtf_path = "Homo_sapiens.GRCh38.108.gtf.gz"
        self.genes = self._load_genes(gtf_path)
        
    def _load_genes(self, gtf_path: str) -> pd.DataFrame:
        """Load genes from GTF file"""
        logger.info(f"Loading genes from {gtf_path}")
        
        genes = []
        open_fn = gzip.open if gtf_path.endswith('.gz') else open
        
        with open_fn(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                    
                # Parse attributes
                attrs = dict(x.strip().split(' ', 1) for x in 
                           fields[8].rstrip(';').split('; '))
                attrs = {k: v.strip('"') for k, v in attrs.items()}
                
                genes.append({
                    'chrom': fields[0].replace('chr', ''),
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'strand': fields[6],
                    'gene_name': attrs.get('gene_name', 'Unknown'),
                    'gene_id': attrs.get('gene_id', '0'),
                    'gene_type': attrs.get('gene_biotype', 'unknown')
                })
        
        return pd.DataFrame(genes)
    
    def get_gene_info(self, chrom: str, position: int) -> Dict[str, str]:
        """Get gene information for a given position"""
        try:
            # Remove 'chr' prefix if present
            chrom = str(chrom).replace('chr', '')
            position = int(position)
            
            # Find overlapping genes
            matching_genes = self.genes[
                (self.genes['chrom'] == chrom) &
                (self.genes['start'] <= position) &
                (self.genes['end'] >= position)
            ]
            
            if not matching_genes.empty:
                gene = matching_genes.iloc[0]
                return {
                    'Hugo_Symbol': gene['gene_name'],
                    'Entrez_Gene_Id': gene['gene_id'],
                    'Gene_Strand': gene['strand'],
                    'Variant_Classification': self._get_variant_classification(gene['gene_type']),
                    'Gene_Type': gene['gene_type'],
                    'Impact': self._get_impact(gene['gene_type'])
                }
            
            return self._get_igr_annotation()
            
        except Exception as e:
            logger.warning(f"Gene annotation failed for {chrom}:{position} - {e}")
            return self._get_igr_annotation()
    
    def _get_variant_classification(self, gene_type: str) -> str:
        """Determine variant classification based on gene type"""
        if gene_type == 'protein_coding':
            return 'Missense_Mutation'
        elif gene_type == 'pseudogene':
            return 'RNA'
        elif gene_type.endswith('RNA'):
            return 'RNA'
        elif gene_type == 'lincRNA':
            return 'LINC'
        return 'IGR'
    
    def _get_impact(self, gene_type: str) -> str:
        """Determine impact based on gene type"""
        if gene_type == 'protein_coding':
            return 'MODERATE'
        elif gene_type in ['pseudogene', 'lincRNA']:
            return 'LOW'
        return 'MODIFIER'
    
    def _get_igr_annotation(self) -> Dict[str, str]:
        """Return standard IGR (intergenic region) annotation"""
        return {
            'Hugo_Symbol': 'IGR',
            'Entrez_Gene_Id': '0',
            'Gene_Strand': '+',
            'Variant_Classification': 'IGR',
            'Gene_Type': 'intergenic_region',
            'Impact': 'MODIFIER'
        }

def download_gtf(output_dir: str = ".") -> str:
    """Download the GTF file if not present"""
    url = "https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz"
    output_path = os.path.join(output_dir, "Homo_sapiens.GRCh38.108.gtf.gz")
    
    if not os.path.exists(output_path):
        logger.info(f"Downloading GTF file from {url}")
        # Create directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Download with progress indication
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        
        with open(output_path, 'wb') as f:
            if total_size == 0:
                f.write(response.content)
            else:
                downloaded = 0
                total_size_mb = total_size / (1024 * 1024)
                for data in response.iter_content(chunk_size=4096):
                    downloaded += len(data)
                    f.write(data)
                    done = int(50 * downloaded / total_size)
                    downloaded_mb = downloaded / (1024 * 1024)
                    print(f'\rDownloading GTF: [{"=" * done}{" " * (50-done)}] '
                          f'{downloaded_mb:.1f}/{total_size_mb:.1f} MB', 
                          end='', flush=True)
        print()  # New line after progress bar
        
        logger.info(f"Downloaded GTF file to {output_path}")
    else:
        logger.info(f"GTF file already exists at {output_path}")
    
    return output_path

def main():
    """Example usage"""
    # Create a data directory
    data_dir = "annotation_data"
    
    # Download GTF if needed
    gtf_path = download_gtf(data_dir)
    
    # Initialize annotator
    annotator = SimpleGeneAnnotator(gtf_path)
    
    # Test some known positions
    test_positions = [
        ('17', 7675000),  # TP53
        ('13', 32315000),  # BRCA2
        ('7', 140753336)  # BRAF
    ]
    
    for chrom, pos in test_positions:
        info = annotator.get_gene_info(chrom, pos)
        print(f"\nAnnotation for {chrom}:{pos}")
        for key, value in info.items():
            print(f"{key}: {value}")

if __name__ == "__main__":
    main()
