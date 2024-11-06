# Cancer Variant Simulation and Validation Pipeline

A comprehensive framework for validating variant calling methods in cancer genomics, focusing on addressing the challenges of low-read-count false negatives in next-generation sequencing (NGS) data.

## Overview

This pipeline validates a multinomial-based variant calling approach that improves the detection of true cancer-related mutations, particularly in low-coverage regions. It was designed to simulate and validate variants across different coverage depths and variant allele frequencies (VAFs).

### Key Features
- Simulation of cancer variants using SomatoSim
- Configurable coverage levels from ultra-low to medium depth
- Realistic cancer mutation profiles (e.g., C>T transitions common in cancer)
- Quality filtering based on mapping and base qualities
- Statistical validation using multinomial confidence intervals
- Comprehensive performance metrics and visualizations

### Pipeline Components
- BAM file generation with controlled variant allele frequencies
- Pileup generation and analysis
- MAF file conversion for compatibility with standard tools
- Comprehensive validation and comparison of variant calling methods
- Statistical analysis and visualization of results

## Installation

### Prerequisites
```
- Conda or Miniconda
- Git
- Python 3.7 or higher
```

### Setup Steps

1. Clone the pipeline repository:
```bash
git clone [repository-url]
cd [repository-name]
```

2. Set up conda channels:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. Download and install SomatoSim:

*Make sure that you install SomatoSim inside this cloned repository the directory structure is shown below

```bash
git clone https://github.com/BieseckerLab/SomatoSim.git
cd SomatoSim
```

4. Create and configure environments:
```bash
chmod +x make_separate_envs.sh
./make_separate_envs.sh
```

The setup script creates:
- Four conda environments for different components
- Wrapper scripts in ~/bin directory
- Required environment variables

### System Structure

#### Pipeline Directory
```
.
├── SomatoSim/               # SomatoSim installation
├── scripts/                 # Pipeline scripts
│   ├── validate_multinomial.py
│   └── convert_truth_to_maf.py
├── output/                  # Created when pipeline runs
├── Snakefile               # Pipeline workflow
└── make_separate_envs.sh   # Environment setup script
```

#### Wrapper Scripts Location
```
~/bin/                      # Wrapper scripts directory
   ├── somatosim_wrapper.sh
   ├── samtools_wrapper.sh
   └── snakemake_wrapper.sh
```

### Environment Setup

The pipeline uses four separate conda environments:

1. **bam_sim_base**: Core Python packages
   ```yaml
   python=3.9
   numpy
   pandas
   scipy
   scikit-learn
   matplotlib=3.9.2
   seaborn
   pytest
   pysam=0.19.1
   biopython
   ```

2. **bam_sim_somatosim**: SomatoSim environment
   ```yaml
   python=3.7
   matplotlib=3.2.0
   numpy=1.18.5
   pandas=1.0.2
   pysam=0.15.4
   samtools
   scikit-learn=1.0.2
   ```

3. **bam_sim_samtools**: Samtools operations
   ```yaml
   samtools=1.9
   htslib=1.9
   ```

4. **bam_sim_snake**: Snakemake workflow
   ```yaml
   python=3.9
   snakemake-minimal=7.32.4
   ```

### Installation Verification

1. Check your PATH configuration:
```bash
echo $PATH | grep -q "$HOME/bin" && echo "Path is set correctly" || echo "Add ~/bin to PATH"
```

2. If needed, add to your ~/.bashrc or ~/.zshrc:
```bash
export PATH="$HOME/bin:$PATH"
```

3. Verify tool installations:
```bash
somatosim_wrapper.sh --help
samtools_wrapper.sh --version
snakemake_wrapper.sh --version
```

## Usage

### Configuration

Edit the config section in Snakefile to set your parameters:
you can also change the number of SNVs in the Snakefile easily
at "number_snv": 10000
```python
config = {
    "samples": ["tumor1"],
    "coverage_levels": {
        "ultra_low": {
            "vaf_low": 0.01,
            "vaf_high": 0.05,
            "depth_multiplier": 0.5
        },
        "very_low": {
            "vaf_low": 0.05,
            "vaf_high": 0.15,
            "depth_multiplier": 0.75
        },
        "low": {
            "vaf_low": 0.15,
            "vaf_high": 0.3,
            "depth_multiplier": 1
        },
        "medium": {
            "vaf_low": 0.3,
            "vaf_high": 0.5,
            "depth_multiplier": 1.5
        }
    }
}
```

### Running the Pipeline

1. Start the pipeline:
```bash
snakemake_wrapper.sh --cores 4 -p
```

2. Monitor progress in the log files:
```bash
tail -f logs/validate_multinomial.log
```

### Output Structure

The pipeline generates several output directories:
```
output/
├── bams/         # Simulated BAM files
├── pileups/      # Pileup files
├── mafs/         # MAF format files
├── validation/   # Validation reports
└── comparison/   # Method comparisons
```

### Key Outputs

1. **Validation Reports**
   - `validation/multinomial_validation_report.txt`
   - `validation/validation_plots.pdf`

2. **Comparison Results**
   - `comparison/method_comparison.tsv`
   - `comparison/error_profiles.tsv`
   - `comparison/coverage_impact.tsv`

## Troubleshooting

### Common Issues

1. **SomatoSim Installation**
   - Use Python 3.7 environment
   - Try: `pip install -e .`

2. **Wrapper Scripts**
   - Ensure ~/bin is in PATH
   - Make scripts executable: `chmod +x ~/bin/*.sh`
   - Source config: `source ~/.bashrc`

3. **Conda Environments**
   - Update conda: `conda update -n base conda`
   - Check channel order
   - Reinstall specific environment: `conda env create -f environment.yml`

## Citation

If you use this pipeline in your research, please cite:
[Citation information to be added]

## License

MIT
