# File: Snakefile

import os
from pathlib import Path

# Output directories
OUTPUT_DIR = "output"
BAM_DIR = os.path.join(OUTPUT_DIR, "bams")
PILEUP_DIR = os.path.join(OUTPUT_DIR, "pileups")
MAF_DIR = os.path.join(OUTPUT_DIR, "mafs")
VALIDATION_DIR = os.path.join(OUTPUT_DIR, "validation")
COMPARISON_DIR = os.path.join(OUTPUT_DIR, "comparison")

# Create output directories
for d in [OUTPUT_DIR, BAM_DIR, PILEUP_DIR, MAF_DIR, VALIDATION_DIR, COMPARISON_DIR]:
    Path(d).mkdir(parents=True, exist_ok=True)

# Enhanced Configuration
config = {
    "samples": ["tumor1"],
    "somatosim_data": {
        "input_bam": "SomatoSim/test_data/test_BAM.bam",
        "input_bed": "SomatoSim/test_data/test_BED.bed"
    },
    "coverage_levels": {
        "ultra_low": {
            "vaf_low": 0.01,  # More realistic VAF range for ultra-low coverage
            "vaf_high": 0.05,
            "depth_multiplier": 0.5  # 50% of original depth
        },
        "very_low": {
            "vaf_low": 0.05,
            "vaf_high": 0.15,
            "depth_multiplier": 0.75  # 75% of original depth
        },
        "low": {
            "vaf_low": 0.15,
            "vaf_high": 0.3,
            "depth_multiplier": 1  # 100% of original depth
        },
        "medium": {
            "vaf_low": 0.3,
            "vaf_high": 0.5,
            "depth_multiplier": 1.5  # 150% Original depth
        }
    },
    "mutation_profiles": {
        "tumor1": {
            "number_snv": 10000,  # Increased number of variants
            "random_seed": 42,
            "mutation_types": {
                "C>T": 0.5,  # Common in cancer
                "C>A": 0.2,
                "T>C": 0.15,
                "T>G": 0.15
            }
        }
    },
    "mapping_quality": {
        "coverage_MQ": 20,  # Increased mapping quality threshold
        "coverage_BQ": 20,  # Increased base quality threshold
        "read_min_MQ": 20,
        "position_min_BQ": 20
    },
    "confidence_intervals": {
        "alpha": 0.01,  # 99% confidence level
        "minimum_depth": 10  # Minimum read depth for reliable calls
    },
    "variant_callers": {
        "standard": {
            "min_depth": 10,
            "min_vaf": 0.05
        },
        "multinomial": {
            "alpha": 0.01,
            "min_depth": 10
        }
    }
}

# Global wildcard constraints
wildcard_constraints:
    sample="|".join(config["samples"]),
    coverage="|".join(config["coverage_levels"].keys())

rule all:
    input:
        # Basic outputs
        expand(
            os.path.join(BAM_DIR, "{sample}_{coverage}.bam"),
            sample=config["samples"],
            coverage=config["coverage_levels"].keys()
        ),
        expand(
            os.path.join(PILEUP_DIR, "{sample}_{coverage}.pileup"),
            sample=config["samples"],
            coverage=config["coverage_levels"].keys()
        ),
        expand(
            os.path.join(MAF_DIR, "{sample}_{coverage}.maf"),
            sample=config["samples"],
            coverage=config["coverage_levels"].keys()
        ),
        # Enhanced validation outputs
        os.path.join(VALIDATION_DIR, "multinomial_validation_report.txt"),
        os.path.join(VALIDATION_DIR, "validation_plots.pdf"),
        # Comparison outputs
        os.path.join(COMPARISON_DIR, "method_comparison.tsv"),
        os.path.join(COMPARISON_DIR, "error_profiles.tsv"),
        os.path.join(COMPARISON_DIR, "coverage_impact.tsv")

rule simulate_variants:
    input:
        input_bam=config["somatosim_data"]["input_bam"],
        input_bed=config["somatosim_data"]["input_bed"]
    output:
        bam=os.path.join(BAM_DIR, "{sample}_{coverage}.bam"),
        bai=os.path.join(BAM_DIR, "{sample}_{coverage}.bam.bai"),
        sim_output=os.path.join(BAM_DIR, "{sample}_{coverage}.simulation_output.txt"),
        sim_log=os.path.join(BAM_DIR, "{sample}_{coverage}.simulation_log.txt")
    params:
        output_dir=BAM_DIR,
        output_prefix=lambda wildcards: f"{wildcards.sample}_{wildcards.coverage}",
        vaf_low=lambda wildcards: config["coverage_levels"][wildcards.coverage]["vaf_low"],
        vaf_high=lambda wildcards: config["coverage_levels"][wildcards.coverage]["vaf_high"],
        depth_multiplier=lambda wildcards: config["coverage_levels"][wildcards.coverage]["depth_multiplier"],
        number_snv=lambda wildcards: config["mutation_profiles"][wildcards.sample]["number_snv"],
        random_seed=lambda wildcards: config["mutation_profiles"][wildcards.sample]["random_seed"],
        coverage_MQ=config["mapping_quality"]["coverage_MQ"],
        coverage_BQ=config["mapping_quality"]["coverage_BQ"],
        read_min_MQ=config["mapping_quality"]["read_min_MQ"],
        position_min_BQ=config["mapping_quality"]["position_min_BQ"]
    log:
        os.path.join("logs", "{sample}_{coverage}_somatosim.log")
    shell:
        """
        set -x  # Print commands for debugging
        
        # Create directories
        mkdir -p logs
        mkdir -p {params.output_dir}
        
        # Basic simulation without mutation types for debugging
        somatosim_wrapper.sh \
            -i {input.input_bam} \
            -b {input.input_bed} \
            -o {params.output_dir} \
            --output-prefix {params.output_prefix} \
            --vaf-low {params.vaf_low} \
            --vaf-high {params.vaf_high} \
            --number-snv {params.number_snv} \
            --random-seed {params.random_seed} \
            --coverage-MQ {params.coverage_MQ} \
            --coverage-BQ {params.coverage_BQ} \
            --read-min-MQ {params.read_min_MQ} \
            --position-min-BQ {params.position_min_BQ} \
            2>&1 | tee -a {log}
            
        # Print debug info
        echo "Debug: Checking for output files..." >> {log}
        ls -l {params.output_dir} >> {log}
        
        # Copy files if they exist
        if [ -f "{params.output_dir}/{params.output_prefix}.somatosim.bam" ]; then
            echo "Debug: Found BAM file, copying..." >> {log}
            cp -v "{params.output_dir}/{params.output_prefix}.somatosim.bam" "{output.bam}"
            cp -v "{params.output_dir}/{params.output_prefix}.somatosim.bam.bai" "{output.bai}"
        else
            echo "Error: Expected BAM file not found at {params.output_dir}/{params.output_prefix}.somatosim.bam" >> {log}
            echo "Directory contents:" >> {log}
            ls -l {params.output_dir} >> {log}
            exit 1
        fi
        
        if [ -f "{params.output_dir}/{params.output_prefix}_simulation_output.txt" ]; then
            echo "Debug: Found simulation output files, copying..." >> {log}
            cp -v "{params.output_dir}/{params.output_prefix}_simulation_output.txt" "{output.sim_output}"
            cp -v "{params.output_dir}/{params.output_prefix}_simulation_log.txt" "{output.sim_log}"
        else
            echo "Error: Expected output file not found at {params.output_dir}/{params.output_prefix}_simulation_output.txt" >> {log}
            echo "Directory contents:" >> {log}
            ls -l {params.output_dir} >> {log}
            exit 1
        fi
        """

rule bam_to_pileup:
    input:
        bam=os.path.join(BAM_DIR, "{sample}_{coverage}.bam")
    output:
        pileup=os.path.join(PILEUP_DIR, "{sample}_{coverage}.pileup")
    params:
        coverage_MQ=config["mapping_quality"]["coverage_MQ"],
        coverage_BQ=config["mapping_quality"]["coverage_BQ"]
    log:
        os.path.join("logs", "{sample}_{coverage}_pileup.log")
    shell:
        """
        mkdir -p logs
        samtools_wrapper.sh mpileup \
            -q {params.coverage_MQ} \
            -Q {params.coverage_BQ} \
            {input.bam} > {output.pileup} 2> {log}
        """

rule sim_output_to_maf:
    input:
        sim_output=os.path.join(BAM_DIR, "{sample}_{coverage}.simulation_output.txt")
    output:
        maf=os.path.join(MAF_DIR, "{sample}_{coverage}.maf")
    log:
        os.path.join("logs", "{sample}_{coverage}_maf.log")
    script:
        "scripts/convert_truth_to_maf.py"


rule validate_multinomial:
    input:
        pileups=expand(
            os.path.join(PILEUP_DIR, "{sample}_{coverage}.pileup"),
            sample=config["samples"],
            coverage=config["coverage_levels"].keys()
        ),
        sim_outputs=expand(
            os.path.join(BAM_DIR, "{sample}_{coverage}.simulation_output.txt"),
            sample=config["samples"],
            coverage=config["coverage_levels"].keys()
        ),
        mafs=expand(
            os.path.join(MAF_DIR, "{sample}_{coverage}.maf"),
            sample=config["samples"],
            coverage=config["coverage_levels"].keys()
        )
    output:
        report=os.path.join(VALIDATION_DIR, "multinomial_validation_report.txt"),
        plots=os.path.join(VALIDATION_DIR, "validation_plots.pdf"),
        method_comparison=os.path.join(COMPARISON_DIR, "method_comparison.tsv"),
        error_profiles=os.path.join(COMPARISON_DIR, "error_profiles.tsv"),
        coverage_impact=os.path.join(COMPARISON_DIR, "coverage_impact.tsv")
    params:
        standard_params=config["variant_callers"]["standard"],
        multinomial_params=config["variant_callers"]["multinomial"],
        joined_pileups=" ".join(expand(
            os.path.join(PILEUP_DIR, "{sample}_{coverage}.pileup"),
            sample=config["samples"],
            coverage=config["coverage_levels"].keys()
        )),
        joined_sim_outputs=" ".join(expand(
            os.path.join(BAM_DIR, "{sample}_{coverage}.simulation_output.txt"),
            sample=config["samples"],
            coverage=config["coverage_levels"].keys()
        )),
        joined_mafs=" ".join(expand(
            os.path.join(MAF_DIR, "{sample}_{coverage}.maf"),
            sample=config["samples"],
            coverage=config["coverage_levels"].keys()
        )),
        coverage_levels=" ".join(config["coverage_levels"].keys())
    resources:
        mem_mb=4000
    threads: 3
    log:
        os.path.join("logs", "validate_multinomial.log")
    shell:
        """
        somatosim_wrapper.sh scripts/validate_multinomial.py \
            --pileups {params.joined_pileups} \
            --sim-outputs {params.joined_sim_outputs} \
            --mafs {params.joined_mafs} \
            --coverage-levels {params.coverage_levels} \
            --report {output.report} \
            --plots {output.plots} \
            --method-comparison {output.method_comparison} \
            --error-profiles {output.error_profiles} \
            --coverage-impact {output.coverage_impact} \
            --standard-min-depth {params.standard_params[min_depth]} \
            --standard-min-vaf {params.standard_params[min_vaf]} \
            --multinomial-alpha {params.multinomial_params[alpha]} \
            --multinomial-min-depth {params.multinomial_params[min_depth]} 2> {log}
        """
