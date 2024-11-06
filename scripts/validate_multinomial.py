import os
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, List, Counter as CounterType
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
import argparse
import scipy.stats as st  # Import for scipy.stats used in MultinomialValidator

class ReadProcessor:
    """Process read data from pileup format."""
    
    def __init__(self, bases: str, quals: str, mq: str, ref: str, min_depth: int = 10):
        self.bases = bases
        self.quals = quals
        self.mq = mq
        self.ref = ref.upper()
        self.min_depth = min_depth
        self.filtered_counter = self._process_bases()
        self.call, self.genotype = self._make_call()
        self.vaf = self._calculate_vaf()

    def _process_bases(self) -> CounterType:
        """Process bases string and return filtered base counts."""
        base_counter = Counter()
        for base in self.bases.upper():
            if base in 'ACGT':
                base_counter[base] += 1
        return base_counter

    def _make_call(self) -> Tuple[str, Optional[str]]:
        """Make variant call based on processed bases."""
        if sum(self.filtered_counter.values()) < self.min_depth:
            return 'MISSING', None
            
        # Sort bases by count
        ordered_bases = sorted(self.filtered_counter.items(), 
                             key=lambda x: x[1], reverse=True)
        if not ordered_bases:
            return 'MISSING', None
            
        # Get ordered bases and counts
        ord_bases, ord_counts = zip(*ordered_bases)
        
        # Get multinomial validator result
        validator = MultinomialValidator()
        validation_result = validator.validate_call(self.filtered_counter, 
                                                  sum(self.filtered_counter.values()))
        
        if validation_result == "MISSING":
            return 'MISSING', None
        elif validation_result == "VALID":
            top_base = ord_bases[0]
            if top_base == self.ref:
                return 'REF', None
            else:
                return 'ALT', top_base
        else:  # CHECK
            return 'CHECK', '?'

    def _calculate_vaf(self) -> float:
        """Calculate variant allele fraction."""
        total = sum(self.filtered_counter.values())
        if total == 0:
            return 0.0
            
        if self.genotype and self.genotype in self.filtered_counter:
            return self.filtered_counter[self.genotype] / total
        return 0.0

    def get_mutation_type(self) -> str:
        """Determine if SNV."""
        if not self.genotype or self.call != 'ALT':
            return 'NONE'
        if len(self.genotype) == 1:  # Single base change
            return 'SNV'
        return 'NONE'

class MultinomialValidator:
    """Enhanced validator for variant calls using multinomial statistics."""
    
    def __init__(self, alpha: float = 0.01, min_depth: int = 10):
        """Initialize validator with specified parameters."""
        self.alpha = alpha
        self.min_depth = min_depth
        self.results = defaultdict(dict)
        self.mutation_types = ['SNV']  # Only SNVs
        logging.debug(f"Initialized MultinomialValidator with alpha={alpha}, min_depth={min_depth}")

    def validate_call(self, base_counts, total_depth):
        """Validation matching the main pipeline's ReadProcessor logic"""
        # Sort bases by count to match process_pileup.py's ordering
        ord_bases, ord_counts = zip(*sorted(base_counts.items(), 
                                          key=lambda x: x[1], 
                                          reverse=True))
        
        # Match the main pipeline's CI calculation
        counts3 = (ord_counts[:2] + (sum(ord_counts[2:]),))
        CI_tuples = self._get_multinom_CIs(counts3)
        CI_mins = [int(np.floor(i[0])) for i in CI_tuples]
        
        # Exactly match process_pileup.py's genotype logic
        if CI_mins[0] < 1:
            return "MISSING"  # A touches zero
        elif CI_mins[1] < 1:
            return "VALID"    # Only B touches zero
        elif CI_mins[2] < 1:
            return "VALID"    # C+D touches zero
        else:
            return "CHECK"    # Nothing touches zero

    def _get_multinom_CIs(self, vals):
        """Exact copy of the CI calculation from process_pileup.py"""
        CI_list = list()
        N = sum(vals)
        for val in vals:
            if N > 0:
                z_a2 = st.norm.ppf(1 - self.alpha / 2)
                p = np.float64(val) / N
                this_CI_lt = p + z_a2 ** 2 / (2 * N)
                this_CI_rt = z_a2 * np.sqrt((p * (1 - p) + z_a2 ** 2 / 4 / N) / N)
                this_CI_b = 1 + z_a2 ** 2 / N
                this_CI_max = (this_CI_lt + this_CI_rt) / this_CI_b
                this_CI_min = (this_CI_lt - this_CI_rt) / this_CI_b
                CI_tuple = (this_CI_min * N, this_CI_max * N)
            else:
                CI_tuple = (0, 0)
            CI_list.append(CI_tuple)
        return tuple(CI_list)

    def _evaluate_call(self, call: str, genotype: str, true_alt: Optional[str]) -> bool:
        """Evaluate if a variant call matches the truth data."""
        if true_alt is None:
            return call == 'REF' or call == 'MISSING'
        if call == 'MISSING':
            return False  # This is where a false negative can occur
        if call == 'REF':
            return False
        return genotype == true_alt

    def analyze_all_files(self, pileup_files: Dict[str, str], 
                         truth_files: Dict[str, str]) -> None:
        """Analyze multiple pileup files with corresponding truth files."""
        for coverage, pileup_file in pileup_files.items():
            truth_file = truth_files.get(coverage)
            self.results[coverage] = self.analyze_pileup(pileup_file, truth_file)

    def analyze_pileup(self, pileup_file: str, truth_file: Optional[str] = None) -> Dict:
        """Analyze pileup file using both standard and multinomial methods."""
        if not self._validate_input_files(pileup_file):
            return self._get_empty_results()
            
        true_variants = self._load_truth_data(truth_file)
        
        metrics = defaultdict(list)
        standard_metrics = defaultdict(list)
        positions_processed = 0
        mutation_counts = Counter()
        depth_accuracy = defaultdict(list)
        vaf_accuracy = defaultdict(list)
        
        with open(pileup_file) as f:
            for line_count, line in enumerate(f, 1):
                if line_count % 1000 == 0:
                    logging.debug(f"Processed {line_count} positions")
                
                try:
                    result = self._process_pileup_line(line, true_variants)
                    if result is None:
                        continue
                    
                    # Track metrics for both methods
                    self._update_metrics(metrics, result, 'multinomial')
                    self._update_metrics(standard_metrics, result, 'standard')
                    
                    # Track additional analysis metrics
                    self._track_mutation_metrics(result, mutation_counts)
                    self._track_depth_metrics(result, depth_accuracy)
                    self._track_vaf_metrics(result, vaf_accuracy)
                    
                    positions_processed += 1
                        
                except Exception as e:
                    logging.error(f"Error processing line {line_count}: {str(e)}")
                    continue
        
        if positions_processed == 0:
            return self._get_empty_results()
        
        summary = self._calculate_summary_statistics(
            metrics, mutation_counts, depth_accuracy,
            vaf_accuracy, positions_processed
        )
        
        standard_summary = self._calculate_summary_statistics(
            standard_metrics, mutation_counts, depth_accuracy,
            vaf_accuracy, positions_processed, method='standard'
        )
            
        return {
            'metrics': metrics,
            'summary': summary,
            'standard_metrics': standard_metrics,
            'standard_summary': standard_summary,
            'mutation_profile': dict(mutation_counts),
            'depth_accuracy': dict(depth_accuracy),
            'vaf_accuracy': dict(vaf_accuracy)
        }

    def _process_pileup_line(self, line: str, true_variants: Dict) -> Optional[Dict]:
        """Process single pileup line and return results."""
        fields = line.strip().split('\t')
        if len(fields) < 6:
            logging.warning(f"Invalid line format: {line}")
            return None
            
        chrom, pos, ref, depth, bases, quals = fields[:6]
        position = f"{chrom}:{pos}"
        
        # Get mapping quality if available
        mq = fields[7] if len(fields) > 7 else quals
        
        # Process reads
        processor = ReadProcessor(bases, quals, mq, ref, self.min_depth)
        standard_call = self._make_standard_call(processor.filtered_counter, ref)
        
        # Get truth data
        true_var = true_variants.get(position, {})
        true_alt = true_var.get('alt_allele')
        true_vaf = true_var.get('vaf')
        
        # Evaluate calls
        multi_correct = self._evaluate_call(processor.call, processor.genotype, true_alt)
        standard_correct = self._evaluate_call(standard_call[0], standard_call[1], true_alt)
        
        return {
            'position': position,
            'total_depth': sum(processor.filtered_counter.values()),
            'base_counts': processor.filtered_counter,
            'multinomial_call': processor.call,
            'standard_call': standard_call[0],
            'genotype': processor.genotype,
            'vaf': processor.vaf,
            'true_vaf': true_vaf,
            'ref_base': ref.upper(),
            'is_true_variant': true_alt is not None,
            'multinomial_correct': multi_correct,
            'standard_correct': standard_correct,
            'mutation_type': processor.get_mutation_type()
        }

    def _make_standard_call(self, base_counts, ref_base):
        """Make a standard variant call based on base counts and a reference base."""
        sorted_bases = sorted(base_counts.items(), key=lambda x: x[1], reverse=True)
        top_base, top_count = sorted_bases[0] if sorted_bases else (ref_base, 0)
        
        # Calculate VAF and total depth
        total = sum(base_counts.values())
        vaf = top_count / total if total > 0 else 0
        
        # Depth-based filtering
        if total < self.min_depth:
            return "MISSING", None
        
        # VAF threshold-based call
        if vaf < 0.2:
            return "MISSING", None
        elif top_base == ref_base:
            return "REF", None
        else:
            return "ALT", top_base

    def _update_metrics(self, metrics: Dict, result: Dict, method: str = 'multinomial'):
        """Update metrics dictionary with results from processed line."""
        key_prefix = '' if method == 'multinomial' else 'standard_'
        call_key = f"{key_prefix}call"
        correct_key = f"{key_prefix}correct"
    
        # Existing metrics
        metrics['total_depth'].append(result['total_depth'])
        metrics['vaf'].append(result['vaf'])
        metrics['true_vaf'].append(result.get('true_vaf'))
        metrics['is_true_variant'].append(result['is_true_variant'])
        metrics[call_key].append(result[f"{method}_call"])
        metrics[correct_key].append(result[f"{method}_correct"])
        
        # Add false negative tracking
        is_false_negative = (result['is_true_variant'] and 
                            result[f"{method}_call"] == 'MISSING')
        metrics.setdefault('false_negatives', []).append(is_false_negative)
        
        # Track depth for false negatives
        if is_false_negative:
            depth_bin = f"{result['total_depth'] // 10 * 10}-{(result['total_depth'] // 10 + 1) * 10}"
            metrics.setdefault('fn_by_depth', {}).setdefault(depth_bin, 0)
            metrics['fn_by_depth'][depth_bin] += 1

    def _track_mutation_metrics(self, result: Dict, mutation_counts: Counter):
        """Track mutation-specific metrics."""
        if result['is_true_variant']:
            mutation_type = result['mutation_type']
            if mutation_type in self.mutation_types:
                mutation_counts[mutation_type] += 1

    def _track_depth_metrics(self, result: Dict, depth_accuracy: Dict):
        """Track depth-related accuracy metrics."""
        depth = result['total_depth']
        depth_bin = f"{depth // 10 * 10}-{(depth // 10 + 1) * 10}"
        depth_accuracy[depth_bin].append(result['multinomial_correct'])

    def _track_vaf_metrics(self, result: Dict, vaf_accuracy: Dict):
        """Track VAF-related accuracy metrics."""
        if result['true_vaf'] is not None:
            vaf = result['true_vaf']
            vaf_bin = f"{vaf:.1f}-{(vaf + 0.1):.1f}"
            vaf_accuracy[vaf_bin].append(result['multinomial_correct'])

    def _calculate_summary_statistics(self, metrics: Dict, mutation_counts: Counter,
                                    depth_accuracy: Dict, vaf_accuracy: Dict,
                                    positions_processed: int, method: str = 'multinomial') -> Dict:
        """Calculate summary statistics from collected metrics."""
        key_prefix = '' if method == 'multinomial' else 'standard_'
        correct_key = f"{key_prefix}correct"
    
        total_depth = metrics.get('total_depth', [])
        correct_calls = sum(metrics.get(correct_key, []))
        
        # Calculate false negative rate
        true_variants = sum(metrics['is_true_variant'])
        false_negatives = sum(metrics.get('false_negatives', []))
        fn_rate = false_negatives / true_variants if true_variants > 0 else 0
        
        return {
            'total_positions': positions_processed,
            'mean_depth': np.mean(total_depth) if total_depth else 0,
            'median_depth': np.median(total_depth) if total_depth else 0,
            'calls': Counter(metrics.get(f"{key_prefix}call", [])),
            'accuracy': correct_calls / positions_processed if positions_processed > 0 else 0,
            'positions_evaluated': positions_processed,
            'mutation_counts': dict(mutation_counts),
            'depth_accuracy': {k: np.mean(v) for k, v in depth_accuracy.items()},
            'vaf_accuracy': {k: np.mean(v) for k, v in vaf_accuracy.items()},
            'false_negative_rate': fn_rate,
            'false_negatives': false_negatives,
            'fn_by_depth': metrics.get('fn_by_depth', {})
        }

    def _write_method_comparison(self, output_prefix: str):
        """Write method comparison data to file."""
        comparison_file = os.path.join(output_prefix, 'method_comparison.tsv')
        data = []
        for coverage in self.results:
            data.append({
                'Coverage': coverage,
                'Multinomial_Accuracy': self.results[coverage]['summary']['accuracy'],
                'Standard_Accuracy': self.results[coverage]['standard_summary']['accuracy'],
                'Mean_Depth': self.results[coverage]['summary']['mean_depth']
            })
    
        df = pd.DataFrame(data)
        df.to_csv(comparison_file, sep='\t', index=False)

    def _write_error_profiles(self, output_prefix: str):
        """Write error profile data to file."""
        error_file = os.path.join(output_prefix, 'error_profiles.tsv')
        data = []
        for coverage, result in self.results.items():
            for mut_type, count in result['summary']['mutation_counts'].items():
                data.append({
                    'Coverage': coverage,
                    'Mutation_Type': mut_type,
                    'Error_Count': count
                })
    
        df = pd.DataFrame(data)
        df.to_csv(error_file, sep='\t', index=False)

    def _write_coverage_impact(self, output_prefix: str):
        """Write coverage impact data to file."""
        impact_file = os.path.join(output_prefix, 'coverage_impact.tsv')
        data = []
        for coverage, result in self.results.items():
            data.append({
                'Coverage': coverage,
                'Mean_Depth': result['summary']['mean_depth'],
                'Accuracy': result['summary']['accuracy'],
                'Total_Positions': result['summary']['total_positions']
            })
    
        df = pd.DataFrame(data)
        df.to_csv(impact_file, sep='\t', index=False)

    def generate_report(self, report_file: str, plot_file: str):
        """Generate comprehensive validation report and plots."""
        logging.info(f"Generating report: {report_file}")
        os.makedirs(os.path.dirname(report_file), exist_ok=True)
        
        self._write_text_report(report_file)
        self._generate_plots(plot_file)
        
        logging.info(f"Report generation completed")

    def write_comparison_outputs(self, output_prefix: str):
        """Write comprehensive comparison analysis outputs."""
        os.makedirs(output_prefix, exist_ok=True)
        
        # Method comparison
        self._write_method_comparison(output_prefix)
        
        # Error profiles
        self._write_error_profiles(output_prefix)
        
        # Coverage impact
        self._write_coverage_impact(output_prefix)

    def _plot_method_comparison(self, ax):
        """Plot comparison between multinomial and standard methods."""
        coverages = list(self.results.keys())
        multi_acc = [res['summary'].get('accuracy', 0) for res in self.results.values()]
        std_acc = [res['standard_summary'].get('accuracy', 0) for res in self.results.values()]
        
        x = np.arange(len(coverages))
        width = 0.35
        
        ax.bar(x - width/2, multi_acc, width, label='Multinomial')
        ax.bar(x + width/2, std_acc, width, label='Standard')
        
        ax.set_ylabel('Accuracy')
        ax.set_title('Method Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels(coverages)
        ax.legend()
        plt.setp(ax.get_xticklabels(), rotation=45)

    def _plot_coverage_impact(self, ax):
        """Plot impact of coverage on accuracy."""
        coverages = list(self.results.keys())
        accuracies = [res['summary'].get('accuracy', 0) for res in self.results.values()]
        depths = [res['summary'].get('mean_depth', 0) for res in self.results.values()]
        
        color = 'tab:blue'
        ax.set_xlabel('Mean Coverage Depth')
        ax.set_ylabel('Accuracy', color=color)
        ax.plot(depths, accuracies, color=color, marker='o')
        ax.tick_params(axis='y', labelcolor=color)
        
        ax.set_title('Coverage Impact on Accuracy')

    def _plot_error_profiles(self, ax):
        """Plot error profiles across different coverage levels."""
        data = []
        for coverage, result in self.results.items():
            for mut_type, count in result['summary'].get('mutation_counts', {}).items():
                data.append({
                    'Coverage': coverage,
                    'Mutation Type': mut_type,
                    'Count': count
                })
        
        if data:
            df = pd.DataFrame(data)
            sns.barplot(data=df, x='Coverage', y='Count', hue='Mutation Type', ax=ax)
            ax.set_title('Error Profiles by Coverage')
            plt.setp(ax.get_xticklabels(), rotation=45)

    def _plot_vaf_analysis(self, ax):
        """Plot VAF distribution analysis."""
        all_vafs = []
        all_true_vafs = []
        
        for result in self.results.values():
            metrics = result['metrics']
            all_vafs.extend(metrics.get('vaf', []))
            all_true_vafs.extend(metrics.get('true_vaf', []))
        
        if all_vafs and all_true_vafs:
            ax.scatter(all_true_vafs, all_vafs, alpha=0.5)
            ax.plot([0, 1], [0, 1], 'r--')  # Diagonal line
            ax.set_xlabel('True VAF')
            ax.set_ylabel('Observed VAF')
            ax.set_title('VAF Correlation')

    def _plot_roc_curves(self, ax):
        """Plot ROC curves for different coverage levels."""
        for coverage, result in self.results.items():
            if 'metrics' in result:
                y_true = result['metrics'].get('is_true_variant', [])
                y_score = result['metrics'].get('vaf', [])
                
                if y_true and y_score:
                    fpr, tpr, _ = roc_curve(y_true, y_score)
                    roc_auc = auc(fpr, tpr)
                    ax.plot(fpr, tpr, label=f'{coverage} (AUC = {roc_auc:.2f})')
        
        ax.plot([0, 1], [0, 1], 'k--')
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title('ROC Curves by Coverage')
        ax.legend()

    def _plot_false_negatives(self, ax):
        """Plot false negative rates across coverage levels and depths."""
        coverages = list(self.results.keys())
        fn_rates = [res['summary'].get('false_negative_rate', 0) for res in self.results.values()]
        
        # Bar plot for overall FN rates
        ax.bar(coverages, fn_rates)
        ax.set_ylabel('False Negative Rate')
        ax.set_title('False Negative Rates by Coverage')
        
        # Add text annotations for actual counts
        for i, coverage in enumerate(coverages):
            fn_count = self.results[coverage]['summary'].get('false_negatives', 0)
            ax.text(i, fn_rates[i], f'n={fn_count}', 
                    ha='center', va='bottom')
        
        plt.setp(ax.get_xticklabels(), rotation=45)

    def _plot_fn_depth_distribution(self, ax):
        """Plot distribution of false negatives across read depths."""
        # Combine FN depth data across coverage levels
        all_fn_depths = defaultdict(int)
        for result in self.results.values():
            fn_by_depth = result['summary'].get('fn_by_depth', {})
            for depth_bin, count in fn_by_depth.items():
                all_fn_depths[depth_bin] += count
        
        # Create bar plot
        depths = sorted(all_fn_depths.keys())
        counts = [all_fn_depths[d] for d in depths]
        
        ax.bar(np.arange(len(depths)), counts)
        ax.set_xticks(np.arange(len(depths)))
        ax.set_xticklabels(depths)
        ax.set_xlabel('Read Depth')
        ax.set_ylabel('Number of False Negatives')
        ax.set_title('False Negatives Distribution by Read Depth')
        plt.setp(ax.get_xticklabels(), rotation=45)

    def _generate_error_plot(self, plot_file: str, error_message: str):
        """Generate a simple error plot when main plotting fails."""
        plt.figure(figsize=(8, 6))
        plt.text(0.5, 0.5, f'Error generating plots:\n{error_message}',
                ha='center', va='center', wrap=True)
        plt.axis('off')
        plt.savefig(plot_file)
        plt.close()

    def _generate_plots(self, plot_file: str):
        """Generate visualization plots."""
        try:
            plt.style.use('seaborn')
            fig = plt.figure(figsize=(20, 20))
            
            gs = plt.GridSpec(4, 2, figure=fig)
            ax1 = fig.add_subplot(gs[0, 0])  # Method Comparison
            ax2 = fig.add_subplot(gs[0, 1])  # Coverage Impact
            ax3 = fig.add_subplot(gs[1, :])  # Error Profiles
            ax4 = fig.add_subplot(gs[2, 0])  # VAF Analysis
            ax5 = fig.add_subplot(gs[2, 1])  # ROC Curves
            ax6 = fig.add_subplot(gs[3, 0])  # False Negative Rates
            ax7 = fig.add_subplot(gs[3, 1])  # FN Depth Distribution
            
            self._plot_method_comparison(ax1)
            self._plot_coverage_impact(ax2)
            self._plot_error_profiles(ax3)
            self._plot_vaf_analysis(ax4)
            self._plot_roc_curves(ax5)
            self._plot_false_negatives(ax6)
            self._plot_fn_depth_distribution(ax7)
            
            plt.tight_layout()
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logging.error(f"Error generating plots: {str(e)}")
            self._generate_error_plot(plot_file, str(e))

    def _write_text_report(self, report_file: str):
        """Write validation summary to a text report file."""
        with open(report_file, 'w') as f:
            f.write("Validation Report\n")
            f.write("=" * 20 + "\n\n")
            
            f.write("Summary of Validation Metrics:\n")
            for coverage, metrics in self.results.items():
                f.write(f"\nCoverage Level: {coverage}\n")
                f.write("-" * 20 + "\n")
                for key, value in metrics['summary'].items():
                    f.write(f"{key}: {value}\n")
            
        logging.info(f"Report written to {report_file}")

    def _validate_input_files(self, pileup_file: str) -> bool:
        """Validate input file existence and size."""
        if not os.path.exists(pileup_file):
            logging.error(f"Pileup file not found: {pileup_file}")
            return False
        if os.path.getsize(pileup_file) == 0:
            logging.error(f"Pileup file is empty: {pileup_file}")
            return False
        return True

    def _load_truth_data(self, truth_file: Optional[str]) -> Dict:
        """Load and parse truth data file."""
        true_variants = {}
        if truth_file and os.path.exists(truth_file):
            try:
                df = pd.read_csv(truth_file, sep='\s+')
                for _, row in df.iterrows():
                    pos = str(row['position'])
                    true_variants[f"{row['chromosome']}:{pos}"] = {
                        'alt_allele': row['alt_allele'],
                        'vaf': row.get('vaf', None),
                        'mutation_type': row.get('mutation_type', None)
                    }
            except Exception as e:
                logging.error(f"Error reading truth file: {str(e)}")
        return true_variants

    def _get_empty_results(self) -> Dict:
        """Return empty results structure."""
        return {
            'metrics': defaultdict(list),
            'summary': {
                'total_positions': 0,
                'mean_depth': 0,
                'median_depth': 0,
                'calls': Counter(),
                'accuracy': 0,
                'positions_evaluated': 0,
                'mutation_counts': {},
                'depth_accuracy': {},
                'vaf_accuracy': {},
                'false_negative_rate': 0,
                'false_negatives': 0,
                'fn_by_depth': {}
            },
            'standard_metrics': defaultdict(list),
            'standard_summary': {
                'total_positions': 0,
                'mean_depth': 0,
                'median_depth': 0,
                'calls': Counter(),
                'accuracy': 0,
                'false_negative_rate': 0,
                'false_negatives': 0,
                'fn_by_depth': {}
            }
        }

def main():
    """Main execution function."""
    try:
        logging.basicConfig(level=logging.INFO)
        logging.info("Starting validation script")
        
        # Set up argument parser
        parser = argparse.ArgumentParser()
        parser.add_argument('--pileups', nargs='+', help='Input pileup files')
        parser.add_argument('--sim-outputs', nargs='+', help='Input simulation output files')
        parser.add_argument('--mafs', nargs='+', help='Input MAF files')
        parser.add_argument('--report', help='Output report file')
        parser.add_argument('--plots', help='Output plots file')
        parser.add_argument('--method-comparison', help='Method comparison output file')
        parser.add_argument('--error-profiles', help='Error profiles output file')
        parser.add_argument('--coverage-impact', help='Coverage impact output file')
        parser.add_argument('--standard-min-depth', type=int, default=10)
        parser.add_argument('--standard-min-vaf', type=float, default=0.05)
        parser.add_argument('--multinomial-alpha', type=float, default=0.01)
        parser.add_argument('--multinomial-min-depth', type=int, default=10)
        parser.add_argument('--coverage-levels', nargs='+', help='Coverage level names')
        args = parser.parse_args()
        
        if not args.report or not args.plots or not args.method_comparison or not args.error_profiles or not args.coverage_impact:
            raise ValueError("One or more output arguments are missing. Please ensure all output files are specified.")

        print("Pileups:", args.pileups)
        print("Simulation Outputs:", args.sim_outputs)
        print("Coverage Levels:", args.coverage_levels)
        
        # Create output directories
        Path(args.report).parent.mkdir(parents=True, exist_ok=True)
        Path(os.path.dirname(args.method_comparison)).mkdir(parents=True, exist_ok=True)
        
        # Get input files with original coverage level mapping
        pileup_files = {
            coverage: str(pileup)
            for coverage, pileup in zip(
                args.coverage_levels,
                args.pileups
            )
        }
        
        truth_files = {
            coverage: str(truth)
            for coverage, truth in zip(
                args.coverage_levels,
                args.sim_outputs
            )
        }
        
        # Initialize and run validator
        validator = MultinomialValidator(
            alpha=args.multinomial_alpha,
            min_depth=args.multinomial_min_depth
        )
        
        validator.analyze_all_files(pileup_files, truth_files)
        
        # Generate reports and outputs
        validator.generate_report(
            str(args.report),
            str(args.plots)
        )
        
        validator.write_comparison_outputs(
            os.path.dirname(str(args.method_comparison))
        )
        
        logging.info("Validation completed successfully")
        
    except Exception as e:
        logging.error(f"Error in validation: {str(e)}", exc_info=True)
        raise

if __name__ == '__main__':
    main()
