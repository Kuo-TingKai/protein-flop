#!/usr/bin/env python3
"""
Main script for protein conformation transition analysis using TDA.
This script demonstrates the analysis of beta-hairpin peptide folding/unfolding
using topological data analysis methods.
"""

import numpy as np
import argparse
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.data_processing import generate_synthetic_conformation_data
from src.analysis import ProteinConformationAnalyzer


def main():
    """
    Main function to run protein conformation analysis.
    """
    parser = argparse.ArgumentParser(
        description='Analyze protein conformation transitions using TDA'
    )
    parser.add_argument(
        '--n-samples',
        type=int,
        default=1000,
        help='Number of conformations to generate (default: 1000)'
    )
    parser.add_argument(
        '--n-residues',
        type=int,
        default=6,
        help='Number of residues in turn region (default: 6)'
    )
    parser.add_argument(
        '--max-dim',
        type=int,
        default=2,
        help='Maximum homology dimension (default: 2)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='output',
        help='Output directory for plots (default: output)'
    )
    parser.add_argument(
        '--no-show',
        action='store_true',
        help='Do not display plots interactively'
    )
    parser.add_argument(
        '--use-rmsd',
        action='store_true',
        default=True,
        help='Use RMSD distance matrix (default: True)'
    )
    
    args = parser.parse_args()
    
    print("Generating synthetic conformation data...")
    print(f"  - Folded state: {args.n_samples} samples")
    print(f"  - Unfolded state: {args.n_samples} samples")
    print(f"  - Turn region: {args.n_residues} residues")
    print()
    
    # Generate synthetic data for demonstration
    # In real usage, these would come from MD simulations
    folded_coords = generate_synthetic_conformation_data(
        n_samples=args.n_samples,
        n_residues=args.n_residues,
        state='folded',
        noise_level=0.1
    )
    
    unfolded_coords = generate_synthetic_conformation_data(
        n_samples=args.n_samples,
        n_residues=args.n_residues,
        state='unfolded',
        noise_level=0.15
    )
    
    print("Initializing TDA analyzer...")
    analyzer = ProteinConformationAnalyzer(max_dim=args.max_dim)
    
    print("Analyzing conformations...")
    comparison_results = analyzer.compare_states(
        folded_coords,
        unfolded_coords,
        labels=('Folded State', 'Unfolded State'),
        use_rmsd=args.use_rmsd
    )
    
    print("Computing comparison statistics...")
    analyzer.print_comparison_summary(comparison_results)
    
    print("Generating visualizations...")
    analyzer.visualize_comparison(
        comparison_results,
        output_dir=args.output_dir,
        show=not args.no_show
    )
    
    print(f"\nAnalysis complete! Results saved to '{args.output_dir}'")
    print("\nKey findings:")
    print("  - Compare persistence diagrams to see topological differences")
    print("  - Betti curves show how topological features change with filtration")
    print("  - Differences in persistence indicate structural stability variations")


if __name__ == '__main__':
    main()

