#!/usr/bin/env python3
"""
Example usage script demonstrating how to use the protein-flop analysis tools.
This script shows how to analyze protein conformation transitions step by step.
"""

import numpy as np
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.data_processing import (
    generate_synthetic_conformation_data,
    flatten_coordinates,
    compute_rmsd_matrix
)
from src.tda_computation import (
    compute_persistent_homology,
    compute_persistence_statistics,
    compare_persistence_diagrams
)
from src.visualization import (
    plot_persistence_diagram,
    compare_persistence_diagrams_plot
)
from src.analysis import ProteinConformationAnalyzer


def example_basic_analysis():
    """
    Basic example: Analyze a single conformation state.
    """
    print("=" * 80)
    print("Example 1: Basic Analysis of Single Conformation State")
    print("=" * 80)
    
    # Generate synthetic data
    coordinates = generate_synthetic_conformation_data(
        n_samples=500,
        n_residues=6,
        state='folded',
        noise_level=0.1
    )
    
    # Flatten for point cloud
    points = flatten_coordinates(coordinates)
    
    # Compute persistent homology
    print("Computing persistent homology...")
    tda_results = compute_persistent_homology(
        points,
        max_dim=2,
        metric='euclidean'
    )
    
    # Print statistics
    stats = compute_persistence_statistics(tda_results['diagrams'])
    print("\nPersistence Statistics:")
    for dim, stat in stats.items():
        print(f"  Dimension {dim}:")
        print(f"    - Feature count: {stat['count']}")
        print(f"    - Mean persistence: {stat['mean_persistence']:.4f}")
        print(f"    - Max persistence: {stat['max_persistence']:.4f}")
    
    # Plot persistence diagram
    print("\nPlotting persistence diagram...")
    plot_persistence_diagram(
        tda_results['diagrams'],
        max_dim=2,
        title="Persistence Diagram - Folded State"
    )


def example_comparison_analysis():
    """
    Example: Compare two conformational states.
    """
    print("\n" + "=" * 80)
    print("Example 2: Comparison of Folded vs Unfolded States")
    print("=" * 80)
    
    # Generate data for both states
    folded_coords = generate_synthetic_conformation_data(
        n_samples=1000,
        n_residues=6,
        state='folded',
        noise_level=0.1
    )
    
    unfolded_coords = generate_synthetic_conformation_data(
        n_samples=1000,
        n_residues=6,
        state='unfolded',
        noise_level=0.15
    )
    
    # Use the analyzer class
    analyzer = ProteinConformationAnalyzer(max_dim=2)
    
    print("Analyzing both states...")
    comparison_results = analyzer.compare_states(
        folded_coords,
        unfolded_coords,
        labels=('Folded', 'Unfolded'),
        use_rmsd=True
    )
    
    # Print summary
    analyzer.print_comparison_summary(comparison_results)
    
    # Visualize
    print("\nGenerating comparison plots...")
    analyzer.visualize_comparison(
        comparison_results,
        output_dir='output',
        show=True
    )


def example_custom_analysis():
    """
    Example: Custom analysis with specific parameters.
    """
    print("\n" + "=" * 80)
    print("Example 3: Custom Analysis with RMSD Distance Matrix")
    print("=" * 80)
    
    # Generate data
    coordinates = generate_synthetic_conformation_data(
        n_samples=800,
        n_residues=5,
        state='folded',
        noise_level=0.12
    )
    
    # Compute RMSD matrix
    print("Computing RMSD distance matrix...")
    rmsd_matrix = compute_rmsd_matrix(coordinates)
    print(f"  RMSD matrix shape: {rmsd_matrix.shape}")
    print(f"  Mean RMSD: {np.mean(rmsd_matrix):.4f}")
    print(f"  Max RMSD: {np.max(rmsd_matrix):.4f}")
    
    # Use RMSD matrix for TDA
    points = flatten_coordinates(coordinates)
    tda_results = compute_persistent_homology(
        points,
        distance_matrix=rmsd_matrix,
        max_dim=2,
        metric='precomputed',
        max_filtration=np.max(rmsd_matrix) * 0.8
    )
    
    # Focus on dimension 1 (loops/holes)
    dim1_diagram = tda_results['diagrams'][1]
    if dim1_diagram.size > 0:
        persistences = dim1_diagram[:, 1] - dim1_diagram[:, 0]
        print(f"\nDimension 1 (H_1) - Loops/Holes:")
        print(f"  - Number of features: {len(dim1_diagram)}")
        print(f"  - Mean persistence: {np.mean(persistences):.4f}")
        print(f"  - Most persistent feature: {np.max(persistences):.4f}")
    
    # Plot
    plot_persistence_diagram(
        tda_results['diagrams'],
        max_dim=2,
        title="Persistence Diagram with RMSD Distance"
    )


if __name__ == '__main__':
    print("Protein Flop - Example Usage Script")
    print("=" * 80)
    print()
    
    # Run examples
    try:
        example_basic_analysis()
        example_comparison_analysis()
        example_custom_analysis()
        
        print("\n" + "=" * 80)
        print("All examples completed successfully!")
        print("=" * 80)
        
    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\n\nError during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

