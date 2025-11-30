"""
Main analysis module that combines data processing, TDA computation, and visualization.
Implements the complete workflow for analyzing protein conformation transitions.
"""

import numpy as np
from typing import Dict, Tuple, Optional
import matplotlib.pyplot as plt

from .data_processing import (
    extract_turn_region_coordinates,
    flatten_coordinates,
    compute_rmsd_matrix,
    generate_synthetic_conformation_data
)
from .tda_computation import (
    compute_persistent_homology,
    compute_betti_numbers,
    compute_persistence_statistics,
    compare_persistence_diagrams
)
from .visualization import (
    plot_persistence_diagram,
    compare_persistence_diagrams_plot,
    plot_betti_curve
)


class ProteinConformationAnalyzer:
    """
    Main class for analyzing protein conformation transitions using TDA.
    """
    
    def __init__(self, max_dim: int = 2):
        """
        Initialize analyzer.
        
        Args:
            max_dim: Maximum homology dimension to compute
        """
        self.max_dim = max_dim
        self.results_folded = None
        self.results_unfolded = None
    
    def analyze_conformation(self,
                           coordinates: np.ndarray,
                           use_rmsd: bool = True,
                           max_filtration: Optional[float] = None) -> Dict:
        """
        Analyze a single conformation state.
        
        Args:
            coordinates: Array of shape (n_frames, n_residues, 3)
            use_rmsd: Whether to use RMSD distance matrix
            max_filtration: Maximum filtration value
        
        Returns:
            Dictionary containing TDA results
        """
        # Flatten coordinates for point cloud representation
        points = flatten_coordinates(coordinates)
        
        # Compute distance matrix
        if use_rmsd:
            distance_matrix = compute_rmsd_matrix(coordinates)
        else:
            distance_matrix = None
        
        # Compute persistent homology
        tda_results = compute_persistent_homology(
            points,
            distance_matrix=distance_matrix,
            max_dim=self.max_dim,
            metric='precomputed' if use_rmsd else 'euclidean',
            max_filtration=max_filtration
        )
        
        # Compute statistics
        stats = compute_persistence_statistics(tda_results['diagrams'])
        
        return {
            'tda_results': tda_results,
            'statistics': stats,
            'coordinates': coordinates,
            'points': points
        }
    
    def compare_states(self,
                      folded_coordinates: np.ndarray,
                      unfolded_coordinates: np.ndarray,
                      labels: Tuple[str, str] = ('Folded', 'Unfolded'),
                      use_rmsd: bool = True) -> Dict:
        """
        Compare two conformational states.
        
        Args:
            folded_coordinates: Coordinates for folded state
            unfolded_coordinates: Coordinates for unfolded state
            labels: Labels for the two states
            use_rmsd: Whether to use RMSD distance matrix
        
        Returns:
            Dictionary containing comparison results
        """
        # Analyze both states
        self.results_folded = self.analyze_conformation(
            folded_coordinates, use_rmsd=use_rmsd
        )
        self.results_unfolded = self.analyze_conformation(
            unfolded_coordinates, use_rmsd=use_rmsd
        )
        
        # Compare persistence diagrams for each dimension
        comparisons = {}
        for dim in range(self.max_dim + 1):
            diagram1 = self.results_folded['tda_results']['diagrams'][dim]
            diagram2 = self.results_unfolded['tda_results']['diagrams'][dim]
            
            comparisons[dim] = compare_persistence_diagrams(
                diagram1, diagram2, dimension=dim
            )
        
        return {
            'folded_results': self.results_folded,
            'unfolded_results': self.results_unfolded,
            'comparisons': comparisons,
            'labels': labels
        }
    
    def visualize_comparison(self,
                           comparison_results: Dict,
                           output_dir: Optional[str] = None,
                           show: bool = True) -> None:
        """
        Visualize comparison results.
        
        Args:
            comparison_results: Results from compare_states()
            output_dir: Optional directory to save plots
            show: Whether to display plots
        """
        labels = comparison_results['labels']
        folded_results = comparison_results['folded_results']
        unfolded_results = comparison_results['unfolded_results']
        
        # Plot persistence diagrams for each dimension
        for dim in range(self.max_dim + 1):
            diagram1 = folded_results['tda_results']['diagrams'][dim]
            diagram2 = unfolded_results['tda_results']['diagrams'][dim]
            
            # Side-by-side comparison
            fig = compare_persistence_diagrams_plot(
                diagram1, diagram2,
                dimension=dim,
                labels=labels,
                title=f'Persistence Diagram Comparison - Dimension {dim}',
                show=show
            )
            
            if output_dir:
                import os
                os.makedirs(output_dir, exist_ok=True)
                fig.savefig(f'{output_dir}/persistence_diagram_dim{dim}.png', 
                           dpi=300, bbox_inches='tight')
                plt.close(fig)
        
        # Plot Betti curves
        for dim in range(1, min(3, self.max_dim + 1)):  # Focus on dim 1 and 2
            diagram1 = folded_results['tda_results']['diagrams'][dim]
            diagram2 = unfolded_results['tda_results']['diagrams'][dim]
            
            # Determine filtration range
            max_filt = max(
                folded_results['tda_results']['max_filtration'],
                unfolded_results['tda_results']['max_filtration']
            )
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
            
            plot_betti_curve(
                folded_results['tda_results']['diagrams'],
                dimension=dim,
                filtration_range=(0, max_filt),
                title=f'{labels[0]} - Betti {dim}',
                ax=ax1,
                show=False
            )
            
            plot_betti_curve(
                unfolded_results['tda_results']['diagrams'],
                dimension=dim,
                filtration_range=(0, max_filt),
                title=f'{labels[1]} - Betti {dim}',
                ax=ax2,
                show=False
            )
            
            plt.tight_layout()
            if show:
                plt.show()
            
            if output_dir:
                fig.savefig(f'{output_dir}/betti_curve_dim{dim}.png',
                           dpi=300, bbox_inches='tight')
                plt.close(fig)
    
    def print_comparison_summary(self, comparison_results: Dict) -> None:
        """
        Print summary of comparison results.
        
        Args:
            comparison_results: Results from compare_states()
        """
        print("=" * 80)
        print("PROTEIN CONFORMATION TRANSITION ANALYSIS")
        print("=" * 80)
        print()
        
        labels = comparison_results['labels']
        comparisons = comparison_results['comparisons']
        
        for dim in range(self.max_dim + 1):
            comp = comparisons[dim]
            print(f"Dimension {dim} (H_{dim}):")
            print(f"  {labels[0]}:")
            print(f"    - Feature count: {comp['diagram1_stats']['count']}")
            print(f"    - Mean persistence: {comp['diagram1_stats']['mean_persistence']:.4f}")
            print(f"    - Max persistence: {comp['diagram1_stats']['max_persistence']:.4f}")
            print(f"  {labels[1]}:")
            print(f"    - Feature count: {comp['diagram2_stats']['count']}")
            print(f"    - Mean persistence: {comp['diagram2_stats']['mean_persistence']:.4f}")
            print(f"    - Max persistence: {comp['diagram2_stats']['max_persistence']:.4f}")
            print(f"  Differences:")
            print(f"    - Count difference: {comp['count_diff']}")
            print(f"    - Mean persistence difference: {comp['mean_persistence_diff']:.4f}")
            print(f"    - Max persistence difference: {comp['max_persistence_diff']:.4f}")
            print()

