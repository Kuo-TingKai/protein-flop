"""
Topological Data Analysis (TDA) computation module.
Implements persistent homology using Vietoris-Rips filtration.
"""

import numpy as np
from typing import Dict, Tuple, List, Optional
from ripser import ripser
from scipy.spatial.distance import pdist, squareform


def compute_persistent_homology(points: np.ndarray,
                               distance_matrix: Optional[np.ndarray] = None,
                               max_dim: int = 2,
                               metric: str = 'euclidean',
                               max_filtration: Optional[float] = None) -> Dict:
    """
    Compute persistent homology using Vietoris-Rips filtration.
    
    Args:
        points: Point cloud of shape (n_points, n_dimensions)
        distance_matrix: Optional precomputed distance matrix.
                        If None, will be computed from points.
        max_dim: Maximum homology dimension to compute (0, 1, 2, ...)
        metric: Distance metric to use ('euclidean', 'rmsd', etc.)
        max_filtration: Maximum filtration value (epsilon).
                       If None, uses maximum distance in point cloud.
    
    Returns:
        Dictionary containing persistence diagrams and other TDA results
    """
    # Compute distance matrix if not provided
    if distance_matrix is None:
        if metric == 'euclidean':
            distance_matrix = squareform(pdist(points, metric='euclidean'))
        else:
            distance_matrix = squareform(pdist(points, metric=metric))
    
    # Set max_filtration if not provided
    if max_filtration is None:
        max_filtration = np.max(distance_matrix) * 1.1
    
    # Compute persistent homology using ripser
    result = ripser(
        distance_matrix,
        maxdim=max_dim,
        thresh=max_filtration,
        metric='precomputed'
    )
    
    return {
        'diagrams': result['dgms'],
        'cocycles': result.get('cocycles', []),
        'distance_matrix': distance_matrix,
        'max_filtration': max_filtration
    }


def compute_betti_numbers(persistence_diagrams: List[np.ndarray],
                         filtration_value: float) -> Dict[int, int]:
    """
    Compute Betti numbers at a specific filtration value.
    
    Args:
        persistence_diagrams: List of persistence diagrams for each dimension
        filtration_value: Filtration value (epsilon) at which to compute Betti numbers
    
    Returns:
        Dictionary mapping dimension -> Betti number
    """
    betti_numbers = {}
    
    for dim, diagram in enumerate(persistence_diagrams):
        if diagram.size == 0:
            betti_numbers[dim] = 0
        else:
            # Count features that are born before filtration_value
            # and die after filtration_value
            alive_features = np.sum(
                (diagram[:, 0] <= filtration_value) & 
                (diagram[:, 1] > filtration_value)
            )
            betti_numbers[dim] = int(alive_features)
    
    return betti_numbers


def compute_persistence_statistics(persistence_diagrams: List[np.ndarray]) -> Dict:
    """
    Compute statistics about persistence features.
    
    Args:
        persistence_diagrams: List of persistence diagrams for each dimension
    
    Returns:
        Dictionary containing statistics for each dimension
    """
    stats = {}
    
    for dim, diagram in enumerate(persistence_diagrams):
        if diagram.size == 0:
            stats[dim] = {
                'count': 0,
                'mean_persistence': 0.0,
                'max_persistence': 0.0,
                'min_persistence': 0.0,
                'total_persistence': 0.0
            }
        else:
            persistences = diagram[:, 1] - diagram[:, 0]  # death - birth
            stats[dim] = {
                'count': len(diagram),
                'mean_persistence': float(np.mean(persistences)),
                'max_persistence': float(np.max(persistences)),
                'min_persistence': float(np.min(persistences)),
                'total_persistence': float(np.sum(persistences))
            }
    
    return stats


def compare_persistence_diagrams(diagram1: np.ndarray,
                                diagram2: np.ndarray,
                                dimension: int = 1) -> Dict:
    """
    Compare two persistence diagrams and compute differences.
    
    Args:
        diagram1: First persistence diagram
        diagram2: Second persistence diagram
        dimension: Homology dimension being compared
    
    Returns:
        Dictionary containing comparison metrics
    """
    if diagram1.size == 0 and diagram2.size == 0:
        return {
            'count_diff': 0,
            'mean_persistence_diff': 0.0,
            'max_persistence_diff': 0.0
        }
    
    if diagram1.size == 0:
        persistences1 = np.array([0.0])
    else:
        persistences1 = diagram1[:, 1] - diagram1[:, 0]
    
    if diagram2.size == 0:
        persistences2 = np.array([0.0])
    else:
        persistences2 = diagram2[:, 1] - diagram2[:, 0]
    
    comparison = {
        'count_diff': len(diagram1) - len(diagram2),
        'mean_persistence_diff': float(np.mean(persistences1) - np.mean(persistences2)),
        'max_persistence_diff': float(np.max(persistences1) - np.max(persistences2)),
        'diagram1_stats': {
            'count': len(diagram1),
            'mean_persistence': float(np.mean(persistences1)),
            'max_persistence': float(np.max(persistences1))
        },
        'diagram2_stats': {
            'count': len(diagram2),
            'mean_persistence': float(np.mean(persistences2)),
            'max_persistence': float(np.max(persistences2))
        }
    }
    
    return comparison


def filter_persistence_diagram(diagram: np.ndarray,
                              min_persistence: float = 0.0) -> np.ndarray:
    """
    Filter persistence diagram to keep only features with persistence >= min_persistence.
    
    Args:
        diagram: Persistence diagram
        min_persistence: Minimum persistence threshold
    
    Returns:
        Filtered persistence diagram
    """
    if diagram.size == 0:
        return diagram
    
    persistences = diagram[:, 1] - diagram[:, 0]
    mask = persistences >= min_persistence
    return diagram[mask]

