"""
Flop computation module for implementing Atiyah Flop and protein
conformation Flop operations.

This module provides computational implementations of MMP flip/flop
operations and their analogies to protein conformation transitions.
"""

import numpy as np
from typing import Dict, Optional
from scipy.spatial.distance import pdist, squareform


def sample_atiyah_flop_before(n_samples: int = 1000,
                              radius: float = 1.0,
                              noise_level: float = 0.05) -> np.ndarray:
    """
    Sample points from the resolved variety X~ before Atiyah Flop.

    The variety X = {(x,y,z,w) in C^4 | xy = zw} is resolved to X~
    which contains an exceptional curve E ≅ P^1 (topologically S^2)
    that contracts to the singularity P.
    
    Args:
        n_samples: Number of points to sample
        radius: Radius of the sampling neighborhood around the singularity
        noise_level: Gaussian noise level for sampling
    
    Returns:
        Array of shape (n_samples, 6) representing points in R^6 (3D complex = 6D real)
    """
    np.random.seed(42)
    
    # Sample from the exceptional curve E (P^1, topologically S^2)
    # In the resolution X~, E is parameterized by [u:v] in P^1
    n_exceptional = n_samples // 2
    n_regular = n_samples - n_exceptional
    
    points = []
    
    # Sample points on the exceptional curve E
    # E is a sphere S^2 embedded in the resolution
    for _ in range(n_exceptional):
        # Sample uniformly on S^2
        u = np.random.randn(3)
        u = u / np.linalg.norm(u)
        
        # Map to the resolution X~ near the singularity
        # The exceptional curve is at the origin in the base, but extends in P^1
        # We represent this as a small sphere around the origin
        point = radius * 0.1 * u + np.random.randn(6) * noise_level * radius
        points.append(point)
    
    # Sample regular points near the singularity (satisfying xy = zw)
    for _ in range(n_regular):
        # Sample (x, y, z, w) satisfying xy = zw
        # Use parameterization: x = t, y = s, z = t, w = s (one solution)
        t = np.random.randn() * radius * 0.5
        s = np.random.randn() * radius * 0.5
        
        # Map to real coordinates (treating complex as R^2)
        x_re, x_im = t, 0.0
        y_re, y_im = s, 0.0
        z_re, z_im = t, 0.0

        # Add noise and map to 6D real space
        point = (np.array([x_re, x_im, y_re, y_im, z_re, z_im]) +
                 np.random.randn(6) * noise_level * radius)
        points.append(point)
    
    return np.array(points)


def sample_atiyah_flop_after(n_samples: int = 1000,
                             radius: float = 1.0,
                             noise_level: float = 0.05) -> np.ndarray:
    """
    Sample points from the resolved variety X~' after Atiyah Flop.
    
    After the Flop, we have a different resolution X~' with exceptional curve E'
    that has different algebraic properties (self-intersection +1 vs -1).
    
    Args:
        n_samples: Number of points to sample
        radius: Radius of the sampling neighborhood
        noise_level: Gaussian noise level for sampling
    
    Returns:
        Array of shape (n_samples, 6) representing points in R^6
    """
    np.random.seed(43)  # Different seed for different sampling
    
    # Similar structure but with different embedding
    n_exceptional = n_samples // 2
    n_regular = n_samples - n_exceptional
    
    points = []
    
    # Sample points on the exceptional curve E' (after Flop)
    # The key difference is in how the curve is embedded
    for _ in range(n_exceptional):
        # Sample on S^2 but with different orientation/embedding
        u = np.random.randn(3)
        u = u / np.linalg.norm(u)
        
        # The Flop changes the embedding, which we model as a rotation
        # and slight change in the local geometry
        rotation = np.array([[1, 0, 0, 0, 0, 0],
                             [0, 1, 0, 0, 0, 0],
                             [0, 0, -1, 0, 0, 0],  # Sign change represents Flop
                             [0, 0, 0, -1, 0, 0],
                             [0, 0, 0, 0, 1, 0],
                             [0, 0, 0, 0, 0, 1]])
        
        point_base = radius * 0.1 * np.concatenate([u, u[:3]])  # 6D point
        point = rotation @ point_base + np.random.randn(6) * noise_level * radius
        points.append(point)
    
    # Sample regular points (satisfying xy = zw, but with different parameterization)
    for _ in range(n_regular):
        # Different parameterization after Flop: x = t, y = -s, z = -t, w = s
        t = np.random.randn() * radius * 0.5
        s = np.random.randn() * radius * 0.5
        
        x_re, x_im = t, 0.0
        y_re, y_im = -s, 0.0  # Sign change
        z_re, z_im = -t, 0.0  # Sign change

        point = (np.array([x_re, x_im, y_re, y_im, z_re, z_im]) +
                 np.random.randn(6) * noise_level * radius)
        points.append(point)
    
    return np.array(points)


def protein_flop_analogy_folded(n_samples: int = 1000,
                                n_residues: int = 6,
                                noise_level: float = 0.1) -> np.ndarray:
    """
    Generate protein conformation data representing the "before Flop" state.
    For beta-hairpin peptide, this represents the folded state.
    
    This is analogous to X~ before the Flop operation.
    
    Args:
        n_samples: Number of conformations to generate
        n_residues: Number of residues in the turn region
        noise_level: Noise level for sampling
    
    Returns:
        Array of shape (n_samples, n_residues, 3) containing C-alpha coordinates
    """
    np.random.seed(42)
    
    # Create a compact beta-hairpin turn structure (folded state)
    # This represents the stable, well-structured state before "Flop"
    base_coords = np.array([
        [0.0, 0.0, 0.0],
        [3.8, 0.0, 0.0],
        [3.8, 3.8, 0.0],
        [0.0, 3.8, 0.0],
        [-3.8, 3.8, 0.0],
        [-3.8, 0.0, 0.0]
    ])
    
    # Add curvature to form a tight turn (like the exceptional curve E)
    angles = np.linspace(0, np.pi, n_residues)
    base_coords[:, 2] = np.sin(angles) * 2.0  # Z-coordinate for turn
    
    # Generate samples
    samples = []
    for _ in range(n_samples):
        noise = np.random.randn(n_residues, 3) * noise_level
        sample = base_coords + noise
        samples.append(sample)
    
    return np.array(samples)


def protein_flop_analogy_unfolded(n_samples: int = 1000,
                                  n_residues: int = 6,
                                  noise_level: float = 0.15) -> np.ndarray:
    """
    Generate protein conformation data representing the "after Flop" state.
    For beta-hairpin peptide, this represents an unfolded/intermediate state.
    
    This is analogous to X~' after the Flop operation.
    
    Args:
        n_samples: Number of conformations to generate
        n_residues: Number of residues in the turn region
        noise_level: Noise level for sampling
    
    Returns:
        Array of shape (n_samples, n_residues, 3) containing C-alpha coordinates
    """
    np.random.seed(43)  # Different seed
    
    # Create a more extended, less structured conformation (unfolded state)
    # This represents the state after "Flop" - topologically equivalent but
    # with different local geometry (like E' with different self-intersection)
    base_coords = np.array([
        [0.0, 0.0, 0.0],
        [4.5, 0.5, 0.5],
        [9.0, 1.0, 1.0],
        [13.5, 1.5, 0.5],
        [18.0, 2.0, 0.0],
        [22.5, 2.5, -0.5]
    ])
    
    # Less curvature, more extended (like E' with different embedding)
    base_coords[:, 2] = np.random.randn(n_residues) * 1.0
    
    # Generate samples
    samples = []
    for _ in range(n_samples):
        noise = np.random.randn(n_residues, 3) * noise_level
        sample = base_coords + noise
        samples.append(sample)
    
    return np.array(samples)


def compute_flop_difference(points_before: np.ndarray,
                           points_after: np.ndarray,
                           metric: str = 'euclidean') -> Dict:
    """
    Compute topological and geometric differences between before and after Flop.
    
    This function quantifies the effect of the Flop operation by comparing
    the two point clouds using various metrics.
    
    Args:
        points_before: Point cloud before Flop (shape: n_samples, n_dim)
        points_after: Point cloud after Flop (shape: n_samples, n_dim)
        metric: Distance metric to use
    
    Returns:
        Dictionary containing comparison metrics
    """
    # Flatten if needed
    if len(points_before.shape) == 3:
        points_before = points_before.reshape(points_before.shape[0], -1)
    if len(points_after.shape) == 3:
        points_after = points_after.reshape(points_after.shape[0], -1)
    
    # Compute distance matrices
    dist_before = squareform(pdist(points_before, metric=metric))
    dist_after = squareform(pdist(points_after, metric=metric))
    
    # Statistical comparisons
    comparison = {
        'mean_distance_before': float(np.mean(dist_before)),
        'mean_distance_after': float(np.mean(dist_after)),
        'std_distance_before': float(np.std(dist_before)),
        'std_distance_after': float(np.std(dist_after)),
        'max_distance_before': float(np.max(dist_before)),
        'max_distance_after': float(np.max(dist_after)),
        'dimension': points_before.shape[1],
        'n_samples_before': points_before.shape[0],
        'n_samples_after': points_after.shape[0]
    }
    
    # Compute centroid distances
    centroid_before = np.mean(points_before, axis=0)
    centroid_after = np.mean(points_after, axis=0)
    comparison['centroid_distance'] = float(np.linalg.norm(centroid_before - centroid_after))
    
    # Compute spread (variance)
    spread_before = np.var(points_before, axis=0).sum()
    spread_after = np.var(points_after, axis=0).sum()
    comparison['spread_before'] = float(spread_before)
    comparison['spread_after'] = float(spread_after)
    comparison['spread_ratio'] = float(spread_after / spread_before) if spread_before > 0 else 0.0
    
    return comparison


def apply_protein_flop_transformation(coordinates: np.ndarray,
                                     flop_center: Optional[np.ndarray] = None,
                                     flop_axis: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Apply a Flop-like transformation to protein coordinates.
    
    This simulates the effect of a Flop operation on protein conformation,
    which involves local reorganization while preserving overall topology.
    
    Args:
        coordinates: Input coordinates (n_frames, n_residues, 3)
        flop_center: Center of the Flop operation (default: centroid)
        flop_axis: Axis for the Flop transformation (default: principal axis)
    
    Returns:
        Transformed coordinates
    """
    if flop_center is None:
        flop_center = np.mean(coordinates, axis=(0, 1))
    
    # Center coordinates
    centered = coordinates - flop_center
    
    # Apply Flop transformation: local rotation + reflection
    # This mimics the "cut and paste" operation of Flop
    if flop_axis is None:
        # Use first principal component as axis
        coords_flat = centered.reshape(-1, 3)
        cov = np.cov(coords_flat.T)
        _, eigenvecs = np.linalg.eigh(cov)
        flop_axis = eigenvecs[:, -1]  # Principal axis
    
    # Create rotation matrix for Flop (180 degree rotation around axis)
    # Combined with a reflection to change orientation
    angle = np.pi
    axis = flop_axis / np.linalg.norm(flop_axis)
    
    # Rodrigues' rotation formula
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
    
    # Apply transformation
    transformed = np.zeros_like(centered)
    for i in range(centered.shape[0]):
        for j in range(centered.shape[1]):
            transformed[i, j] = R @ centered[i, j]
    
    # Add back center
    return transformed + flop_center

