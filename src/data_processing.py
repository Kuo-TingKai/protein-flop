"""
Data processing module for protein conformation data.
Handles loading, preprocessing, and feature extraction from MD simulation data.
"""

import numpy as np
from typing import Tuple, Optional, List
import mdtraj as md
from scipy.spatial.distance import pdist, squareform


def extract_ca_coordinates(trajectory: md.Trajectory, 
                          residue_indices: Optional[List[int]] = None) -> np.ndarray:
    """
    Extract C-alpha coordinates from MD trajectory.
    
    Args:
        trajectory: MDTraj trajectory object
        residue_indices: Optional list of residue indices to extract.
                        If None, extracts all residues.
    
    Returns:
        Array of shape (n_frames, n_residues, 3) containing C-alpha coordinates
    """
    if residue_indices is None:
        # Extract all C-alpha atoms
        ca_atoms = trajectory.topology.select('name CA')
        ca_coords = trajectory.xyz[:, ca_atoms, :]
    else:
        # Extract specific residues
        ca_coords_list = []
        for idx in residue_indices:
            residue = trajectory.topology.residue(idx)
            ca_atom = [atom for atom in residue.atoms if atom.name == 'CA'][0]
            atom_idx = ca_atom.index
            ca_coords_list.append(trajectory.xyz[:, atom_idx, :])
        ca_coords = np.stack(ca_coords_list, axis=1)
    
    return ca_coords


def extract_turn_region_coordinates(trajectory: md.Trajectory,
                                    turn_start: int,
                                    turn_end: int) -> np.ndarray:
    """
    Extract coordinates from turn region (4-6 residues) for TDA analysis.
    
    Args:
        trajectory: MDTraj trajectory object
        turn_start: Starting residue index of turn region
        turn_end: Ending residue index of turn region (exclusive)
    
    Returns:
        Array of shape (n_frames, n_residues_in_turn, 3) containing coordinates
    """
    residue_indices = list(range(turn_start, turn_end))
    return extract_ca_coordinates(trajectory, residue_indices)


def compute_rmsd_matrix(coordinates: np.ndarray) -> np.ndarray:
    """
    Compute RMSD (Root Mean Square Deviation) distance matrix between conformations.
    
    Args:
        coordinates: Array of shape (n_frames, n_atoms, 3)
    
    Returns:
        Symmetric distance matrix of shape (n_frames, n_frames)
    """
    n_frames = coordinates.shape[0]
    rmsd_matrix = np.zeros((n_frames, n_frames))
    
    for i in range(n_frames):
        for j in range(i+1, n_frames):
            # Center both conformations
            coords_i = coordinates[i] - coordinates[i].mean(axis=0)
            coords_j = coordinates[j] - coordinates[j].mean(axis=0)
            
            # Compute RMSD
            diff = coords_i - coords_j
            rmsd = np.sqrt(np.mean(diff**2))
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
    
    return rmsd_matrix


def flatten_coordinates(coordinates: np.ndarray) -> np.ndarray:
    """
    Flatten coordinate array for TDA analysis.
    
    Args:
        coordinates: Array of shape (n_frames, n_atoms, 3)
    
    Returns:
        Array of shape (n_frames, n_atoms * 3) for point cloud representation
    """
    return coordinates.reshape(coordinates.shape[0], -1)


def load_trajectory_from_file(filepath: str, 
                              top_file: Optional[str] = None) -> md.Trajectory:
    """
    Load MD trajectory from file.
    
    Args:
        filepath: Path to trajectory file (e.g., .xtc, .dcd, .pdb)
        top_file: Optional topology file if needed
    
    Returns:
        MDTraj trajectory object
    """
    if top_file:
        return md.load(filepath, top=top_file)
    else:
        return md.load(filepath)


def generate_synthetic_conformation_data(n_samples: int = 1000,
                                         n_residues: int = 6,
                                         state: str = 'folded',
                                         noise_level: float = 0.1) -> np.ndarray:
    """
    Generate synthetic conformation data for testing.
    Creates point clouds representing different conformational states.
    
    Args:
        n_samples: Number of conformations to generate
        n_residues: Number of residues in turn region
        state: 'folded' or 'unfolded' - determines base structure
        noise_level: Standard deviation of Gaussian noise
    
    Returns:
        Array of shape (n_samples, n_residues, 3) containing synthetic coordinates
    """
    np.random.seed(42)
    
    if state == 'folded':
        # Create a compact beta-hairpin turn structure
        # Residues form a tight turn with specific geometry
        base_coords = np.array([
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [3.8, 3.8, 0.0],
            [0.0, 3.8, 0.0],
            [-3.8, 3.8, 0.0],
            [-3.8, 0.0, 0.0]
        ])
        # Add slight curvature to form turn
        base_coords[:, 2] = np.sin(np.linspace(0, np.pi, n_residues)) * 2.0
    else:  # unfolded
        # Create a more extended, less structured conformation
        base_coords = np.array([
            [0.0, 0.0, 0.0],
            [4.5, 0.5, 0.5],
            [9.0, 1.0, 1.0],
            [13.5, 1.5, 0.5],
            [18.0, 2.0, 0.0],
            [22.5, 2.5, -0.5]
        ])
        # Less curvature, more random
        base_coords[:, 2] = np.random.randn(n_residues) * 1.0
    
    # Generate samples by adding noise
    samples = []
    for _ in range(n_samples):
        noise = np.random.randn(n_residues, 3) * noise_level
        sample = base_coords + noise
        samples.append(sample)
    
    return np.array(samples)

