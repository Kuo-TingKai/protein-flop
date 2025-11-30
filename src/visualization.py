"""
Visualization module for TDA results.
Creates persistence diagrams, Betti curves, and comparison plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from typing import List, Dict, Optional, Tuple
import matplotlib.patches as mpatches


def plot_persistence_diagram(persistence_diagrams: List[np.ndarray],
                            max_dim: int = 2,
                            title: str = "Persistence Diagram",
                            ax: Optional[plt.Axes] = None,
                            show: bool = True) -> plt.Axes:
    """
    Plot persistence diagram for multiple dimensions.
    
    Args:
        persistence_diagrams: List of persistence diagrams for each dimension
        max_dim: Maximum dimension to plot
        title: Plot title
        ax: Optional matplotlib axes to plot on
        show: Whether to display the plot
    
    Returns:
        Matplotlib axes object
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    
    colors = ['blue', 'red', 'green', 'orange', 'purple']
    
    max_birth = 0
    max_death = 0
    
    for dim in range(min(max_dim + 1, len(persistence_diagrams))):
        diagram = persistence_diagrams[dim]
        if diagram.size > 0:
            births = diagram[:, 0]
            deaths = diagram[:, 1]
            max_birth = max(max_birth, np.max(births))
            max_death = max(max_death, np.max(deaths))
            
            ax.scatter(births, deaths, 
                      c=colors[dim % len(colors)],
                      label=f'$H_{dim}$',
                      alpha=0.6,
                      s=50)
    
    # Plot diagonal line (y = x)
    max_val = max(max_birth, max_death) * 1.1
    ax.plot([0, max_val], [0, max_val], 
           'k--', alpha=0.5, linewidth=1, label='Diagonal')
    
    ax.set_xlabel('Birth', fontsize=12)
    ax.set_ylabel('Death', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    if show:
        plt.tight_layout()
        plt.show()
    
    return ax


def plot_betti_curve(persistence_diagrams: List[np.ndarray],
                    dimension: int,
                    filtration_range: Optional[Tuple[float, float]] = None,
                    n_points: int = 100,
                    title: Optional[str] = None,
                    ax: Optional[plt.Axes] = None,
                    show: bool = True) -> plt.Axes:
    """
    Plot Betti number as a function of filtration value.
    
    Args:
        persistence_diagrams: List of persistence diagrams
        dimension: Homology dimension to plot
        filtration_range: Optional (min, max) range for filtration values
        n_points: Number of points to sample
        title: Plot title
        ax: Optional matplotlib axes
        show: Whether to display the plot
    
    Returns:
        Matplotlib axes object
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    diagram = persistence_diagrams[dimension]
    
    if diagram.size == 0:
        ax.plot([0, 1], [0, 0], 'b-', linewidth=2)
        ax.set_ylabel(f'$\\beta_{dimension}$', fontsize=12)
        ax.set_xlabel('Filtration Value (ε)', fontsize=12)
        if title:
            ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        if show:
            plt.tight_layout()
            plt.show()
        return ax
    
    # Determine filtration range
    if filtration_range is None:
        min_val = np.min(diagram[:, 0])
        max_val = np.max(diagram[:, 1])
        filtration_range = (0, max_val * 1.1)
    
    epsilon_values = np.linspace(filtration_range[0], filtration_range[1], n_points)
    betti_values = []
    
    for eps in epsilon_values:
        # Count features alive at this filtration value
        alive = np.sum((diagram[:, 0] <= eps) & (diagram[:, 1] > eps))
        betti_values.append(alive)
    
    ax.plot(epsilon_values, betti_values, 'b-', linewidth=2, label=f'$\\beta_{dimension}$')
    ax.set_ylabel(f'Betti Number $\\beta_{dimension}$', fontsize=12)
    ax.set_xlabel('Filtration Value (ε)', fontsize=12)
    
    if title is None:
        title = f'Betti Curve for Dimension {dimension}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    if show:
        plt.tight_layout()
        plt.show()
    
    return ax


def compare_persistence_diagrams_plot(diagram1: np.ndarray,
                                      diagram2: np.ndarray,
                                      dimension: int = 1,
                                      labels: Tuple[str, str] = ('State 1', 'State 2'),
                                      title: Optional[str] = None,
                                      show: bool = True) -> plt.Figure:
    """
    Compare two persistence diagrams side by side.
    
    Args:
        diagram1: First persistence diagram
        diagram2: Second persistence diagram
        dimension: Homology dimension
        labels: Labels for the two states
        title: Plot title
        show: Whether to display the plot
    
    Returns:
        Matplotlib figure object
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    colors = ['blue', 'red', 'green', 'orange', 'purple']
    color = colors[dimension % len(colors)]
    
    # Plot first diagram
    if diagram1.size > 0:
        births1 = diagram1[:, 0]
        deaths1 = diagram1[:, 1]
        max_val1 = max(np.max(births1), np.max(deaths1)) * 1.1
        ax1.scatter(births1, deaths1, c=color, alpha=0.6, s=50)
        ax1.plot([0, max_val1], [0, max_val1], 'k--', alpha=0.5, linewidth=1)
        ax1.set_xlim(0, max_val1)
        ax1.set_ylim(0, max_val1)
    else:
        max_val1 = 1.0
        ax1.set_xlim(0, max_val1)
        ax1.set_ylim(0, max_val1)
    
    ax1.set_xlabel('Birth', fontsize=12)
    ax1.set_ylabel('Death', fontsize=12)
    ax1.set_title(f'{labels[0]} - $H_{dimension}$', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal', adjustable='box')
    
    # Plot second diagram
    if diagram2.size > 0:
        births2 = diagram2[:, 0]
        deaths2 = diagram2[:, 1]
        max_val2 = max(np.max(births2), np.max(deaths2)) * 1.1
        ax2.scatter(births2, deaths2, c=color, alpha=0.6, s=50)
        ax2.plot([0, max_val2], [0, max_val2], 'k--', alpha=0.5, linewidth=1)
        ax2.set_xlim(0, max_val2)
        ax2.set_ylim(0, max_val2)
    else:
        max_val2 = 1.0
        ax2.set_xlim(0, max_val2)
        ax2.set_ylim(0, max_val2)
    
    ax2.set_xlabel('Birth', fontsize=12)
    ax2.set_ylabel('Death', fontsize=12)
    ax2.set_title(f'{labels[1]} - $H_{dimension}$', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal', adjustable='box')
    
    if title:
        fig.suptitle(title, fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    if show:
        plt.show()
    
    return fig


def plot_betti_comparison(betti_curves1: Dict[int, np.ndarray],
                         betti_curves2: Dict[int, np.ndarray],
                         dimension: int,
                         filtration_values: np.ndarray,
                         labels: Tuple[str, str] = ('State 1', 'State 2'),
                         title: Optional[str] = None,
                         show: bool = True) -> plt.Figure:
    """
    Compare Betti curves for two states.
    
    Args:
        betti_curves1: Dictionary mapping dimension -> Betti values for state 1
        betti_curves2: Dictionary mapping dimension -> Betti values for state 2
        dimension: Homology dimension to compare
        filtration_values: Array of filtration values
        labels: Labels for the two states
        title: Plot title
        show: Whether to display the plot
    
    Returns:
        Matplotlib figure object
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    if dimension in betti_curves1:
        ax.plot(filtration_values, betti_curves1[dimension], 
               'b-', linewidth=2, label=labels[0], alpha=0.7)
    
    if dimension in betti_curves2:
        ax.plot(filtration_values, betti_curves2[dimension], 
               'r--', linewidth=2, label=labels[1], alpha=0.7)
    
    ax.set_ylabel(f'Betti Number $\\beta_{dimension}$', fontsize=12)
    ax.set_xlabel('Filtration Value (ε)', fontsize=12)
    
    if title is None:
        title = f'Betti Number Comparison for Dimension {dimension}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if show:
        plt.show()
    
    return fig

