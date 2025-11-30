"""
Protein Flop: Topological Data Analysis for Protein Conformation Transitions
"""

__version__ = "0.1.0"

# Import main modules
from . import data_processing
from . import tda_computation
from . import visualization
from . import analysis
from . import flop_computation

__all__ = [
    'data_processing',
    'tda_computation',
    'visualization',
    'analysis',
    'flop_computation'
]

