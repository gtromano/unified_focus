"""
FOCuS

This package provides fast changepoint detection methods for univariate,
multivariate, and nonparametric data using the FOCuS algorithm.
"""

from .python_wrapper import (
    Detector,
    generate_projection_indexes,
    focus_offline,
)

__version__ = "0.1.9"

__all__ = [
    "Detector",
    "generate_projection_indexes",
    "focus_offline",
]
