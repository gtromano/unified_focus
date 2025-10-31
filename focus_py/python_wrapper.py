"""
Python wrapper for FOCuS changepoint detection library.

This provides a more Pythonic interface to the C++ pybind11 extension.
"""

import numpy as np
from typing import Optional, Union, List, Dict, Any
import focus as _focus


class Detector:
    """
    Changepoint detector wrapper.
    
    This class provides a Pythonic interface to the FOCuS changepoint detection
    methods, automatically handling array conversions and providing better
    type hints and documentation.
    
    Parameters
    ----------
    type : str
        Type of detector:
        - 'multivariate': For multivariate data with dimension reduction
        - 'univariate': For univariate data (two-sided)
        - 'univariate_one_sided': For univariate data (one-sided)
        - 'npfocus': For nonparametric detection
    dim_indexes : list of array-like, optional
        For multivariate detectors: list of dimension index arrays for projections.
        Each element should be an array of integers representing dimensions to project onto.
    quantiles : array-like, optional
        For npfocus detectors: quantiles to use for nonparametric detection.
    pruning_mult : int, default=2
        Pruning multiplier for candidate segment management.
    pruning_offset : int, default=1
        Pruning offset for candidate segment management.
    side : str, default='right'
        For one-sided detectors: 'right' for increases, 'left' for decreases.
    
    Examples
    --------
    >>> # Univariate detector
    >>> detector = Detector(type='univariate')
    >>> detector.update([1.0])
    >>> result = detector.get_statistics(family='gaussian', theta0=0.0)
    
    >>> # Multivariate detector with projection
    >>> detector = Detector(type='multivariate', 
    ...                     dim_indexes=[[0, 1], [1, 2]])
    >>> detector.update([1.0, 2.0, 3.0])
    """
    
    def __init__(
        self,
        type: str,
        dim_indexes: Optional[List[Union[List[int], np.ndarray]]] = None,
        quantiles: Optional[Union[List[float], np.ndarray]] = None,
        pruning_mult: int = 2,
        pruning_offset: int = 1,
        side: str = "right"
    ):
        # Convert dim_indexes to list of numpy arrays if provided
        if dim_indexes is not None:
            dim_indexes = [np.asarray(idx, dtype=np.int32) for idx in dim_indexes]
        
        # Convert quantiles to numpy array if provided
        if quantiles is not None:
            quantiles = np.asarray(quantiles, dtype=np.float64)
        
        self._detector = _focus.Detector(
            type=type,
            dim_indexes=dim_indexes,
            quantiles=quantiles,
            pruning_mult=pruning_mult,
            pruning_offset=pruning_offset,
            side=side
        )
        self._type = type
    
    def update(self, y: Union[float, List[float], np.ndarray]) -> None:
        """
        Update detector with new observation(s).
        
        Parameters
        ----------
        y : float, list, or array
            New observation(s). For univariate detectors, can be a scalar or 1D array.
            For multivariate detectors, should be a 1D array matching the dimension.
        """
        y = np.atleast_1d(np.asarray(y, dtype=np.float64))
        self._detector.update(y)
    
    def get_statistics(
        self,
        family: str,
        theta0: Optional[Union[float, List[float], np.ndarray]] = None,
        shape: Optional[float] = None
    ) -> Dict[str, Any]:
        """
        Compute changepoint statistics.
        
        Parameters
        ----------
        family : str
            Distribution family:
            - 'gaussian': Gaussian (normal) distribution
            - 'poisson': Poisson distribution
            - 'bernoulli': Bernoulli (binary) distribution
            - 'gamma': Gamma distribution (requires `shape` parameter)
            - 'npfocus': Nonparametric detection
        theta0 : float, list, or array, optional
            Null hypothesis parameter. For univariate: scalar.
            For multivariate: vector matching dimension.
        shape : float, optional
            Shape parameter for gamma distribution. Required if family='gamma'.
        
        Returns
        -------
        dict
            Dictionary with keys:
            - 'stopping_time': int, current time index
            - 'changepoint': int or None, detected changepoint location
            - 'stat': float or array, test statistic(s)
        """
        if theta0 is not None:
            theta0 = np.atleast_1d(np.asarray(theta0, dtype=np.float64))
        
        if shape is not None:
            shape = np.array([shape], dtype=np.float64)
        
        result = self._detector.get_statistics(
            family=family,
            theta0=theta0,
            shape=shape
        )
        
        # Convert stat to scalar if it's a 0-d or 1-element array
        if isinstance(result['stat'], np.ndarray):
            if result['stat'].size == 1:
                result['stat'] = float(result['stat'].flat[0])
        
        return result
    
    def get_n_candidates(self) -> int:
        """Get number of candidate segments currently tracked."""
        return self._detector.get_n_candidates()
    
    def get_n(self) -> int:
        """Get number of observations processed."""
        return self._detector.get_n()
    
    def get_sn(self) -> np.ndarray:
        """Get cumulative sum statistic."""
        return self._detector.get_sn()
    
    def get_candidates(self) -> Dict[str, Any]:
        """
        Get candidate segments.
        
        Returns
        -------
        dict
            Dictionary with keys:
            - 'tau': array of int, candidate segment boundaries
            - 'st': list of arrays, sufficient statistics for each candidate
            - 'side': list of str, side indicator for each candidate
        """
        return self._detector.get_candidates()
    
    @property
    def type(self) -> str:
        """Get detector type."""
        return self._type
    
    def __repr__(self) -> str:
        return f"Detector(type='{self._type}', n={self.get_n()}, n_candidates={self.get_n_candidates()})"


def generate_projection_indexes(D: int, p: int) -> List[np.ndarray]:
    """
    Generate circular combinations for multivariate projections.
    
    Parameters
    ----------
    D : int
        Total number of dimensions.
    p : int
        Size of each projection subset.
    
    Returns
    -------
    list of arrays
        List of index arrays for projections.
    
    Examples
    --------
    >>> indexes = generate_projection_indexes(D=5, p=2)
    >>> len(indexes)
    5
    """
    return _focus.generate_projection_indexes(D, p)


def focus_offline(
    Y: Union[List[float], np.ndarray],
    threshold: Union[float, List[float], np.ndarray],
    type: str,
    family: str,
    theta0: Optional[Union[float, List[float], np.ndarray]] = None,
    dim_indexes: Optional[List[Union[List[int], np.ndarray]]] = None,
    quantiles: Optional[Union[List[float], np.ndarray]] = None,
    pruning_mult: int = 2,
    pruning_offset: int = 1,
    side: str = "right",
    shape: Optional[float] = None
) -> Dict[str, Any]:
    """
    Run complete offline changepoint detection.
    
    This function processes the entire dataset at once and returns detection results.
    
    Parameters
    ----------
    Y : array-like
        Data array. Can be 1D (univariate) or 2D (multivariate, observations x dimensions).
    threshold : float or array-like
        Detection threshold(s). Can be a scalar (applied to all statistics) or
        an array matching the number of statistics.
    type : str
        Detector type: 'multivariate', 'univariate', 'univariate_one_sided', 'npfocus'.
    family : str
        Distribution family: 'gaussian', 'poisson', 'bernoulli', 'gamma', 'npfocus'.
    theta0 : float or array-like, optional
        Null hypothesis parameter.
    dim_indexes : list of array-like, optional
        For multivariate: dimension indexes for projections.
    quantiles : array-like, optional
        For npfocus: quantiles to use.
    pruning_mult : int, default=2
        Pruning multiplier.
    pruning_offset : int, default=1
        Pruning offset.
    side : str, default='right'
        For one-sided: 'right' or 'left'.
    shape : float, optional
        Shape parameter for gamma distribution.
    
    Returns
    -------
    dict
        Detection results with keys:
        - 'stat': array, test statistics over time (n_obs x n_stats)
        - 'changepoint': array, detected changepoints at each time (or -1 for None)
        - 'detection_time': int or None, time of first detection
        - 'detected_changepoint': int or None, changepoint at detection time
        - 'candidates': dict, candidate segments
        - 'threshold': array, threshold(s) used
        - 'n': int, number of observations processed
        - 'type': str, detector type
        - 'family': str, distribution family
        - 'shape': float or None, shape parameter (for gamma)
    
    Examples
    --------
    >>> # Univariate Gaussian detection
    >>> Y = np.concatenate([np.random.randn(100), np.random.randn(100) + 2])
    >>> result = focus_offline(Y, threshold=10.0, type='univariate', family='gaussian')
    >>> print(f"Detection at time: {result['detection_time']}")
    >>> print(f"Changepoint at: {result['detected_changepoint']}")
    
    >>> # Multivariate detection
    >>> Y = np.random.randn(200, 3)
    >>> Y[100:, :] += 1.0  # Add mean shift
    >>> result = focus_offline(Y, threshold=15.0, type='multivariate', family='gaussian')
    """
    # Convert inputs to numpy arrays
    Y = np.asarray(Y, dtype=np.float64)
    threshold = np.atleast_1d(np.asarray(threshold, dtype=np.float64))
    
    if theta0 is not None:
        theta0 = np.atleast_1d(np.asarray(theta0, dtype=np.float64))
    
    if dim_indexes is not None:
        dim_indexes = [np.asarray(idx, dtype=np.int32) for idx in dim_indexes]
    
    if quantiles is not None:
        quantiles = np.asarray(quantiles, dtype=np.float64)
    
    if shape is not None:
        shape = np.array([shape], dtype=np.float64)
    
    result = _focus.focus_offline(
        Y=Y,
        threshold=threshold,
        type=type,
        family=family,
        theta0=theta0,
        dim_indexes=dim_indexes,
        quantiles=quantiles,
        pruning_mult=pruning_mult,
        pruning_offset=pruning_offset,
        side=side,
        shape=shape
    )
    
    # Convert -1 to None in changepoint array for better Pythonic interface
    cp_array = result['changepoint']
    result['changepoint'] = np.where(cp_array == -1, None, cp_array)
    
    return result


__all__ = [
    'Detector',
    'generate_projection_indexes',
    'focus_offline',
]
