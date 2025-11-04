"""
focus_py: Python bindings for the FOCuS changepoint detection library.

Provides a clean, Pythonic interface to the high-performance C++ implementation
of the FOCuS and md-FOCuS algorithms for online and offline changepoint detection.

Main interfaces:
- `focus_offline`: full C++ offline detector (fast, batch mode)
- `Detector`: online/sequential interface (step-by-step updates)

See the README or documentation for detailed usage and examples.
"""

import numpy as np
from typing import Optional, Union, List, Dict, Any

# Import compiled backend
try:
    from . import _focus
except ImportError:
    import _focus


class Detector:
    """
    Online (sequential) changepoint detector.

    Provides a step-by-step interface to the FOCuS algorithm.
    Each `update()` call adds new data, and `get_statistics()` computes
    the current test statistic and detection result.

    Parameters
    ----------
    type : str
        One of {"univariate", "univariate_one_sided", "multivariate", "npfocus"}.
    dim_indexes : list of array-like, optional
        Projection index sets for multivariate detectors.
    quantiles : array-like, optional
        Quantiles for nonparametric ("npfocus") detectors.
    pruning_mult, pruning_offset : int, optional
        Candidate pruning parameters (default: 2, 1).
    side : str, optional
        For one-sided detectors, "right" (increases) or "left" (decreases).

    Examples
    --------
    >>> det = Detector(type="univariate")
    >>> for y in [0.1, 0.2, 2.5]:
    ...     det.update(y)
    ...     res = det.get_statistics(family="gaussian")
    ...     if res["stat"] > 20:
    ...         print("Detection:", res)
    """

    def __init__(
        self,
        type: str,
        dim_indexes: Optional[List[Union[List[int], np.ndarray]]] = None,
        quantiles: Optional[Union[List[float], np.ndarray]] = None,
        pruning_mult: int = 2,
        pruning_offset: int = 1,
        side: str = "right",
    ):
        if dim_indexes is not None:
            dim_indexes = [np.asarray(idx, dtype=np.int32) for idx in dim_indexes]
        if quantiles is not None:
            quantiles = np.asarray(quantiles, dtype=np.float64)

        self._detector = _focus.Detector(
            type=type,
            dim_indexes=dim_indexes,
            quantiles=quantiles,
            pruning_mult=pruning_mult,
            pruning_offset=pruning_offset,
            side=side,
        )
        self._type = type

    def update(self, y: Union[float, List[float], np.ndarray]) -> None:
        """
        Update the detector with new observation(s).

        Parameters
        ----------
        y : float, list, or array
            New observation(s). Must match the detector's dimensionality.
        """
        y = np.ascontiguousarray(y, dtype=np.float64)
        self._detector.update(y)

    def get_statistics(
        self,
        family: str,
        theta0: Optional[Union[float, List[float], np.ndarray]] = None,
        shape: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Compute current changepoint statistic.

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

        result = self._detector.get_statistics(family=family, theta0=theta0, shape=shape)

        # Convert singleton stats to scalar
        stat = result.get("stat")
        if isinstance(stat, np.ndarray) and stat.size == 1:
            result["stat"] = float(stat)
        return result

    def get_n_candidates(self) -> int:
        """Return the number of candidate segments currently tracked."""
        return self._detector.get_n_candidates()

    def get_n(self) -> int:
        """Return the number of observations processed."""
        return self._detector.get_n()

    def get_sn(self) -> np.ndarray:
        """Return the current cumulative sum statistic."""
        return self._detector.get_sn()

    def get_candidates(self) -> Dict[str, Any]:
        """
        Return current candidate segments.

        Returns
        -------
        dict
            Dictionary with keys:
            - 'tau': array of int, candidate changepoints
            - 'st': list of arrays, sufficient statistics for each candidate (e.g., cumulative sums of the data)
            - 'side': list of str, side indicator for each candidate
        """
        return self._detector.get_candidates()

    @property
    def type(self) -> str:
        """Return the detector type."""
        return self._type

    def __repr__(self) -> str:
        return f"Detector(type='{self._type}', n={self.get_n()}, n_candidates={self.get_n_candidates()})"


def generate_projection_indexes(D: int, p: int) -> List[np.ndarray]:
    """
    Generate projection index sets for high-dimentional multivariate detectors.

    Parameters
    ----------
    D : int
        Number of total dimensions.
    p : int
        Projection subset size.

    Returns
    -------
    list of arrays
        Circular combinations of indices.
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
    shape: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Run the FOCuS detector in batch/offline mode (entirely in C++).

    Processes all data at once and returns detection results and trajectories.

    Parameters
    ----------
    Y : array-like
        Data array. Can be 1D (univariate) or 2D (multivariate, observations x dimensions).
    threshold : float or array-like
        Detection threshold(s). Can be a scalar (applied to all statistics) or
        an array matching the number of statistics. Use np.inf for no thresholding.
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
    Y = np.ascontiguousarray(Y, dtype=np.float64)
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
        shape=shape,
    )

    cp = result["changepoint"]
    result["changepoint"] = np.where(cp == -1, None, cp)
    return result


__all__ = ["Detector", "generate_projection_indexes", "focus_offline"]
