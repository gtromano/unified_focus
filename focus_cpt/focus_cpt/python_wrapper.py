"""
focus_cpt: Python bindings for the FOCuS changepoint detection library.

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
    """Online (sequential) changepoint detector.

    A thin, Pythonic wrapper around the compiled FOCuS C++ detector. Use this
    class to run the online (sequential) detector: call ``update(y)`` for new
    observation(s) and ``get_statistics(...)`` to compute the current test
    statistic and (optional) detection result.

    Parameters
    ----------
    type : str
        Detector type; one of:
        - ``"univariate"``: two-sided univariate detection
        - ``"univariate_one_sided"``: one-sided univariate detection
        - ``"multivariate"``: multivariate detection using projection sets
        - ``"npfocus"``: nonparametric NP-FOCuS detector (requires ``quantiles``)
        - ``"arp"``: AutoRegressive Process detection (requires ``rho``)
    dim_indexes : list of array-like, optional
        For ``type='multivariate'``, a list of index arrays (0-based) giving
        which dimensions to use for each projection. See
        :func:`generate_projection_indexes` for a helper to build these.
    quantiles : array-like, optional
        For ``type='npfocus'`` this numeric vector of quantiles is required.
    pruning_mult : int, default 2
        Candidate pruning multiplier (controls how aggressively candidates are pruned).
    pruning_offset : int, default 1
        Candidate pruning offset.
    side : str, default 'right'
        For one-sided univariate detectors: ``'right'`` (detect increases) or
        ``'left'`` (detect decreases).
    anomaly_intensity : float, optional
        Optional anomaly intensity threshold used to discard candidate
        segments with too-small signal magnitude. ``None`` (default) disables
        this behaviour.
    rho : array-like, optional
        For ``type='arp'``, this numeric vector of AR coefficients (lag-1, lag-2, ..., lag-p)
        is required.
    mu0_arp : float, optional
        For ``type='arp'``, optional pre-change mean parameter. When provided,
        enables more efficient pruning by filtering candidates based on the known
        pre-change parameter. ``None`` (default) disables this behaviour.

    Notes
    -----
    Multivariate detection with many dimensions can be approximated by using
    projections (subsets of dimensions) supplied via ``dim_indexes``.
    
    
    
    AutoRegressive Process (ARP):
    When ``type = "arp"``, the ``rho`` parameter must be provided as a numeric
    vector of AR coefficients. The detector then computes statistics optimal for
    detecting changepoints in AR(p) processes. Use ``get_statistics(family = "arp")``
    to retrieve the test statistics. The optional ``mu0_arp`` parameter specifies
    the pre-change mean and is tied to the pruning logic: if provided, it enables
    more efficient pruning by allowing the algorithm to filter candidates based on
    the known pre-change parameter.

    Examples
    --------
    >>> # Simple univariate online usage
    >>> det = Detector(type='univariate')
    >>> for y in [0.1, 0.2, 2.5]:
    ...     det.update(y)
    ...     r = det.get_statistics(family='gaussian')
    ...     if r['stat'] is not None and r['stat'] > 20:
    ...         print('Detection at', r['stopping_time'], 'cp=', r['changepoint'])

    >>> # Multivariate detector with projection sets
    >>> proj = generate_projection_indexes(D=3, p=2)
    >>> det_mv = Detector(type='multivariate', dim_indexes=proj)
    >>> det_mv.update([0.5, 1.2, -0.3])

    >>> # Nonparametric NP-FOCuS (quantiles required)
    >>> import numpy as np
    >>> quants = np.percentile(np.random.randn(1000), [25, 50, 75])
    >>> det_np = Detector(type='npfocus', quantiles=quants)
    
    >>> # AutoRegressive Process (ARP) detector
    >>> rho = np.array([0.7, 0.2])  # AR(2) coefficients
    >>> det_arp = Detector(type='arp', rho=rho)
    """

    def __init__(
        self,
        type: str,
        dim_indexes: Optional[List[Union[List[int], np.ndarray]]] = None,
        quantiles: Optional[Union[List[float], np.ndarray]] = None,
        pruning_mult: int = 2,
        pruning_offset: int = 1,
        side: str = "right",
        anomaly_intensity: Optional[float] = None,
        rho: Optional[Union[List[float], np.ndarray]] = None,
        mu0_arp: Optional[float] = None,
    ):
        if dim_indexes is not None:
            dim_indexes = [np.asarray(idx, dtype=np.int32) for idx in dim_indexes]
        if quantiles is not None:
            quantiles = np.asarray(quantiles, dtype=np.float64)
        if anomaly_intensity is not None:
            anomaly_intensity = np.array([anomaly_intensity], dtype=np.float64)
        if rho is not None:
            rho = np.asarray(rho, dtype=np.float64)
        
        # ARP detector requires mu0_arp to be specified
        if type == "arp" and mu0_arp is None:
            raise ValueError(
                "The mu0_arp unknown case is still under development, "
                "please specify a parameter `mu0_arp` when using type='arp'."
            )
        
        # if mu0_arp is not None:
        #     mu0_arp = np.array([mu0_arp], dtype=np.float64)

        self._detector = _focus.Detector(
            type=type,
            dim_indexes=dim_indexes,
            quantiles=quantiles,
            pruning_mult=pruning_mult,
            pruning_offset=pruning_offset,
            side=side,
            anomaly_intensity=anomaly_intensity,
            rho=rho,
            mu0_arp=mu0_arp,
        )
        self._type = type

    def update(self, y: Union[float, List[float], np.ndarray], lambda_: float = 1.0) -> None:
        """
        Update the detector with new observation(s).

        Parameters
        ----------
        y : float, list, or array
            New observation(s). Must match the detector's dimensionality.
        lambda_ : float, default 1.0
            Observation weight (scaling factor for the observation count increment).
            Allows for non-fixed background rates by weighting observations. When
            ``lambda_ = 1.0`` (default), behaves as standard CUSUM with unit
            increments. For ``lambda_ < 1.0``, observations are weighted less;
            for ``lambda_ > 1.0``, observations are weighted more heavily.
            Useful for modeling variable background rates or importance weighting
            of observations.
        """
        y = np.ascontiguousarray(y, dtype=np.float64)
        self._detector.update(y, lambda_)

    def get_statistics(
        self,
        family: str,
        theta0: Optional[Union[float, List[float], np.ndarray]] = None,
        shape: Optional[float] = None,
    ) -> "DetectorStatistics":
        """Compute the current changepoint statistic.

        Parameters
        ----------
        family : str
            Distribution family. Supported values:
            - ``'gaussian'``: Gaussian (normal) log-likelihood ratio
            - ``'poisson'``: Poisson cost
            - ``'bernoulli'``: Bernoulli (binary) cost
            - ``'gamma'``: Gamma cost (requires positive ``shape``)
            - ``'npfocus'``: Nonparametric NP-FOCuS (only valid when the detector
              was created with ``type='npfocus'`` and ``quantiles`` provided)
            - ``'arp'``: AutoRegressive Process cost (only valid when the detector
              was created with ``type='arp'`` and ``rho`` provided)
        theta0 : float, list, or array, optional
            Null-hypothesis parameter(s). For univariate detectors supply a
            scalar; for multivariate detectors supply a vector matching the
            number of dimensions. Note: for ``family='arp'``, the pre-change
            mean should be specified at detector creation time via ``mu0_arp``,
            not here.
        shape : float, optional
            Shape parameter for the gamma family. Must be provided and
            positive when ``family=='gamma'``; otherwise it is ignored.

        Returns
        -------
        DetectorStatistics
            A dictionary, with a compact printed representation, with keys:
            - ``'stopping_time'``: float, number of observations processed
            - ``'changepoint'``: float or ``None``, estimated changepoint location
              (index of the last pre-change observation), or ``None``
            - ``'stat'``: float, numpy array or ``None``, the test statistic(s)
            The values are also available as attributes, e.g. ``result.stat``.

        Notes
        -----
        - When ``family=='gamma'`` a positive ``shape`` must be supplied.
        - For ``family=='npfocus'`` the detector must have been created with
          ``type='npfocus'`` and given ``quantiles``; NP-FOCuS returns two
          statistics (sum and max across quantiles) for which the returned
          ``'stat'`` can be an array.
        - For ``family=='arp'`` the detector must have been created with
          ``type='arp'`` and given ``rho``; the ARP cost function computes
          the test statistic optimal for AR(p) processes. The ``theta0``
          parameter is ignored for ARP; use ``mu0_arp`` at creation time instead.
        """
        if theta0 is not None:
            theta0 = np.atleast_1d(np.asarray(theta0, dtype=np.float64))
        if shape is not None:
            shape = np.array([shape], dtype=np.float64)

        result = self._detector.get_statistics(family=family, theta0=theta0, shape=shape)

        # Convert singleton stats to scalar
        stat = result.get("stat")
        if isinstance(stat, np.ndarray) and stat.size == 1:
            result["stat"] = float(stat.item())
        return DetectorStatistics(result, family=family)

    def get_n_candidates(self) -> int:
        """Return the number of candidate segments currently tracked."""
        return self._detector.get_n_candidates()

    def get_n(self) -> int:
        """Return the number of observations processed."""
        return self._detector.get_n()

    def get_sn(self) -> np.ndarray:
        """Return the current cumulative-sum statistic.

        For univariate detectors this will be a length-1 array; for multivariate
        detectors it returns an array of length equal to the number of
        dimensions.
        """
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


def _format_number(x: Any, missing: str = "not available") -> str:
    """Format a scalar compactly: integers without decimals, else 4 significant digits."""
    if x is None:
        return missing
    x = float(x)
    if np.isfinite(x) and x.is_integer():
        return str(int(x))
    return f"{x:.4g}"


def _field(label: str, value: Any) -> str:
    """One ``label: value`` line, with the values aligned."""
    return f"  {label + ':':<15}{value}"


def _stat_names(n_stats: int, family: Optional[str]) -> List[str]:
    if family == "npfocus" and n_stats == 2:
        return ["sum", "max"]
    if n_stats == 1:
        return ["stat"]
    return [f"stat{i + 1}" for i in range(n_stats)]


def _family_label(family: str, shape: Optional[float]) -> str:
    if family == "gamma" and shape is not None:
        return f"gamma (shape = {_format_number(shape)})"
    return family


def _detection_label(detection_time: Optional[float], changepoint: Optional[float]) -> str:
    if detection_time is None:
        return "none"
    label = f"at time {_format_number(detection_time)}"
    if changepoint is not None:
        label += f" (changepoint estimate: {_format_number(changepoint)})"
    return label


def _item(key: str, doc: str) -> property:
    """Read-only attribute access to a dictionary item."""
    return property(lambda self: self[key], doc=doc)


class DetectorStatistics(dict):
    """Changepoint statistics returned by :meth:`Detector.get_statistics`.

    A dictionary with keys ``'stopping_time'``, ``'changepoint'`` and
    ``'stat'``, whose values are also available as attributes (e.g.
    ``result.stat``). The family used to compute the statistics is stored in
    ``family``. Printing the object gives a compact summary, as the ``print``
    method of the R class ``focus_statistics``.
    """

    stopping_time = _item("stopping_time", "Current time index (number of observations processed).")
    changepoint = _item("changepoint", "Estimated changepoint location, or ``None``.")
    stat = _item("stat", "Test statistic(s): float, numpy array (``npfocus``) or ``None``.")

    def __init__(self, values=(), family: Optional[str] = None):
        super().__init__(values)
        self.family = family

    def __repr__(self) -> str:
        header = "focus statistics"
        if self.family is not None:
            header += f" (family: {self.family})"
        if self.stat is None or np.ndim(self.stat) == 0:
            label, value = "statistic", _format_number(self.stat)
        else:
            values = np.ravel(self.stat)
            names = _stat_names(values.size, self.family)
            label = "statistics"
            value = ", ".join(f"{k} = {_format_number(v)}" for k, v in zip(names, values))
        return "\n".join([
            header,
            _field("stopping time", _format_number(self.stopping_time)),
            _field("changepoint", _format_number(self.changepoint)),
            _field(label, value),
        ])


class OfflineResult(dict):
    """Result of :func:`focus_offline`.

    A dictionary with keys ``'stat'``, ``'changepoint'``, ``'detection_time'``,
    ``'detected_changepoint'``, ``'candidates'``, ``'threshold'``, ``'n'``,
    ``'type'``, ``'family'`` and ``'shape'`` (see :func:`focus_offline`), whose
    values are also available as attributes. Printing the object describes the
    detection, :meth:`summary` summarises the test statistics and :meth:`plot`
    draws their trace, as the ``print``, ``summary`` and ``plot`` methods of the
    R class ``focus_offline``.
    """

    stat = _item("stat", "Test statistics over time, array of shape (n, n_stats).")
    changepoint = _item("changepoint", "Changepoint estimate at each time (``None`` if not available).")
    detection_time = _item("detection_time", "Time of the first detection, or ``None``.")
    detected_changepoint = _item("detected_changepoint", "Changepoint estimate at detection, or ``None``.")
    candidates = _item("candidates", "Candidate segments at the end of the run.")
    threshold = _item("threshold", "Threshold(s) used for detection.")
    n = _item("n", "Number of observations processed.")
    type = _item("type", "Detector type.")
    family = _item("family", "Distribution family.")
    shape = _item("shape", "Shape parameter (gamma family), or ``None``.")

    def __repr__(self) -> str:
        return "\n".join([
            "focus offline detection",
            _field("detector type", self.type),
            _field("family", _family_label(self.family, self.shape)),
            _field("observations", _format_number(self.n)),
            _field("threshold", ", ".join(_format_number(t) for t in np.ravel(self.threshold))),
            _field("detection", _detection_label(self.detection_time, self.detected_changepoint)),
        ])

    def summary(self) -> "OfflineSummary":
        """Summarise the test statistics.

        Returns
        -------
        OfflineSummary
            A dictionary with the detector type, family, shape, number of
            observations, detection time, detected changepoint and number of
            final candidates (``'n_candidates'``), and ``'statistics'``: a list
            with, for each test statistic, its maximum over time, the (1-based)
            time of the maximum, the changepoint estimate at that time and the
            threshold.
        """
        stat = np.asarray(self.stat, dtype=np.float64)
        if stat.ndim == 1:
            stat = stat[:, None]
        n_obs, n_stats = stat.shape
        threshold = np.ravel(self.threshold)
        if threshold.size == 1:
            threshold = np.repeat(threshold, n_stats)
        rows = []
        for j, name in enumerate(_stat_names(n_stats, self.family)):
            row = {"statistic": name, "max": None, "time_of_max": None,
                   "changepoint_at_max": None, "threshold": float(threshold[j])}
            if n_obs > 0:
                t = int(np.argmax(stat[:, j]))
                cp = self.changepoint[t]
                row.update(max=float(stat[t, j]), time_of_max=t + 1,
                           changepoint_at_max=None if cp is None else float(cp))
            rows.append(row)
        return OfflineSummary({
            "type": self.type,
            "family": self.family,
            "shape": self.shape,
            "n": self.n,
            "detection_time": self.detection_time,
            "detected_changepoint": self.detected_changepoint,
            "n_candidates": len(self.candidates["tau"]),
            "statistics": rows,
        })

    def plot(self, ax=None, **kwargs):
        """Plot the trace of the test statistic(s) over time.

        Finite thresholds are drawn as dashed horizontal lines and, if a
        detection occurred, the detection time and the estimated changepoint as
        dotted vertical lines. Requires ``matplotlib``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on. By default, a new figure is created.
        **kwargs
            Passed to ``ax.plot`` for the statistic traces.

        Returns
        -------
        matplotlib.axes.Axes
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError as err:
            raise ImportError("OfflineResult.plot() requires matplotlib "
                              "(pip install matplotlib)") from err
        if ax is None:
            _, ax = plt.subplots()
        stat = np.asarray(self.stat, dtype=np.float64)
        if stat.ndim == 1:
            stat = stat[:, None]
        lines = ax.plot(np.arange(1, stat.shape[0] + 1), stat, **kwargs)
        threshold = np.ravel(self.threshold)
        if threshold.size == 1:
            if np.isfinite(threshold[0]):
                ax.axhline(threshold[0], linestyle="--", color="black")
        else:
            for line, thr in zip(lines, threshold):
                if np.isfinite(thr):
                    ax.axhline(thr, linestyle="--", color=line.get_color())
        if self.detection_time is not None:
            ax.axvline(self.detection_time, linestyle=":", color="black")
        if self.detected_changepoint is not None:
            ax.axvline(self.detected_changepoint, linestyle=":", color="grey")
        if len(lines) > 1:
            for line, name in zip(lines, _stat_names(len(lines), self.family)):
                line.set_label(name)
            ax.legend(frameon=False)
        ax.set_xlabel("Time")
        ax.set_ylabel("Statistic")
        return ax


class OfflineSummary(dict):
    """Summary of an :class:`OfflineResult`, as returned by :meth:`OfflineResult.summary`."""

    def __repr__(self) -> str:
        lines = [
            "focus offline detection: summary",
            _field("detector type", self["type"]),
            _field("family", _family_label(self["family"], self["shape"])),
            _field("observations", _format_number(self["n"])),
            _field("detection", _detection_label(self["detection_time"], self["detected_changepoint"])),
            _field("candidates", self["n_candidates"]),
            "",
            "Statistics:",
        ]
        header = ["", "max", "time of max", "changepoint at max", "threshold"]
        table = [header] + [
            [row["statistic"]] + [_format_number(row[key], "None")
                                  for key in ("max", "time_of_max", "changepoint_at_max", "threshold")]
            for row in self["statistics"]
        ]
        widths = [max(len(row[i]) for row in table) for i in range(len(header))]
        for row in table:
            cells = [row[0].ljust(widths[0])] + [c.rjust(w) for c, w in zip(row[1:], widths[1:])]
            lines.append(" ".join(cells))
        return "\n".join(lines)


def generate_projection_indexes(d: int, p: int) -> List[np.ndarray]:
    """
    Generate projection index sets for high-dimensional multivariate detectors.

    Parameters
    ----------
    d : int
        Number of total dimensions.
    p : int
        Projection subset size.

    Returns
    -------
    list of arrays
        Circular combinations of indices.
    
    Examples
    --------
    >>> # 2-dim projections from 5 dimensions
    >>> generate_projection_indexes(d=5, p=2)
    [array([0,1]), array([1,2]), array([2,3]), array([3,4]), array([4,0])]
    """
    return _focus.generate_projection_indexes(d, p)


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
    anomaly_intensity: Optional[float] = None,
    rho: Optional[Union[List[float], np.ndarray]] = None,
    mu0_arp: Optional[float] = None,
) -> OfflineResult:
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
        Detector type: 'multivariate', 'univariate', 'univariate_one_sided', 'npfocus', 'arp'.
    family : str
        Distribution family: 'gaussian', 'poisson', 'bernoulli', 'gamma', 'npfocus', 'arp'.
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
    anomaly_intensity : float, optional
        Anomaly intensity threshold for pruning candidates. Only candidates with
        sufficient signal magnitude are retained. Default is None (disabled).
    rho : array-like, optional
        For arp detector: AR coefficients (lag-1, lag-2, ..., lag-p).
    mu0_arp : float, optional
        For arp detector: Pre-change mean parameter. When provided,
        enables more efficient pruning based on the known pre-change parameter.
    
    Notes
    -----
    - For multivariate data the detector computes a statistic per projection.
        Detection occurs when any statistic exceeds its corresponding threshold.
    - If ``quantiles`` are provided but ``type!='npfocus'``, a ``ValueError``
        is raised.
    - If ``rho`` is provided but ``type!='arp'``, a ``ValueError`` is raised.
    
    Returns
    -------
    OfflineResult
        Detection results: a dictionary (with a compact printed representation
        and ``summary()`` and ``plot()`` methods) with keys:
        - 'stat': array, test statistics over time (n_obs x n_stats)
        - 'changepoint': array, detected changepoints at each time (or None for None)
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
    
    >>> # ARP detection
    >>> Y = np.concatenate([np.random.randn(100), np.random.randn(100) + 1.0])
    >>> rho = np.array([0.7, 0.2])  # AR(2) coefficients
    >>> result = focus_offline(Y, threshold=10.0, type='arp', family='arp', rho=rho)
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
    if anomaly_intensity is not None:
        anomaly_intensity = np.array([anomaly_intensity], dtype=np.float64)
    if rho is not None:
        rho = np.asarray(rho, dtype=np.float64)
    
    

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
        anomaly_intensity=anomaly_intensity,
        rho=rho,
        mu0_arp=mu0_arp,
    )

    cp = result["changepoint"]
    result["changepoint"] = np.where(cp == -1, None, cp)
    return OfflineResult(result)


__all__ = [
    "Detector",
    "DetectorStatistics",
    "OfflineResult",
    "OfflineSummary",
    "generate_projection_indexes",
    "focus_offline",
]
