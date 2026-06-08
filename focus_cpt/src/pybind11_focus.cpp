// focus_cptbind11.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "Info.h"
#include "MultivariateInfo.h"
#include "ChangepointResult.h"
#include "Costs.h"
#include "NonparametricInfo.h"
#include "CostsNonparametric.h"
#include "ARpInfo.h"
#include "CostsArp.h"

namespace py = pybind11;
using namespace changepoint;

// Python-facing Detector class that wraps Info
class Detector {
private:
    std::shared_ptr<Info> info_;
    std::string type_;
    
public:
    Detector(const std::string& type,
             const py::object& dim_indexes = py::none(),
             const py::object& quantiles = py::none(),
             int pruning_mult = 2,
             int pruning_offset = 1,
             const std::string& side = "right",
             const py::object& anomaly_intensity = py::none(),
             const py::object& rho = py::none(),
             const py::object& mu0_arp = py::none())
        : type_(type) {
        
        // Warn if quantiles provided but not npfocus
        if (!quantiles.is_none() && type != "npfocus") {
            py::print("Warning: `quantiles` parameter provided but will be ignored unless type == 'npfocus'");
        }
        
        // Warn if rho provided but not arp
        if (!rho.is_none() && type != "arp") {
            py::print("Warning: `rho` parameter provided but will be ignored unless type == 'arp'");
        }
        
        // Parse anomaly_intensity (scalar)
        double anomaly_intensity_val = std::numeric_limits<double>::quiet_NaN();
        if (!anomaly_intensity.is_none()) {
            py::array_t<double> ai = anomaly_intensity.cast<py::array_t<double>>();
            auto ai_buf = ai.request();
            if (ai_buf.size >= 1) {
                anomaly_intensity_val = static_cast<double*>(ai_buf.ptr)[0];
            }
        }
        
        if (type == "multivariate") {
            // Parse dim_indexes
            std::vector<std::vector<int>> dim_idx;
            if (!dim_indexes.is_none()) {
                py::list l = dim_indexes.cast<py::list>();
                for (size_t i = 0; i < l.size(); ++i) {
                    py::array_t<int> arr = l[i].cast<py::array_t<int>>();
                    auto buf = arr.request();
                    if (buf.size == 0) {
                        throw std::runtime_error("Each dim_indexes element must be a non-empty array");
                    }
                    
                    std::vector<int> proj;
                    proj.reserve(buf.size);
                    int* ptr = static_cast<int*>(buf.ptr);
                    for (ssize_t j = 0; j < buf.size; ++j) {
                        proj.push_back(ptr[j]);
                    }
                    dim_idx.push_back(std::move(proj));
                }
            }
            
            info_ = std::make_shared<MultivariateInfo>(
                std::vector<double>{0.0},
                0,
                dim_idx,
                pruning_mult,
                pruning_offset,
                anomaly_intensity_val
            );
            
        } else if (type == "univariate") {
            info_ = std::make_shared<UnivariateInfo>(0.0, 0, anomaly_intensity_val);
            
        } else if (type == "univariate_one_sided") {
            if (side != "right" && side != "left") {
                throw std::invalid_argument("side must be 'right' or 'left'");
            }
            info_ = std::make_shared<OneSideUnivariateInfo>(0.0, 0, side, anomaly_intensity_val);
            
        } else if (type == "npfocus") {
            if (quantiles.is_none()) {
                throw std::invalid_argument("type == 'npfocus' requires a `quantiles` argument");
            }
            
            py::array_t<double> qv = quantiles.cast<py::array_t<double>>();
            auto buf = qv.request();
            if (buf.size < 1) {
                throw std::invalid_argument("`quantiles` must be a non-empty array");
            }
            
            std::vector<double> quants(static_cast<size_t>(buf.size));
            double* ptr = static_cast<double*>(buf.ptr);
            for (ssize_t i = 0; i < buf.size; ++i) {
                quants[i] = ptr[i];
            }
            
            info_ = std::make_shared<NonparametricInfo>(quants);
            
        } else if (type == "arp") {
            // ARP (AutoRegressive Process) detector: requires rho vector
            if (rho.is_none()) {
                throw std::invalid_argument("type == 'arp' requires a `rho` argument (AR coefficients).");
            }
            
            py::array_t<double> rv = rho.cast<py::array_t<double>>();
            auto buf = rv.request();
            if (buf.size < 1) {
                throw std::invalid_argument("`rho` must be a non-empty array (AR order p >= 1)");
            }
            
            std::vector<double> rho_vec(static_cast<size_t>(buf.size));
            double* ptr = static_cast<double*>(buf.ptr);
            for (ssize_t i = 0; i < buf.size; ++i) {
                rho_vec[i] = ptr[i];
            }
            
            // Determine known_prechange based on whether mu0_arp is provided
            bool known_prechange = !mu0_arp.is_none();

            if (known_prechange) {
                // create a info_ with mu0_arp as the known pre-change mean
                info_ = std::make_shared<ARpInfo>(rho_vec, known_prechange, mu0_arp.cast<double>());
            } else {
                // create a info_ without known pre-change mean (mu0_arp will be ignored)
                info_ = std::make_shared<ARpInfo>(rho_vec, known_prechange);
            }
            
        } else {
            throw std::invalid_argument(
                "type must be one of: 'multivariate', 'univariate', 'univariate_one_sided', 'npfocus', 'arp'"
            );
        }
    }
    
    void update(const py::array_t<double>& y, double lambda = 1.0) {
        auto buf = y.request();
        std::vector<double> y_vec(static_cast<size_t>(buf.size));
        double* ptr = static_cast<double*>(buf.ptr);
        for (ssize_t i = 0; i < buf.size; ++i) {
            y_vec[i] = ptr[i];
        }
        info_->update(y_vec, lambda);
    }
    
    py::dict get_statistics(const std::string& family,
                           const py::object& theta0 = py::none(),
                           const py::object& shape = py::none()) {
        
        // Parse theta0
        std::vector<double> theta0_vec;
        if (!theta0.is_none()) {
            py::array_t<double> th = theta0.cast<py::array_t<double>>();
            auto buf = th.request();
            if (buf.size >= 1) {
                theta0_vec.resize(buf.size);
                double* ptr = static_cast<double*>(buf.ptr);
                for (ssize_t i = 0; i < buf.size; ++i) {
                    theta0_vec[i] = ptr[i];
                }
            }
        }
        
        // Parse shape
        double shape_scalar = std::numeric_limits<double>::quiet_NaN();
        if (!shape.is_none()) {
            py::array_t<double> sh = shape.cast<py::array_t<double>>();
            auto buf = sh.request();
            if (buf.size >= 1) {
                shape_scalar = static_cast<double*>(buf.ptr)[0];
            }
        }
        
        // Validate shape parameter
        if (!std::isnan(shape_scalar) && family != "gamma") {
            py::print("Warning: `shape` parameter provided but will be ignored unless family == 'gamma'");
        }
        
        if (family == "gamma") {
            if (std::isnan(shape_scalar) || !(shape_scalar > 0.0)) {
                throw std::invalid_argument("family == 'gamma' requires a positive scalar `shape` parameter");
            }
        }
        
        // Warn if theta0 provided for ARP
        if (!theta0.is_none() && family == "arp") {
            py::print("Warning: For ARP detection, the pre-change mean should be specified at detector creation time (mu0_arp parameter), not as theta0 in get_statistics. The theta0 parameter is ignored for ARP family.");
        }
        
        // Compute costs
        ChangepointResult result;
        
        if (family == "gaussian") {
            result = compute_costs_gaussian(*info_, theta0_vec);
        } else if (family == "poisson") {
            result = compute_costs_poisson(*info_, theta0_vec);
        } else if (family == "bernoulli") {
            result = compute_costs_bernoulli(*info_, theta0_vec);
        } else if (family == "gamma") {
            result = compute_costs_gamma(*info_, theta0_vec, shape_scalar);
        } else if (family == "npfocus") {
            result = compute_costs_npfocus(*info_, theta0_vec);
        } else if (family == "arp") {
            result = compute_costs_arp(*info_, theta0_vec);
        } else {
            throw std::invalid_argument(
                "Unknown family: must be 'gaussian', 'poisson', 'bernoulli', 'gamma', 'npfocus', or 'arp'"
            );
        }
        
        // Convert to Python dict
        py::dict out;
        out["stopping_time"] = result.stopping_time;
        
        if (result.changepoint.has_value()) {
            out["changepoint"] = result.changepoint.value();
        } else {
            out["changepoint"] = py::none();
        }
        
        if (result.stat.has_value()) {
            std::visit([&out](auto&& x) {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, double>) {
                    out["stat"] = x;
                } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                    out["stat"] = py::array_t<double>(x.size(), x.data());
                }
            }, *result.stat);
        } else {
            out["stat"] = py::none();
        }
        
        return out;
    }
    
    int get_n_candidates() const {
        return static_cast<int>(info_->candidates().size());
    }
    
    int get_n() const {
        return info_->n();
    }
    
    py::array_t<double> get_sn() const {
        const auto& sn = info_->sn();
        return py::array_t<double>(sn.size(), sn.data());
    }
    
    py::dict get_candidates() const {
        const auto& candidates = info_->candidates();
        const size_t K = candidates.size();
        
        std::vector<double> tau_vec;
        std::vector<std::string> side_vec;
        py::list st_list;
        
        tau_vec.reserve(K);
        side_vec.reserve(K);
        
        for (size_t i = 0; i < K; ++i) {
            const Candidate& c = candidates[i];
            tau_vec.push_back(c.tau);
            side_vec.push_back(c.side);
            
            if (!c.st.empty()) {
                st_list.append(py::array_t<double>(c.st.size(), c.st.data()));
            } else {
                st_list.append(py::none());
            }
        }
        
        py::dict out;
        out["tau"] = py::array_t<double>(tau_vec.size(), tau_vec.data());
        out["side"] = side_vec;
        out["st"] = st_list;
        
        return out;
    }
    
    std::string get_type() const {
        return type_;
    }
};


// Standalone function for generating projection indexes
py::list generate_projection_indexes(int D, int p) {
    std::vector<std::vector<int>> combs = generate_circular_combinations(D, p);
    py::list result;
    for (const auto& comb : combs) {
        result.append(py::array_t<int>(comb.size(), comb.data()));
    }
    return result;
}


// Complete offline detection
py::dict focus_offline(const py::array_t<double>& Y,
                      const py::object& threshold,
                      const std::string& type,
                      const std::string& family,
                      const py::object& theta0 = py::none(),
                      const py::object& dim_indexes = py::none(),
                      const py::object& quantiles = py::none(),
                      int pruning_mult = 2,
                      int pruning_offset = 1,
                      const std::string& side = "right",
                      const py::object& shape = py::none(),
                      const py::object& anomaly_intensity = py::none(),
                      const py::object& rho = py::none(),
                      const py::object& mu0_arp = py::none()) {
    
    // Parse input data Y
    auto Y_buf = Y.request();
    int n_obs, p_dim;
    bool is_matrix = (Y_buf.ndim == 2);
    
    if (is_matrix) {
        n_obs = static_cast<int>(Y_buf.shape[0]);
        p_dim = static_cast<int>(Y_buf.shape[1]);
    } else {
        n_obs = static_cast<int>(Y_buf.size);
        p_dim = 1;
    }
    
    double* Y_ptr = static_cast<double*>(Y_buf.ptr);
    
    // Parse threshold
    std::vector<double> threshold_vec;
    py::array_t<double> th = threshold.cast<py::array_t<double>>();
    auto th_buf = th.request();
    if (th_buf.size == 0) {
        throw std::invalid_argument("threshold must be a non-empty array or scalar");
    }
    threshold_vec.resize(th_buf.size);
    double* th_ptr = static_cast<double*>(th_buf.ptr);
    for (ssize_t i = 0; i < th_buf.size; ++i) {
        threshold_vec[i] = th_ptr[i];
    }
    
    // Validate type against data dimensions
    if (type == "multivariate" && p_dim == 1) {
        py::print("Warning: type='multivariate' specified but Y is univariate. Consider using type='univariate' or 'univariate_one_sided'");
    }
    if ((type == "univariate" || type == "univariate_one_sided") && p_dim > 1) {
        py::print("Warning: type='" + type + "' specified but Y is multivariate (" + std::to_string(p_dim) + " columns). Consider using type='multivariate'");
    }
    
    // Validate dim_indexes if provided
    if (!dim_indexes.is_none() && type == "multivariate") {
        py::list l = dim_indexes.cast<py::list>();
        for (size_t i = 0; i < l.size(); ++i) {
            py::array_t<int> iv = l[i].cast<py::array_t<int>>();
            auto iv_buf = iv.request();
            int* iv_ptr = static_cast<int*>(iv_buf.ptr);
            for (ssize_t j = 0; j < iv_buf.size; ++j) {
                int idx = iv_ptr[j];
                if (idx < 0 || idx >= p_dim) {
                    throw std::out_of_range("dim_indexes contains index out of range [0, " + 
                                          std::to_string(p_dim) + ") for column count " + 
                                          std::to_string(p_dim));
                }
            }
        }
    }
    
    // Warn if quantiles provided for non-npfocus
    if (!quantiles.is_none() && type != "npfocus") {
        throw std::invalid_argument("`quantiles` parameter provided but will be ignored unless type == 'npfocus'");
    }

    // Warn if rho provided for non-arp
    if (!rho.is_none() && type != "arp") {
        throw std::invalid_argument("`rho` parameter provided but will be ignored unless type == 'arp'");
    }
    
    // Create detector
    Detector detector(type, dim_indexes, quantiles, pruning_mult, pruning_offset, side, anomaly_intensity, rho, mu0_arp);
    
    // Parse theta0
    std::vector<double> theta0_vec;
    if (!theta0.is_none()) {
        py::array_t<double> th0 = theta0.cast<py::array_t<double>>();
        auto th0_buf = th0.request();
        if (th0_buf.size >= 1) {
            theta0_vec.resize(th0_buf.size);
            double* th0_ptr = static_cast<double*>(th0_buf.ptr);
            for (ssize_t i = 0; i < th0_buf.size; ++i) {
                theta0_vec[i] = th0_ptr[i];
            }
        }
    }
    
    // Parse shape
    double shape_scalar = std::numeric_limits<double>::quiet_NaN();
    if (!shape.is_none()) {
        py::array_t<double> sh = shape.cast<py::array_t<double>>();
        auto sh_buf = sh.request();
        if (sh_buf.size >= 1) {
            shape_scalar = static_cast<double*>(sh_buf.ptr)[0];
        }
    }
    
    // Validate shape
    if (!std::isnan(shape_scalar) && family != "gamma") {
        py::print("Warning: `shape` parameter provided but will be ignored unless family == 'gamma'");
    }
    
    if (family == "gamma") {
        if (std::isnan(shape_scalar) || !(shape_scalar > 0.0)) {
            throw std::invalid_argument("family == 'gamma' requires a positive scalar `shape` parameter");
        }
    }
    
    // Select cost function
    CostFunction cost_fn;
    if (family == "gaussian") {
        cost_fn = compute_costs_gaussian;
    } else if (family == "poisson") {
        cost_fn = compute_costs_poisson;
    } else if (family == "bernoulli") {
        cost_fn = compute_costs_bernoulli;
    } else if (family == "gamma") {
        cost_fn = [shape_scalar](const Info& cs_inner, const std::vector<double>& th0) -> ChangepointResult {
            return compute_costs_gamma(cs_inner, th0, shape_scalar);
        };
    } else if (family == "npfocus") {
        cost_fn = compute_costs_npfocus;
    } else if (family == "arp") {
        cost_fn = compute_costs_arp;
    } else {
        throw std::invalid_argument("Unknown family: must be 'gaussian', 'poisson', 'bernoulli', 'gamma', 'npfocus' or 'arp'");
    }
    
    // Prepare containers
    std::vector<std::vector<double>> stats_per_time;
    std::vector<int> changepoints;
    stats_per_time.reserve(n_obs);
    changepoints.reserve(n_obs);
    
    int n_stats = -1;
    bool threshold_warning_issued = false;
    
    std::vector<double> y_t(p_dim);
    
    int detection_time = -1;  // Python uses -1 for None
    int detected_changepoint = -1;
    int actual_length = 0;
    
    // Run online detection
    for (int t = 0; t < n_obs; ++t) {
        // Fill y_t
        if (is_matrix) {
            for (int j = 0; j < p_dim; ++j) {
                y_t[j] = Y_ptr[t * p_dim + j];
            }
        } else {
            y_t[0] = Y_ptr[t];
        }
        
        // Update detector
        py::array_t<double> y_array(y_t.size(), y_t.data());
        detector.update(y_array);
        
        // Compute statistics
        py::dict result_dict = detector.get_statistics(family, 
                                                       theta0.is_none() ? py::none() : py::cast(theta0_vec),
                                                       shape.is_none() ? py::none() : py::cast(shape_scalar));
        
        // Extract stats
        std::vector<double> stat_vec;
        py::object stat_obj = result_dict["stat"];
        if (!stat_obj.is_none()) {
            if (py::isinstance<py::array_t<double>>(stat_obj)) {
                py::array_t<double> stat_arr = stat_obj.cast<py::array_t<double>>();
                auto stat_buf = stat_arr.request();
                stat_vec.resize(stat_buf.size);
                double* stat_ptr = static_cast<double*>(stat_buf.ptr);
                for (ssize_t i = 0; i < stat_buf.size; ++i) {
                    stat_vec[i] = stat_ptr[i];
                }
            } else if (py::isinstance<py::float_>(stat_obj) || py::isinstance<py::int_>(stat_obj)) {
                stat_vec.push_back(stat_obj.cast<double>());
            }
        }
        
        // First iteration: determine dimensionality
        if (n_stats == -1) {
            n_stats = stat_vec.empty() ? 1 : static_cast<int>(stat_vec.size());
            
            if (threshold_vec.size() != 1 && static_cast<int>(threshold_vec.size()) != n_stats) {
                throw std::invalid_argument("threshold must be a scalar or array of length " + 
                                          std::to_string(n_stats) + " (matching number of statistics)");
            }
            
            if (threshold_vec.size() == 1 && n_stats > 1 && !threshold_warning_issued) {
                py::print("Warning: Single threshold provided for " + std::to_string(n_stats) + 
                         " statistics. Using threshold = " + std::to_string(threshold_vec[0]) + 
                         " for all statistics");
                threshold_warning_issued = true;
            }
        }
        
        // Ensure consistent dimensionality
        if (static_cast<int>(stat_vec.size()) != n_stats) {
            if (stat_vec.empty()) {
                stat_vec.resize(n_stats, 0.0);
            } else {
                throw std::runtime_error("Inconsistent number of statistics returned across time steps");
            }
        }
        
        stats_per_time.push_back(stat_vec);
        
        py::object cp_obj = result_dict["changepoint"];
        if (!cp_obj.is_none()) {
            changepoints.push_back(cp_obj.cast<int>());
        } else {
            changepoints.push_back(-1);  // Use -1 for None in Python
        }
        
        actual_length = t + 1;
        
        // Check for detection
        bool detected = false;
        for (int s = 0; s < n_stats; ++s) {
            double thresh = (threshold_vec.size() == 1) ? threshold_vec[0] : threshold_vec[s];
            if (stat_vec[s] > thresh) {
                detected = true;
                break;
            }
        }
        
        if (detected) {
            detection_time = t + 1;  // 1-based for consistency with R
            if (!cp_obj.is_none()) {
                detected_changepoint = cp_obj.cast<int>();
            }
            break;
        }
    }
    
    // Convert stats to numpy array
    py::array_t<double> stat_mat({actual_length, n_stats});
    auto stat_mat_buf = stat_mat.request();
    double* stat_mat_ptr = static_cast<double*>(stat_mat_buf.ptr);
    for (int t = 0; t < actual_length; ++t) {
        for (int s = 0; s < n_stats; ++s) {
            stat_mat_ptr[t * n_stats + s] = stats_per_time[t][s];
        }
    }
    
    // Convert changepoints to numpy array (replace -1 with None later in Python wrapper if needed)
    py::array_t<int> cp_array(changepoints.size(), changepoints.data());
    
    // Build result dict
    py::dict result;
    result["stat"] = stat_mat;
    result["changepoint"] = cp_array;
    result["detection_time"] = (detection_time == -1) ? py::none() : py::cast(detection_time);
    result["detected_changepoint"] = (detected_changepoint == -1) ? py::none() : py::cast(detected_changepoint);
    result["candidates"] = detector.get_candidates();
    result["threshold"] = py::array_t<double>(threshold_vec.size(), threshold_vec.data());
    result["n"] = actual_length;
    result["type"] = type;
    result["family"] = family;
    result["shape"] = (family == "gamma") ? py::cast(shape_scalar) : py::none();
    
    return result;
}


PYBIND11_MODULE(_focus, m) {
    m.doc() = "FOCuS (Focusing on Candidate Segments) changepoint detection library";
    
    py::class_<Detector>(m, "Detector")
        .def(py::init<const std::string&, const py::object&, const py::object&, int, int, const std::string&, const py::object&, const py::object&, const py::object&>(),
             py::arg("type"),
             py::arg("dim_indexes") = py::none(),
             py::arg("quantiles") = py::none(),
             py::arg("pruning_mult") = 2,
             py::arg("pruning_offset") = 1,
             py::arg("side") = "right",
             py::arg("anomaly_intensity") = py::none(),
             py::arg("rho") = py::none(),
             py::arg("mu0_arp") = py::none(),
             R"pbdoc(
                Create a new changepoint detector.
                
                Parameters
                ----------
                type : str
                    Type of detector: 'multivariate', 'univariate', 'univariate_one_sided', 'npfocus', or 'arp'
                dim_indexes : list of arrays, optional
                    For multivariate: list of dimension index arrays for projections
                quantiles : array, optional
                    For npfocus: quantiles to use
                pruning_mult : int, default=2
                    Pruning multiplier
                pruning_offset : int, default=1
                    Pruning offset
                side : str, default='right'
                    For one-sided: 'right' or 'left'
                anomaly_intensity : float, optional
                    Anomaly intensity threshold for pruning candidates. Only candidates with
                    sufficient signal magnitude are retained. Default is None (disabled).
                rho : array, optional
                    For arp: AR coefficients (lag-1, lag-2, ..., lag-p)
                mu0_arp : float, optional
                    For arp: Pre-change mean for ARP detector. When provided, enables more efficient
                    pruning based on the known pre-change parameter. Default is None (disabled).
             )pbdoc")
        .def("update", &Detector::update, py::arg("y"), py::arg("lambda") = 1.0,
             "Update detector with new observation(s)")
        .def("get_statistics", &Detector::get_statistics,
             py::arg("family"),
             py::arg("theta0") = py::none(),
             py::arg("shape") = py::none(),
             R"pbdoc(
                Compute changepoint statistics.
                
                Parameters
                ----------
                family : str
                    Distribution family: 'gaussian', 'poisson', 'bernoulli', 'gamma', 'npfocus', or 'arp'
                theta0 : array, optional
                    Null hypothesis parameter (not used for 'arp' family)
                shape : float, optional
                    Shape parameter for gamma distribution
                
                Returns
                -------
                dict
                    Dictionary with 'stopping_time', 'changepoint', and 'stat' keys
             )pbdoc")
        .def("get_n_candidates", &Detector::get_n_candidates,
             "Get number of candidate segments")
        .def("get_n", &Detector::get_n,
             "Get number of observations processed")
        .def("get_sn", &Detector::get_sn,
             "Get cumulative sum statistic")
        .def("get_candidates", &Detector::get_candidates,
             "Get candidate segments as a dictionary")
        .def_property_readonly("type", &Detector::get_type,
             "Get detector type");
    
    m.def("generate_projection_indexes", &generate_projection_indexes,
          py::arg("D"), py::arg("p"),
          "Generate circular combinations for multivariate projections");
    
    m.def("focus_offline", &focus_offline,
          py::arg("Y"),
          py::arg("threshold"),
          py::arg("type"),
          py::arg("family"),
          py::arg("theta0") = py::none(),
          py::arg("dim_indexes") = py::none(),
          py::arg("quantiles") = py::none(),
          py::arg("pruning_mult") = 2,
          py::arg("pruning_offset") = 1,
          py::arg("side") = "right",
          py::arg("shape") = py::none(),
          py::arg("anomaly_intensity") = py::none(),
          py::arg("rho") = py::none(),
          py::arg("mu0_arp") = py::none(),
          R"pbdoc(
            Run complete offline changepoint detection.
            
            Parameters
            ----------
            Y : array
                Data array (1D or 2D)
            threshold : array or float
                Detection threshold(s)
            type : str
                Detector type
            family : str
                Distribution family
            theta0 : array, optional
                Null hypothesis parameter
            dim_indexes : list of arrays, optional
                Dimension indexes for multivariate
            quantiles : array, optional
                Quantiles for npfocus
            pruning_mult : int, default=2
                Pruning multiplier
            pruning_offset : int, default=1
                Pruning offset
            side : str, default='right'
                Side for one-sided detection
            shape : float, optional
                Shape parameter for gamma
            anomaly_intensity : float, optional
                Anomaly intensity threshold for pruning candidates. Only candidates with
                sufficient signal magnitude are retained. Default is None (disabled).
            rho : array, optional
                AR coefficients for arp detector
            mu0_arp : float, optional
                Pre-change mean for arp detector
            
            Returns
            -------
            dict
                Detection results including statistics, changepoints, and candidates
          )pbdoc");
}
