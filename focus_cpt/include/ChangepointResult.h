#pragma once

#include <optional>
#include <variant>

namespace changepoint {

struct ChangepointResult {
    double stopping_time;
    std::optional<double> changepoint;
    std::optional<std::variant<double, std::vector<double>>> stat;
};

} // namespace changepoint
