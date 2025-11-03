#pragma once

#include <optional>
#include <variant>

namespace changepoint {

struct ChangepointResult {
    int stopping_time;
    std::optional<int> changepoint;
    std::optional<std::variant<double, std::vector<double>>> stat;
};

} // namespace changepoint
