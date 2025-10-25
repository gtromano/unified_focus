#pragma once

#include <optional>

namespace changepoint {

struct ChangepointResult {
    int stopping_time;
    std::optional<int> changepoint;
    std::optional<double> stat;
};

} // namespace changepoint