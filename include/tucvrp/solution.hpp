#pragma once

#include <vector>

namespace tucvrp {

struct Tour {
    std::vector<int> terminals;
    double demand = 0.0;
    double cost = 0.0;
};

struct SolveResult {
    double cost = 0.0;
    std::vector<Tour> tours;
};

}  // namespace tucvrp
