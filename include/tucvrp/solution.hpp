#pragma once

#include <vector>

namespace tucvrp {

// One vehicle tour together with the terminals it serves.
struct Tour {
    std::vector<int> terminals;
    double demand = 0.0;
    double cost = 0.0;
};

// A complete routing solution made of multiple tours.
struct SolveResult {
    double cost = 0.0;
    std::vector<Tour> tours;
};

}  // namespace tucvrp
