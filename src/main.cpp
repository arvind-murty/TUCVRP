#include "tucvrp/instance.hpp"
#include "tucvrp/preprocessing.hpp"
#include "tucvrp/solver.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            std::cerr << "usage: " << argv[0] << " <instance-file>\n";
            return 1;
        }

        const auto instance = tucvrp::Instance::parse_file(argv[1]);
        const auto exact = tucvrp::ExactSolver::solve(instance);
        const auto lb = tucvrp::Preprocessor::edge_load_lower_bound(instance);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "terminals: " << instance.terminal_count() << "\n";
        std::cout << "total demand: " << instance.total_demand() << "\n";
        std::cout << "edge-load lower bound: " << lb << "\n";
        std::cout << "exact optimum: " << exact.cost << "\n";
        std::cout << "tour count: " << exact.tours.size() << "\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 1;
    }
}
