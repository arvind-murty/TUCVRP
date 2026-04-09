#include "tucvrp/instance.hpp"
#include "tucvrp/preprocessing.hpp"
#include "tucvrp/exact_solver.hpp"
#include "tucvrp/rng.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

// Parse one instance file, report a lower bound, and compute the exact optimum for small inputs.
int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            std::cerr << "usage: " << argv[0] << " <instance-file>\n";
            return 1;
        }

        tucvrp::Rng::seed(12345);

        const auto instance = tucvrp::Instance::parse_file(argv[1]);
        const auto exact = tucvrp::ExactSolver::solve(instance);
        const auto lb = tucvrp::Preprocessor::edge_load_lower_bound(instance);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << instance;
        std::cout << "terminals: " << instance.terminal_count() << "\n";
        std::cout << "total demand: " << instance.total_demand() << "\n";
        std::cout << "edge-load lower bound: " << lb << "\n";
        std::cout << "exact optimum: " << exact.cost << "\n";
        std::cout << "tour count: " << exact.tours.size() << "\n";
        for (std::size_t i = 0; i < exact.tours.size(); ++i) {
            const auto& tour = exact.tours[i];
            std::cout << "tour " << (i + 1) << ": demand=" << tour.demand << ", cost=" << tour.cost << ", terminals=";
            for (std::size_t j = 0; j < tour.terminals.size(); ++j) {
                if (j != 0) {
                    std::cout << ",";
                }
                std::cout << tour.terminals[j];
            }
            std::cout << "\n";
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 1;
    }
}
