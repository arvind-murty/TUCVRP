#include "tucvrp/instance.hpp"
#include "tucvrp/algorithms/one_point_five_approx.hpp"
#include "tucvrp/algorithms/labbe_approx.hpp"
#include "tucvrp/preprocessing.hpp"
#include "tucvrp/exact_solver.hpp"
#include "tucvrp/rng.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

// Parse one instance file, report a lower bound, and compute both exact and approximation tours.
int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            std::cerr << "usage: " << argv[0] << " <instance-file>\n";
            return 1;
        }

        tucvrp::Rng::seed(12345);

        const auto instance = tucvrp::Instance::parse_file(argv[1]);
        const auto approx =
            tucvrp::OnePointFiveApproxSolver::solve(instance, tucvrp::OnePointFiveApproxParams{.epsilon = 0.25});
        const auto labbe = tucvrp::LabbeApproxSolver::solve(instance);
        const auto lb = tucvrp::Preprocessor::edge_load_lower_bound(instance);
        const bool run_exact = instance.terminal_count() <= 20;
        const auto exact = run_exact ? tucvrp::ExactSolver::solve(instance) : tucvrp::SolveResult{};

        std::cout << std::fixed << std::setprecision(6);
        std::cout << instance;
        std::cout << "terminals: " << instance.terminal_count() << "\n";
        std::cout << "total demand: " << instance.total_demand() << "\n";
        std::cout << "edge-load lower bound: " << lb << "\n";
        if (run_exact) {
            std::cout << "exact optimum: " << exact.cost << "\n";
        } else {
            std::cout << "exact optimum: skipped (more than 20 terminals)\n";
        }
        std::cout << "labbe baseline cost: " << labbe.cost << "\n";
        std::cout << "paper solver cost: " << approx.cost << "\n";
        if (run_exact) {
            std::cout << "exact tour count: " << exact.tours.size() << "\n";
            for (std::size_t i = 0; i < exact.tours.size(); ++i) {
                const auto& tour = exact.tours[i];
                std::cout << "exact tour " << (i + 1) << ": demand=" << tour.demand << ", cost=" << tour.cost
                          << ", terminals=";
                for (std::size_t j = 0; j < tour.terminals.size(); ++j) {
                    if (j != 0) {
                        std::cout << ",";
                    }
                    std::cout << tour.terminals[j];
                }
                std::cout << ", walk=";
                for (std::size_t j = 0; j < tour.walk.size(); ++j) {
                    if (j != 0) {
                        std::cout << "->";
                    }
                    std::cout << tour.walk[j];
                }
                std::cout << "\n";
            }
        }
        std::cout << "labbe baseline tour count: " << labbe.tours.size() << "\n";
        for (std::size_t i = 0; i < labbe.tours.size(); ++i) {
            const auto& tour = labbe.tours[i];
            std::cout << "labbe tour " << (i + 1) << ": demand=" << tour.demand << ", cost=" << tour.cost
                      << ", terminals=";
            for (std::size_t j = 0; j < tour.terminals.size(); ++j) {
                if (j != 0) {
                    std::cout << ",";
                }
                std::cout << tour.terminals[j];
            }
            std::cout << ", walk=";
            for (std::size_t j = 0; j < tour.walk.size(); ++j) {
                if (j != 0) {
                    std::cout << "->";
                }
                std::cout << tour.walk[j];
            }
            std::cout << "\n";
        }
        std::cout << "paper tour count: " << approx.tours.size() << "\n";
        for (std::size_t i = 0; i < approx.tours.size(); ++i) {
            const auto& tour = approx.tours[i];
            std::cout << "paper tour " << (i + 1) << ": demand=" << tour.demand << ", cost=" << tour.cost
                      << ", terminals=";
            for (std::size_t j = 0; j < tour.terminals.size(); ++j) {
                if (j != 0) {
                    std::cout << ",";
                }
                std::cout << tour.terminals[j];
            }
            std::cout << ", walk=";
            for (std::size_t j = 0; j < tour.walk.size(); ++j) {
                if (j != 0) {
                    std::cout << "->";
                }
                std::cout << tour.walk[j];
            }
            std::cout << "\n";
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 1;
    }
}
