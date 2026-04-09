#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

namespace tucvrp {

// Shared project-wide random number generator utilities.
class Rng {
  public:
    // Return the shared engine used across the project.
    static std::mt19937_64& engine();
    // Reseed the shared engine deterministically.
    static void seed(std::uint64_t value);
    // Reseed the shared engine from std::random_device and return the chosen seed.
    static std::uint64_t seed_from_device();
    // Return the most recent seed value used for the shared engine.
    static std::uint64_t current_seed();

    // Draw a uniformly random integer in [low, high].
    static int uniform_int(int low, int high);
    // Draw a uniformly random real in [low, high].
    static double uniform_real(double low, double high);
    // Draw a Bernoulli trial with success probability p.
    static bool bernoulli(double p);

    // Shuffle a vector in place using the shared engine.
    template <typename T>
    static void shuffle(std::vector<T>& values) {
        std::shuffle(values.begin(), values.end(), engine());
    }
};

}  // namespace tucvrp
