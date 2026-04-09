#include "tucvrp/rng.hpp"

#include <stdexcept>

namespace tucvrp {

namespace {

std::uint64_t& shared_seed() {
    static std::uint64_t seed = 0;
    return seed;
}

}  // namespace

std::mt19937_64& Rng::engine() {
    static std::mt19937_64 engine(shared_seed());
    return engine;
}

void Rng::seed(std::uint64_t value) {
    shared_seed() = value;
    engine().seed(value);
}

std::uint64_t Rng::seed_from_device() {
    std::random_device device;
    const std::uint64_t value =
        (static_cast<std::uint64_t>(device()) << 32U) ^ static_cast<std::uint64_t>(device());
    seed(value);
    return value;
}

std::uint64_t Rng::current_seed() { return shared_seed(); }

int Rng::uniform_int(int low, int high) {
    if (low > high) {
        throw std::invalid_argument("uniform_int requires low <= high");
    }
    std::uniform_int_distribution<int> distribution(low, high);
    return distribution(engine());
}

double Rng::uniform_real(double low, double high) {
    if (low > high) {
        throw std::invalid_argument("uniform_real requires low <= high");
    }
    std::uniform_real_distribution<double> distribution(low, high);
    return distribution(engine());
}

bool Rng::bernoulli(double p) {
    if (p < 0.0 || p > 1.0) {
        throw std::invalid_argument("bernoulli requires p in [0, 1]");
    }
    std::bernoulli_distribution distribution(p);
    return distribution(engine());
}

}  // namespace tucvrp
