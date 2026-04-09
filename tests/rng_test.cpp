#include "tucvrp/rng.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <stdexcept>
#include <vector>

using tucvrp::Rng;

TEST_CASE("rng seeding is reproducible") {
    Rng::seed(12345);
    const int first_int = Rng::uniform_int(1, 1000);
    const double first_real = Rng::uniform_real(0.0, 1.0);
    const bool first_bool = Rng::bernoulli(0.25);

    Rng::seed(12345);
    REQUIRE(Rng::uniform_int(1, 1000) == first_int);
    REQUIRE(Rng::uniform_real(0.0, 1.0) == Catch::Approx(first_real));
    REQUIRE(Rng::bernoulli(0.25) == first_bool);
}

TEST_CASE("rng shuffle is reproducible under fixed seed") {
    std::vector<int> values_a{1, 2, 3, 4, 5, 6};
    std::vector<int> values_b{1, 2, 3, 4, 5, 6};

    Rng::seed(7);
    Rng::shuffle(values_a);
    Rng::seed(7);
    Rng::shuffle(values_b);

    REQUIRE(values_a == values_b);
}

TEST_CASE("rng validates parameter ranges") {
    REQUIRE_THROWS_AS(Rng::uniform_int(5, 4), std::invalid_argument);
    REQUIRE_THROWS_AS(Rng::uniform_real(2.0, 1.0), std::invalid_argument);
    REQUIRE_THROWS_AS(Rng::bernoulli(-0.1), std::invalid_argument);
    REQUIRE_THROWS_AS(Rng::bernoulli(1.1), std::invalid_argument);
}
