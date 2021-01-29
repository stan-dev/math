/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <vector>
#include <random>
#include <boost/math/statistics/linear_regression.hpp>

using boost::math::statistics::simple_ordinary_least_squares;
using boost::math::statistics::simple_ordinary_least_squares_with_R_squared;

template<typename Real>
void test_line()
{
    std::vector<Real> x(128);
    std::vector<Real> y(128);
    Real expected_c0 = 7;
    Real expected_c1 = 12;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = i;
        y[i] = expected_c0 + expected_c1*x[i];
    }

    auto [computed_c0, computed_c1] = simple_ordinary_least_squares(x, y);

    CHECK_ULP_CLOSE(expected_c0, computed_c0, 0);
    CHECK_ULP_CLOSE(expected_c1, computed_c1, 0);

    auto [computed_c0_R, computed_c1_R, Rsquared] = simple_ordinary_least_squares_with_R_squared(x, y);

    Real expected_Rsquared = 1;
    CHECK_ULP_CLOSE(expected_c0, computed_c0, 0);
    CHECK_ULP_CLOSE(expected_c1, computed_c1, 0);
    CHECK_ULP_CLOSE(expected_Rsquared, Rsquared, 0);

}

template<typename Real>
void test_constant()
{
    std::vector<Real> x(128);
    std::vector<Real> y(128);
    Real expected_c0 = 7;
    Real expected_c1 = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = i;
        y[i] = expected_c0 + expected_c1*x[i];
    }

    auto [computed_c0, computed_c1] = simple_ordinary_least_squares(x, y);

    CHECK_ULP_CLOSE(expected_c0, computed_c0, 0);
    CHECK_ULP_CLOSE(expected_c1, computed_c1, 0);

    auto [computed_c0_R, computed_c1_R, Rsquared] = simple_ordinary_least_squares_with_R_squared(x, y);

    Real expected_Rsquared = 1;
    CHECK_ULP_CLOSE(expected_c0, computed_c0, 0);
    CHECK_ULP_CLOSE(expected_c1, computed_c1, 0);
    CHECK_ULP_CLOSE(expected_Rsquared, Rsquared, 0);

}

template<typename Real>
void test_permutation_invariance()
{
    std::vector<Real> x(256);
    std::vector<Real> y(256);
    std::mt19937_64 gen{123456};
    std::normal_distribution<Real> dis(0, 0.1);

    Real expected_c0 = -7.2;
    Real expected_c1 = -13.5;

    x[0] = 0;
    y[0] = expected_c0 + dis(gen);
    for(size_t i = 1; i < x.size(); ++i) {
        Real t = dis(gen);
        x[i] = x[i-1] + t*t;
        y[i] = expected_c0 + expected_c1*x[i] + dis(gen);
    }

    auto [c0, c1, Rsquared] = simple_ordinary_least_squares_with_R_squared(x, y);
    CHECK_MOLLIFIED_CLOSE(expected_c0, c0, 0.002);
    CHECK_MOLLIFIED_CLOSE(expected_c1, c1, 0.002);

    int j = 0;
    std::mt19937_64 gen1{12345};
    std::mt19937_64 gen2{12345};
    while(j++ < 10) {
        std::shuffle(x.begin(), x.end(), gen1);
        std::shuffle(y.begin(), y.end(), gen2);
        auto [c0_, c1_, Rsquared_] = simple_ordinary_least_squares_with_R_squared(x, y);

        CHECK_ULP_CLOSE(c0, c0_, 100);
        CHECK_ULP_CLOSE(c1, c1_, 100);
        CHECK_ULP_CLOSE(Rsquared, Rsquared_, 65);
    }
}

template<typename Real>
void test_scaling_relations()
{
    std::vector<Real> x(256);
    std::vector<Real> y(256);
    std::mt19937_64 gen{123456};
    std::normal_distribution<Real> dis(0, 0.1);

    Real expected_c0 = 3.2;
    Real expected_c1 = -13.5;

    x[0] = 0;
    y[0] = expected_c0 + dis(gen);
    for(size_t i = 1; i < x.size(); ++i) {
        Real t = dis(gen);
        x[i] = x[i-1] + t*t;
        y[i] = expected_c0 + expected_c1*x[i] + dis(gen);
    }

    auto [c0, c1, Rsquared] = simple_ordinary_least_squares_with_R_squared(x, y);
    CHECK_MOLLIFIED_CLOSE(expected_c0, c0, 0.005);
    CHECK_MOLLIFIED_CLOSE(expected_c1, c1, 0.005);

    // If y -> lambda y, then c0 -> lambda c0 and c1 -> lambda c1.
    Real lambda = 6;

    for (auto& s : y) {
        s *= lambda;
    }

    auto [c0_lambda, c1_lambda, Rsquared_lambda] = simple_ordinary_least_squares_with_R_squared(x, y);

    CHECK_ULP_CLOSE(lambda*c0, c0_lambda, 30);
    CHECK_ULP_CLOSE(lambda*c1, c1_lambda, 30);
    CHECK_ULP_CLOSE(Rsquared, Rsquared_lambda, 3);

    // If x -> lambda x, then c0 -> c0 and c1 -> c1/lambda
    for (auto& s : x) {
        s *= lambda;
    }
    // Put y back into it's original state:
    for (auto& s : y) {
        s /= lambda;
    }
    auto [c0_, c1_, Rsquared_] = simple_ordinary_least_squares_with_R_squared(x, y);

    CHECK_ULP_CLOSE(c0, c0_, 50);
    CHECK_ULP_CLOSE(c1, c1_*lambda, 50);
    CHECK_ULP_CLOSE(Rsquared, Rsquared_, 50);

}


int main()
{
    test_line<float>();
    test_line<double>();
    test_line<long double>();

    test_constant<float>();
    test_constant<double>();
    test_constant<long double>();

    test_permutation_invariance<float>();
    test_permutation_invariance<double>();
    test_permutation_invariance<long double>();

    test_scaling_relations<float>();
    test_scaling_relations<double>();
    test_scaling_relations<long double>();
    return boost::math::test::report_errors();
}
