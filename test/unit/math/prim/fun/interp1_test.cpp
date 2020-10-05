#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <vector>

// This file goes in test/unit/math/prim/fun/
//
// To run the test type:
//  ./runTest.py test/unit/math/prim/fun/interp1_test.cpp

TEST(MathPrim, interp1_test) {
using stan::test::expect_near_rel;
Eigen::VectorXd x(5);
x << 0.00, 1.25, 2.50, 3.75, 5.00;
Eigen::VectorXd y(5);
y << 0.0000000000, 0.9489846194, 0.5984721441, -0.5715613187, -0.9589242747;

Eigen::VectorXd x_test(4);
x_test << 1.0, 1.7, 2.1, 4.0;
Eigen::VectorXd y_ref(4); // Reference values from R
y_ref << 0.7591876955, 0.8228001283, 0.7106361362, -0.6490339099;

Eigen::VectorXd y_test = stan::math::interp1(x, y, x_test);

stan::test::expect_near_rel("test Interp1", y_test, y_ref);
cout << y_test;
}
