#ifndef STAN_TEST_UNIT_MATH_MIX_PROB_UTIL_HPP
#define STAN_TEST_UNIT_MATH_MIX_PROB_UTIL_HPP

#include <gtest/gtest.h>

class Wiener7MixArgs : public ::testing::Test {
 public:
  const char* function = "function";
  void SetUp() {}

  Eigen::VectorXd y{{1.0, 0.05, 3.0, 1.0, 1.0, 2.0, 1.0, 1e-14, 1.0,   1.0,
                     1.0, 1.0,  1.0, 1.0, 1.0, 1.0, 1.0, 1.0,   100.0, 1e10}};
  Eigen::VectorXd a{{2.5, 2.5, 2.5, 0.2, 15.0, 1e-13, 2.5e10, 2.5, 2.5, 2.5,
                     2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5}};
  Eigen::VectorXd t0{{2.2, 0.01, 2.2, 0.2, 0.2, 1.99, 1e-13, 0.0, 0.2, 0.2, 0.2,
                      0.2, 0.2, 0.2, 0.2, 0.2, 1e13, 0.2, 0.2}};
  Eigen::VectorXd w{{1e-13, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1e-13, 0.999999999999,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}};
  Eigen::VectorXd v{{2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1e-12, 0.0, 2.0, 0.0, -5.0,
                     8.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4e10, -4e10}};
  Eigen::VectorXd sv{{0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1e-14, 0.2, 0.2,
                      0.2, 0.02, 9.0, 0.01, 0.2, 4, 0.2, 1e14}};
  Eigen::VectorXd sw{{1, 1e-13, 0.9999999999, 0.1, 0.4, 0.1, 0.1, 0.8, 0.1, 0.1,
                      0.1, 0.1, 0.1, 0.1, 0.1, 0.02, 0.4, 0.1, 0.1}};
  Eigen::VectorXd st0{{0.0, 1e-13, 1e16, 0.3, 0.3, 0.3, 1, 10, 100, 0.3, 0.3,
                       0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.02, 5.0}};
};

class Wiener5MixArgs : public ::testing::Test {
 public:
  const char* function = "function";
  void SetUp() {}

  Eigen::VectorXd y{{1.0, 0.05, 3.0, 1.0, 1.0, 2.0, 1.0, 1e-14, 1.0,   1.0,
                     1.0, 1.0,  1.0, 1.0, 1.0, 1.0, 1.0, 1.0,   100.0, 1e10}};
  Eigen::VectorXd a{{2.5, 2.5, 2.5, 0.2, 15.0, 1e-13, 2.5e10, 2.5, 2.5, 2.5,
                     2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5}};
  Eigen::VectorXd t0{{2.2, 0.01, 2.2, 0.2, 0.2, 1.99, 1e-13, 0.0, 0.2, 0.2, 0.2,
                      0.2, 0.2, 0.2, 0.2, 0.2, 1e13, 0.2, 0.2}};
  Eigen::VectorXd w{{0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1e-13, 0.999999999999,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}};
  Eigen::VectorXd v{{2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1e-12, 0.0, 2.0, 0.0, -5.0,
                     8.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4e10, -4e10}};
  Eigen::VectorXd sv{{0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1e-14, 0.2, 0.2,
                      0.2, 0.02, 9.0, 0.01, 0.2, 4, 0.2, 1e14}};
};

#endif
