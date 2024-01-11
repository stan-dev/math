#ifndef STAN_TEST_UNIT_MATH_MIX_PROB_UTIL_HPP
#define STAN_TEST_UNIT_MATH_MIX_PROB_UTIL_HPP

#include <gtest/gtest.h>

class Wiener7MixArgs : public ::testing::Test {
 public:
  const char* function = "function";
  void SetUp() {}

  Eigen::VectorXd y{{0.1, 0.3}};
  Eigen::VectorXd a{{2.0, 2.5}};
  Eigen::VectorXd t0{{0.2, 0.1}};
  Eigen::VectorXd w{{0.5, 0.2}};
  Eigen::VectorXd v{{2.0, 0.8}};
  Eigen::VectorXd sv{{0.2, 0.1}};
  Eigen::VectorXd sw{{0.1, .99}};
  Eigen::VectorXd st0{{0.3, .8}};
};

class Wiener5MixArgs : public ::testing::Test {
 public:
  const char* function = "function";
  void SetUp() {}

  Eigen::VectorXd y{{1.0, 0.1, 0.3}};
  Eigen::VectorXd a{{2.5, 2.0, 2.5}};
  Eigen::VectorXd t0{{0.2, 0.2, 0.1}};
  Eigen::VectorXd w{{0.5, 0.5, 0.2}};
  Eigen::VectorXd v{{2.0, 2.0, 0.8}};
  Eigen::VectorXd sv{{0.2, 0.2, 0.1}};
};

#endif
