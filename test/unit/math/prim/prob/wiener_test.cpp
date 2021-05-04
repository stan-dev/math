#include <stan/math/prim/prob.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(mathPrimScalProbWienerMat, illegal_tau_gt_y) {
  using stan::math::wiener_log;
  std::vector<double> ys;
  ys.push_back(2.5);
  ys.push_back(0.4);
  EXPECT_THROW(wiener_log<true>(ys, 4.1, 1.9, 0.05, 0.1), std::domain_error);
}

TEST(mathPrimScalProbWienerScal, illegal_tau_gt_y) {
  using stan::math::wiener_log;
  EXPECT_THROW(wiener_log<true>(0.4, 4.1, 1.9, 0.05, 0.1), std::domain_error);
}
