#include <stan/math/prim/prob.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(mathPrimScalProbWienerMat, illegal_tau_gt_y) {
  using stan::math::wiener_lpdf;
  std::vector<double> ys;
  ys.push_back(2.5);
  ys.push_back(0.4);
  EXPECT_THROW(wiener_lpdf<true>(ys, 4.1, 1.9, 0.05, 0.1), std::domain_error);
}

TEST(mathPrimScalProbWienerScal, illegal_tau_gt_y) {
  using stan::math::wiener_lpdf;
  EXPECT_THROW(wiener_lpdf<true>(0.4, 4.1, 1.9, 0.05, 0.1), std::domain_error);
}
