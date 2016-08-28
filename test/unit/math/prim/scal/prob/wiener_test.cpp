#include <stan/math/prim/scal/prob/wiener_log.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(mathPrimScalProbWiener, illegal_tau_gt_y) {
  using stan::math::wiener_log;
  EXPECT_THROW(wiener_log<true>(0.4, 4.1, 1.9, 0.05, 0.1),
               std::domain_error);

  std::vector<double> ys;
  ys.push_back(2.5, 0.4);
  EXPECT_THROW(wiener_log<true>(ys, 4.1, 1.9, 0.05, 0.1),
               std::domain_error);
  
}
