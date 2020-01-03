#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbAgradDistributionsStdNormal, fwd) {
  using stan::math::fvar;
  using stan::math::std_normal_lpdf;

  EXPECT_FLOAT_EQ(-0.918938533204673, std_normal_lpdf<false>(0));
  EXPECT_FLOAT_EQ(-0.918938533204673,
                  std_normal_lpdf<false>(fvar<double>(0)).val());
}
