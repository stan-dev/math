#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <vector>


TEST(ProbAgradDistributionsNormal, fwd) {
  using stan::math::fvar;
  using stan::math::normal_log;

  EXPECT_FLOAT_EQ(-0.918938533204673,
                  normal_log<false>(0, 0, 1));
  EXPECT_FLOAT_EQ(-0.918938533204673,
                  normal_log<false>(0, 0, fvar<double>(1.0)).val());
  EXPECT_FLOAT_EQ(-0.918938533204673,
                  normal_log<false>(0, fvar<double>(0), 1).val());
  EXPECT_FLOAT_EQ(-0.918938533204673,
                  normal_log<false>(0, fvar<double>(0), fvar<double>(1)).val());
  EXPECT_FLOAT_EQ(-0.918938533204673,
                  normal_log<false>(fvar<double>(0), 0, 1).val());
  EXPECT_FLOAT_EQ(-0.918938533204673,
                  normal_log<false>(fvar<double>(0), 0, fvar<double>(1)).val());
  EXPECT_FLOAT_EQ(-0.918938533204673,
                  normal_log<false>(fvar<double>(0), fvar<double>(0), 1).val());
  EXPECT_FLOAT_EQ(-0.918938533204673,
                  normal_log<false>(fvar<double>(0), fvar<double>(0),
                                    fvar<double>(1)).val());
}
