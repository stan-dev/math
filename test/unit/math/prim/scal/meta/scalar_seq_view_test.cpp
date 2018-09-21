#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, ScalarSeqViewDouble) {
  using stan::scalar_seq_view;

  double d = 10;
  scalar_seq_view<double> sv(d);
  EXPECT_FLOAT_EQ(d, sv[0]);
  EXPECT_FLOAT_EQ(d, sv[12]);

  EXPECT_EQ(1, sv.size());
}
