#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>

template <typename C>
void expect_scalar_seq_view_values(C v) {
  using stan::scalar_seq_view;

  v << 1.1, 2.2, 3.3, 4.4;
  scalar_seq_view<C> sv(v);
  EXPECT_FLOAT_EQ(v(0), sv[0]);
  EXPECT_FLOAT_EQ(v(1), sv[1]);
  EXPECT_FLOAT_EQ(v(2), sv[2]);
  EXPECT_FLOAT_EQ(v(3), sv[3]);

  EXPECT_EQ(v.size(), sv.size());
}

TEST(MetaTraits, ScalarSeqViewVector) {
  expect_scalar_seq_view_values(Eigen::VectorXd(4));
}

TEST(MetaTraits, ScalarSeqViewRowVector) {
  expect_scalar_seq_view_values(Eigen::RowVectorXd(4));
}
