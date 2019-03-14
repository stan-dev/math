
#include <stan/math/prim.hpp>
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

TEST(MetaTraits_arr, ScalarSeqViewArray) {
  using stan::scalar_seq_view;
  using std::vector;

  vector<double> v;
  v.push_back(2.2);
  v.push_back(0.0001);
  scalar_seq_view<vector<double> > sv(v);
  EXPECT_FLOAT_EQ(v[0], sv[0]);
  EXPECT_FLOAT_EQ(v[1], sv[1]);

  EXPECT_EQ(v.size(), sv.size());
}

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

TEST(MetaTraits_mat, ScalarSeqViewVector) {
  expect_scalar_seq_view_values(Eigen::VectorXd(4));
}

TEST(MetaTraits_mat, ScalarSeqViewRowVector) {
  expect_scalar_seq_view_values(Eigen::RowVectorXd(4));
}
