#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, ScalarSeqViewDouble) {
  using stan::scalar_seq_view;

  double d = 10;
  scalar_seq_view<double> sv(d);
  EXPECT_FLOAT_EQ(d, sv[0]);
  EXPECT_FLOAT_EQ(d, sv[12]);

  const double d_const = 10;
  scalar_seq_view<const double> sv_const(d_const);
  EXPECT_FLOAT_EQ(d_const, sv_const[0]);
  EXPECT_FLOAT_EQ(d_const, sv_const[12]);

  const double& d_const_ref = 10;
  scalar_seq_view<const double&> sv_const_ref(d_const);
  EXPECT_FLOAT_EQ(d_const_ref, sv_const_ref[0]);
  EXPECT_FLOAT_EQ(d_const_ref, sv_const_ref[12]);

  EXPECT_EQ(1, sv.size());

  double* d_point;
  d_point = static_cast<double*>(malloc(sizeof(double) * 2));
  d_point[0] = 69.0;
  d_point[1] = 420.0;

  scalar_seq_view<decltype(d_point)> d_point_v(d_point);
  EXPECT_FLOAT_EQ(69.0, d_point_v[0]);
  EXPECT_FLOAT_EQ(420.0, d_point_v[1]);
  free(d_point);
}

TEST(MathMetaPrim, ScalarSeqViewArray) {
  using stan::scalar_seq_view;
  using std::vector;

  vector<double> v;
  v.push_back(2.2);
  v.push_back(0.0001);
  scalar_seq_view<vector<double> > sv(v);
  EXPECT_FLOAT_EQ(v[0], sv[0]);
  EXPECT_FLOAT_EQ(v[1], sv[1]);

  const vector<double> v_const{2.2, 0.001};
  scalar_seq_view<const vector<double> > sv_const(v_const);
  EXPECT_FLOAT_EQ(v_const[0], sv_const[0]);
  EXPECT_FLOAT_EQ(v_const[1], sv_const[1]);

  const vector<double>& v_const_ref{2.2, 0.001};
  scalar_seq_view<const vector<double> > sv_const_ref(v_const_ref);
  EXPECT_FLOAT_EQ(v_const_ref[0], sv_const_ref[0]);
  EXPECT_FLOAT_EQ(v_const_ref[1], sv_const_ref[1]);

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

TEST(MathMetaPrim, ScalarSeqViewVector) {
  expect_scalar_seq_view_values(Eigen::VectorXd(4));
}

TEST(MathMetaPrim, ScalarSeqViewRowVector) {
  expect_scalar_seq_view_values(Eigen::RowVectorXd(4));
}
