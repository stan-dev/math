#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, ScalarSeqViewDouble) {
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
