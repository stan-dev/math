#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, VectorSeqView) {
  using Eigen::VectorXd;
  using stan::vector_seq_view;
  using std::vector;

  VectorXd m1(4);
  m1 << 1.1, 2.2, 3.3, 4.4;

  vector_seq_view<VectorXd> vsv(m1);
  EXPECT_FLOAT_EQ(m1(1), vsv[12](1));
  EXPECT_EQ(vsv.size(), 1);

  vector<VectorXd> v;
  v.push_back(m1);
  v.push_back(m1.reverse());
  vector_seq_view<vector<VectorXd> > vsv_vec(v);
  EXPECT_FLOAT_EQ(m1(0), vsv_vec[0][0]);
  EXPECT_FLOAT_EQ(m1(0), vsv_vec[1][3]);
  EXPECT_FLOAT_EQ(m1(2), vsv_vec[1][1]);
  EXPECT_FLOAT_EQ(m1(2), vsv_vec[0][2]);
  EXPECT_EQ(vsv_vec.size(), 2);
}
