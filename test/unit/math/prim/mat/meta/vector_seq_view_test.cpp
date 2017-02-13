#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>

TEST(MetaTraits, VectorSeqView) {
  using Eigen::MatrixXd;
  using std::vector;
  using stan::vector_seq_view;

  MatrixXd m1(2, 2);
  m1 << 1.1, 2.2, 3.3, 4.4;

  vector_seq_view<MatrixXd> vsv(m1);
  EXPECT_FLOAT_EQ(m1(1), vsv[12](1));

  vector<MatrixXd> v;
  v.push_back(m1);
  v.push_back(m1.transpose());
  vector_seq_view<vector<MatrixXd> > vsv_vec(v);
  EXPECT_FLOAT_EQ(m1(1), vsv_vec[0](1));
  EXPECT_FLOAT_EQ(m1(1), vsv_vec[1](2));
  EXPECT_FLOAT_EQ(m1(2), vsv_vec[1](1));
  EXPECT_FLOAT_EQ(m1(2), vsv_vec[0](2));
}
