#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

using Eigen::Matrix;
using Eigen::Dynamic;
using stan::math::fvar;

TEST(AgradFwdMatrixResize, fvar_double) {
  fvar<double> x = 5;
  std::vector<int> dims;
  stan::math::resize(x,dims);
}
TEST(AgradFwdMatrixResize, svec_fvar_double) {
  std::vector<fvar<double> > y;
  std::vector<int> dims;
  EXPECT_EQ(0U, y.size());

  dims.push_back(4U);
  stan::math::resize(y,dims);
  EXPECT_EQ(4U, y.size());

  dims[0] = 2U;
  stan::math::resize(y,dims);
  EXPECT_EQ(2U, y.size());
}
TEST(AgradFwdMatrixResize, vec_fvar_double) {
  Matrix<fvar<double>,Dynamic,1> v(2);
  std::vector<int> dims;
  EXPECT_EQ(2, v.size());

  dims.push_back(17U);
  dims.push_back(1U);
  stan::math::resize(v,dims);
  EXPECT_EQ(17, v.size());

  dims[0] = 3U;
  stan::math::resize(v,dims);
  EXPECT_EQ(3, v.size());
}
TEST(AgradFwdMatrixResize, rvec_fvar_double) {
  Matrix<fvar<double>,1,Dynamic> rv(2);
  std::vector<int> dims;
  EXPECT_EQ(2, rv.size());

  dims.push_back(1U);
  dims.push_back(17U);
  stan::math::resize(rv,dims);
  EXPECT_EQ(17, rv.size());

  dims[1] = 3U;
  stan::math::resize(rv,dims);
  EXPECT_EQ(3, rv.size());
}
TEST(AgradFwdMatrixResize, mat_fvar_double) {
  Matrix<fvar<double>,Dynamic,Dynamic> m(2,3);
  std::vector<int> dims;
  EXPECT_EQ(2, m.rows());
  EXPECT_EQ(3, m.cols());

  dims.push_back(7U);
  dims.push_back(17U);
  stan::math::resize(m,dims);
  EXPECT_EQ(7, m.rows());
  EXPECT_EQ(17, m.cols());
}
TEST(AgradFwdMatrixResize, svec_svec_fvar_double) {
  std::vector<std::vector<fvar<double> > > xx;
  EXPECT_EQ(0U,xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  stan::math::resize(xx,dims);
  EXPECT_EQ(4U,xx.size());
  EXPECT_EQ(5U,xx[0].size());

  dims[0] = 3U;
  dims[1] = 7U;
  stan::math::resize(xx,dims);
  EXPECT_EQ(3U,xx.size());
  EXPECT_EQ(7U,xx[1].size());  
}
TEST(AgradFwdMatrixResize, svec_v_fvar_double) {
  std::vector<Matrix<fvar<double>,Dynamic,1> > xx;
  EXPECT_EQ(0U,xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  dims.push_back(1U);
  stan::math::resize(xx,dims);
  EXPECT_EQ(4U,xx.size());
  EXPECT_EQ(5,xx[0].size());

  dims[0] = 3U;
  dims[1] = 7U;
  stan::math::resize(xx,dims);
  EXPECT_EQ(3U,xx.size());
  EXPECT_EQ(7,xx[1].size());  
}
TEST(AgradFwdMatrixResize, svec_rv_fvar_double) {
  std::vector<Matrix<fvar<double>,1,Dynamic> > xx;
  EXPECT_EQ(0U,xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(1U);
  dims.push_back(5U);
  stan::math::resize(xx,dims);
  EXPECT_EQ(4U,xx.size());
  EXPECT_EQ(5,xx[0].size());

  dims[0] = 3U;
  dims[2] = 7U;
  stan::math::resize(xx,dims);
  EXPECT_EQ(3U,xx.size());
  EXPECT_EQ(7,xx[1].size());  
}
TEST(AgradFwdMatrixResize, svec_svec_matrix_fvar_double) {
  std::vector<std::vector<Matrix<fvar<double>,Dynamic,Dynamic> > > mm;
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  dims.push_back(6U);
  dims.push_back(3U);
  stan::math::resize(mm,dims);
  EXPECT_EQ(4U,mm.size());
  EXPECT_EQ(5U,mm[0].size());
  EXPECT_EQ(6,mm[1][2].rows());
  EXPECT_EQ(3,mm[3][4].cols());
}
TEST(AgradFwdMatrixResize, fvar_fvar_double) {
  fvar<fvar<double> > x = 5;
  std::vector<int> dims;
  stan::math::resize(x,dims);
}
TEST(AgradFwdMatrixResize, svec_fvar_fvar_double) {
  std::vector<fvar<fvar<double> > > y;
  std::vector<int> dims;
  EXPECT_EQ(0U, y.size());

  dims.push_back(4U);
  stan::math::resize(y,dims);
  EXPECT_EQ(4U, y.size());

  dims[0] = 2U;
  stan::math::resize(y,dims);
  EXPECT_EQ(2U, y.size());
}
TEST(AgradFwdMatrixResize, vec_fvar_fvar_double) {
  Matrix<fvar<fvar<double> >,Dynamic,1> v(2);
  std::vector<int> dims;
  EXPECT_EQ(2, v.size());

  dims.push_back(17U);
  dims.push_back(1U);
  stan::math::resize(v,dims);
  EXPECT_EQ(17, v.size());

  dims[0] = 3U;
  stan::math::resize(v,dims);
  EXPECT_EQ(3, v.size());
}
TEST(AgradFwdMatrixResize, rvec_fvar_fvar_double) {
  Matrix<fvar<fvar<double> >,1,Dynamic> rv(2);
  std::vector<int> dims;
  EXPECT_EQ(2, rv.size());

  dims.push_back(1U);
  dims.push_back(17U);
  stan::math::resize(rv,dims);
  EXPECT_EQ(17, rv.size());

  dims[1] = 3U;
  stan::math::resize(rv,dims);
  EXPECT_EQ(3, rv.size());
}
TEST(AgradFwdMatrixResize, mat_fvar_fvar_double) {
  Matrix<fvar<fvar<double> >,Dynamic,Dynamic> m(2,3);
  std::vector<int> dims;
  EXPECT_EQ(2, m.rows());
  EXPECT_EQ(3, m.cols());

  dims.push_back(7U);
  dims.push_back(17U);
  stan::math::resize(m,dims);
  EXPECT_EQ(7, m.rows());
  EXPECT_EQ(17, m.cols());
}
TEST(AgradFwdMatrixResize, svec_svec_fvar_fvar_double) {
  std::vector<std::vector<fvar<fvar<double> > > > xx;
  EXPECT_EQ(0U,xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  stan::math::resize(xx,dims);
  EXPECT_EQ(4U,xx.size());
  EXPECT_EQ(5U,xx[0].size());

  dims[0] = 3U;
  dims[1] = 7U;
  stan::math::resize(xx,dims);
  EXPECT_EQ(3U,xx.size());
  EXPECT_EQ(7U,xx[1].size());  
}
TEST(AgradFwdMatrixResize, svec_v_fvar_fvar_double) {
  std::vector<Matrix<fvar<fvar<double> >,Dynamic,1> > xx;
  EXPECT_EQ(0U,xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  dims.push_back(1U);
  stan::math::resize(xx,dims);
  EXPECT_EQ(4U,xx.size());
  EXPECT_EQ(5,xx[0].size());

  dims[0] = 3U;
  dims[1] = 7U;
  stan::math::resize(xx,dims);
  EXPECT_EQ(3U,xx.size());
  EXPECT_EQ(7,xx[1].size());  
}
TEST(AgradFwdMatrixResize, svec_rv_fvar_fvar_double) {
  std::vector<Matrix<fvar<fvar<double> >,1,Dynamic> > xx;
  EXPECT_EQ(0U,xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(1U);
  dims.push_back(5U);
  stan::math::resize(xx,dims);
  EXPECT_EQ(4U,xx.size());
  EXPECT_EQ(5,xx[0].size());

  dims[0] = 3U;
  dims[2] = 7U;
  stan::math::resize(xx,dims);
  EXPECT_EQ(3U,xx.size());
  EXPECT_EQ(7,xx[1].size());  
}
TEST(AgradFwdMatrixResize, svec_svec_matrix_fvar_fvar_double) {
  std::vector<std::vector<Matrix<fvar<fvar<double> >,Dynamic,Dynamic> > > mm;
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  dims.push_back(6U);
  dims.push_back(3U);
  stan::math::resize(mm,dims);
  EXPECT_EQ(4U,mm.size());
  EXPECT_EQ(5U,mm[0].size());
  EXPECT_EQ(6,mm[1][2].rows());
  EXPECT_EQ(3,mm[3][4].cols());
}
