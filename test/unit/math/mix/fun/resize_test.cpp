#include <stan/math/mix.hpp>
#include <test/unit/math/mix/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST_F(mathMix, resize_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> x = 5;
  std::vector<int> dims;
  stan::math::resize(x, dims);
}
TEST_F(mathMix, resize_svec_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<fvar<var> > y;
  std::vector<int> dims;
  EXPECT_EQ(0U, y.size());

  dims.push_back(4U);
  stan::math::resize(y, dims);
  EXPECT_EQ(4U, y.size());

  dims[0] = 2U;
  stan::math::resize(y, dims);
  EXPECT_EQ(2U, y.size());
}
TEST_F(mathMix, resize_vec_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<var>, Dynamic, 1> v(2);
  std::vector<int> dims;
  EXPECT_EQ(2, v.size());

  dims.push_back(17U);
  dims.push_back(1U);
  stan::math::resize(v, dims);
  EXPECT_EQ(17, v.size());

  dims[0] = 3U;
  stan::math::resize(v, dims);
  EXPECT_EQ(3, v.size());
}
TEST_F(mathMix, resize_rvec_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<var>, 1, Dynamic> rv(2);
  std::vector<int> dims;
  EXPECT_EQ(2, rv.size());

  dims.push_back(1U);
  dims.push_back(17U);
  stan::math::resize(rv, dims);
  EXPECT_EQ(17, rv.size());

  dims[1] = 3U;
  stan::math::resize(rv, dims);
  EXPECT_EQ(3, rv.size());
}
TEST_F(mathMix, resize_mat_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<var>, Dynamic, Dynamic> m(2, 3);
  std::vector<int> dims;
  EXPECT_EQ(2, m.rows());
  EXPECT_EQ(3, m.cols());

  dims.push_back(7U);
  dims.push_back(17U);
  stan::math::resize(m, dims);
  EXPECT_EQ(7, m.rows());
  EXPECT_EQ(17, m.cols());
}
TEST_F(mathMix, resize_svec_svec_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<std::vector<fvar<var> > > xx;
  EXPECT_EQ(0U, xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  stan::math::resize(xx, dims);
  EXPECT_EQ(4U, xx.size());
  EXPECT_EQ(5U, xx[0].size());

  dims[0] = 3U;
  dims[1] = 7U;
  stan::math::resize(xx, dims);
  EXPECT_EQ(3U, xx.size());
  EXPECT_EQ(7U, xx[1].size());
}
TEST_F(mathMix, resize_svec_v_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<Matrix<fvar<var>, Dynamic, 1> > xx;
  EXPECT_EQ(0U, xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  dims.push_back(1U);
  stan::math::resize(xx, dims);
  EXPECT_EQ(4U, xx.size());
  EXPECT_EQ(5, xx[0].size());

  dims[0] = 3U;
  dims[1] = 7U;
  stan::math::resize(xx, dims);
  EXPECT_EQ(3U, xx.size());
  EXPECT_EQ(7, xx[1].size());
}
TEST_F(mathMix, resize_svec_rv_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<Matrix<fvar<var>, 1, Dynamic> > xx;
  EXPECT_EQ(0U, xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(1U);
  dims.push_back(5U);
  stan::math::resize(xx, dims);
  EXPECT_EQ(4U, xx.size());
  EXPECT_EQ(5, xx[0].size());

  dims[0] = 3U;
  dims[2] = 7U;
  stan::math::resize(xx, dims);
  EXPECT_EQ(3U, xx.size());
  EXPECT_EQ(7, xx[1].size());
}
TEST_F(mathMix, resize_svec_svec_matrix_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<std::vector<Matrix<fvar<var>, Dynamic, Dynamic> > > mm;
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  dims.push_back(6U);
  dims.push_back(3U);
  stan::math::resize(mm, dims);
  EXPECT_EQ(4U, mm.size());
  EXPECT_EQ(5U, mm[0].size());
  EXPECT_EQ(6, mm[1][2].rows());
  EXPECT_EQ(3, mm[3][4].cols());
}
TEST_F(mathMix, resize_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  fvar<fvar<var> > x = 5;
  std::vector<int> dims;
  stan::math::resize(x, dims);
}
TEST_F(mathMix, resize_svec_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<fvar<fvar<var> > > y;
  std::vector<int> dims;
  EXPECT_EQ(0U, y.size());

  dims.push_back(4U);
  stan::math::resize(y, dims);
  EXPECT_EQ(4U, y.size());

  dims[0] = 2U;
  stan::math::resize(y, dims);
  EXPECT_EQ(2U, y.size());
}
TEST_F(mathMix, resize_vec_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<fvar<var> >, Dynamic, 1> v(2);
  std::vector<int> dims;
  EXPECT_EQ(2, v.size());

  dims.push_back(17U);
  dims.push_back(1U);
  stan::math::resize(v, dims);
  EXPECT_EQ(17, v.size());

  dims[0] = 3U;
  stan::math::resize(v, dims);
  EXPECT_EQ(3, v.size());
}
TEST_F(mathMix, resize_rvec_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<fvar<var> >, 1, Dynamic> rv(2);
  std::vector<int> dims;
  EXPECT_EQ(2, rv.size());

  dims.push_back(1U);
  dims.push_back(17U);
  stan::math::resize(rv, dims);
  EXPECT_EQ(17, rv.size());

  dims[1] = 3U;
  stan::math::resize(rv, dims);
  EXPECT_EQ(3, rv.size());
}
TEST_F(mathMix, resize_mat_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> m(2, 3);
  std::vector<int> dims;
  EXPECT_EQ(2, m.rows());
  EXPECT_EQ(3, m.cols());

  dims.push_back(7U);
  dims.push_back(17U);
  stan::math::resize(m, dims);
  EXPECT_EQ(7, m.rows());
  EXPECT_EQ(17, m.cols());
}
TEST_F(mathMix, resize_svec_svec_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<std::vector<fvar<fvar<var> > > > xx;
  EXPECT_EQ(0U, xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  stan::math::resize(xx, dims);
  EXPECT_EQ(4U, xx.size());
  EXPECT_EQ(5U, xx[0].size());

  dims[0] = 3U;
  dims[1] = 7U;
  stan::math::resize(xx, dims);
  EXPECT_EQ(3U, xx.size());
  EXPECT_EQ(7U, xx[1].size());
}
TEST_F(mathMix, resize_svec_v_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<Matrix<fvar<fvar<var> >, Dynamic, 1> > xx;
  EXPECT_EQ(0U, xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  dims.push_back(1U);
  stan::math::resize(xx, dims);
  EXPECT_EQ(4U, xx.size());
  EXPECT_EQ(5, xx[0].size());

  dims[0] = 3U;
  dims[1] = 7U;
  stan::math::resize(xx, dims);
  EXPECT_EQ(3U, xx.size());
  EXPECT_EQ(7, xx[1].size());
}
TEST_F(mathMix, resize_svec_rv_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<Matrix<fvar<fvar<var> >, 1, Dynamic> > xx;
  EXPECT_EQ(0U, xx.size());
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(1U);
  dims.push_back(5U);
  stan::math::resize(xx, dims);
  EXPECT_EQ(4U, xx.size());
  EXPECT_EQ(5, xx[0].size());

  dims[0] = 3U;
  dims[2] = 7U;
  stan::math::resize(xx, dims);
  EXPECT_EQ(3U, xx.size());
  EXPECT_EQ(7, xx[1].size());
}
TEST_F(mathMix, resize_svec_svec_matrix_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<std::vector<Matrix<fvar<fvar<var> >, Dynamic, Dynamic> > > mm;
  std::vector<int> dims;
  dims.push_back(4U);
  dims.push_back(5U);
  dims.push_back(6U);
  dims.push_back(3U);
  stan::math::resize(mm, dims);
  EXPECT_EQ(4U, mm.size());
  EXPECT_EQ(5U, mm[0].size());
  EXPECT_EQ(6, mm[1][2].rows());
  EXPECT_EQ(3, mm[3][4].cols());
}
