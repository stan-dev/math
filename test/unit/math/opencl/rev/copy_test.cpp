#ifdef STAN_OPENCL
#include <stan/math/opencl/rev/opencl.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <CL/cl2.hpp>
#include <algorithm>
#include <vector>

TEST(VariCL, var_matrix_to_matrix_cl) {
  using stan::math::var_value;
  Eigen::MatrixXd vals(2, 3);
  vals << 1, 2, 3, 4, 5, 6;
  var_value<Eigen::MatrixXd> a(vals);
  var_value<stan::math::matrix_cl<double>> a_cl = stan::math::to_matrix_cl(a);
  EXPECT_MATRIX_EQ(from_matrix_cl(a_cl.val()), vals);
  a_cl.adj() = stan::math::constant(1, 2, 3);
  a_cl.vi_->chain();
  EXPECT_MATRIX_EQ(a.adj(), Eigen::MatrixXd::Constant(2, 3, 1));
}

TEST(VariCL, matrix_var_to_matrix_cl) {
  using stan::math::var_value;
  Eigen::MatrixXd vals(2, 3);
  vals << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_v vars = vals;
  var_value<stan::math::matrix_cl<double>> a_cl
      = stan::math::to_matrix_cl(vars);
  EXPECT_MATRIX_EQ(from_matrix_cl(a_cl.val()), vals);
  a_cl.adj() = stan::math::constant(1, 2, 3);
  a_cl.vi_->chain();
  EXPECT_MATRIX_EQ(vars.adj(), Eigen::MatrixXd::Constant(2, 3, 1));
}

TEST(VariCL, from_matrix_cl) {
  using stan::math::var_value;
  Eigen::MatrixXd vals(2, 3);
  vals << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_cl<double> vals_cl(vals);
  var_value<stan::math::matrix_cl<double>> a_cl(vals_cl);
  var_value<Eigen::MatrixXd> a = stan::math::from_matrix_cl(a_cl);
  EXPECT_MATRIX_EQ(a.val(), vals);
  a.adj() = Eigen::MatrixXd::Constant(2, 3, 1);
  a.vi_->chain();
  EXPECT_MATRIX_EQ(from_matrix_cl(a_cl.adj()),
                   Eigen::MatrixXd::Constant(2, 3, 1));
}

#endif
