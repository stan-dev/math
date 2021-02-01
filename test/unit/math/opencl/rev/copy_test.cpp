#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
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

TEST(MathMatrixGPU, var_std_vector_to_matrix_cl) {
  using stan::math::var_value;
  std::vector<stan::math::var> vals = {1, 2, 3, 4, 5, 6};
  var_value<stan::math::matrix_cl<double>> a_cl
      = stan::math::to_matrix_cl(vals);
  stan::math::matrix_d a_res = from_matrix_cl(a_cl.val());
  for (size_t i = 0; i < vals.size(); i++) {
    EXPECT_FLOAT_EQ(a_res(i), vals[i].val());
  }
  a_cl.adj() = stan::math::constant(1, vals.size(), 1);
  a_cl.vi_->chain();
  for (size_t i = 0; i < vals.size(); i++) {
    EXPECT_FLOAT_EQ(vals[i].adj(), 1);
  }
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

TEST(MathMatrixGPU, std_vector_matrix_var_to_matrix_cl) {
  using stan::math::var_value;
  stan::math::matrix_v a(2, 3);
  a << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_v b(2, 3);
  b << 7, 8, 9, 0, 1, 2;
  stan::math::matrix_d a_adj(2, 3);
  a_adj << 9, 8, 7, 6, 5, 4;
  stan::math::matrix_d b_adj(2, 3);
  b_adj << 6, 5, 4, 3, 2, 1;
  std::vector<stan::math::matrix_v> vars = {a, b};
  std::vector<stan::math::matrix_d> adjs = {a_adj, b_adj};
  var_value<stan::math::matrix_cl<double>> vars_cl
      = stan::math::to_matrix_cl(vars);

  vars_cl.adj() = stan::math::to_matrix_cl(adjs);
  vars_cl.vi_->chain();
  EXPECT_MATRIX_EQ(vars[0].adj(), a_adj);
  EXPECT_MATRIX_EQ(vars[1].adj(), b_adj);

  stan::math::matrix_d vars_res = stan::math::from_matrix_cl(vars_cl.val());
  for (size_t i = 0; i < a.size(); i++) {
    EXPECT_FLOAT_EQ(vars_res(i), a(i).val());
    EXPECT_FLOAT_EQ(vars_res(i + a.size()), b(i).val());
  }
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
