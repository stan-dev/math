#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_eigen_hierarchy_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen;
  using stan::math::test::all_eigen_dense;
  all_eigen_dense<true, true, true, true, true, is_eigen, Matrix, -1, -1>();
  all_eigen_dense<true, true, true, true, true, is_eigen, Matrix, 1, -1>();
  all_eigen_dense<true, true, true, true, true, is_eigen, Matrix, -1, 1>();
  all_eigen_dense<true, true, true, true, true, is_eigen, Array, -1, -1>();
  all_eigen_dense<true, true, true, true, true, is_eigen, Array, 1, -1>();
  all_eigen_dense<true, true, true, true, true, is_eigen, Array, -1, 1>();
}

TEST(MathMetaPrim, is_eigen_sparse_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen;
  using stan::math::test::all_eigen_sparse;
  all_eigen_sparse<true, true, true, true, is_eigen>();
}

TEST(MathMetaPrim, is_eigen_decomp_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen;
  using stan::math::test::all_eigen_dense_decomp;
  all_eigen_dense_decomp<true, is_eigen>();
}

TEST(MathMetaPrim, is_eigen_expr_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen;
  using stan::math::test::all_eigen_dense_exprs;
  all_eigen_dense_exprs<true, true, true, true, is_eigen, Matrix, -1, -1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen, Matrix, 1, -1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen, Matrix, -1, 1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen, Array, -1, -1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen, Array, 1, -1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen, Array, -1, 1>();
}
