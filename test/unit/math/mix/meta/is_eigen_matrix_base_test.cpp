#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_eigen_matrix_base_hierarchy_tests) {
  using stan::is_eigen_matrix_base;
  using stan::math::test::all_eigen_dense;
  using Eigen::Matrix;
  using Eigen::Array;
  all_eigen_dense<false, true, true, true, true, is_eigen_matrix_base, Matrix, -1, -1>();
  all_eigen_dense<false, true, true, true, true, is_eigen_matrix_base, Matrix, 1, -1>();
  all_eigen_dense<false, true, true, true, true, is_eigen_matrix_base, Matrix, -1, 1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_matrix_base, Array, -1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_matrix_base, Array, 1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_matrix_base, Array, -1, 1>();
}

TEST(MathMetaPrim, is_eigen_matrix_base_sparse_tests) {
  using stan::is_eigen_matrix_base;
  using stan::math::test::all_eigen_sparse;
  using Eigen::Matrix;
  using Eigen::Array;
  all_eigen_sparse<false, false, false, false, is_eigen_matrix_base>();
}

TEST(MathMetaPrim, is_eigen_matrix_base_decomp_tests) {
  using stan::is_eigen_matrix_base;
  using stan::math::test::all_eigen_dense_decomp;
  using Eigen::Matrix;
  using Eigen::Array;
  all_eigen_dense_decomp<false, is_eigen_matrix_base>();

}

TEST(MathMetaPrim, is_eigen_matrix_base_expr_tests) {
  using stan::is_eigen_matrix_base;
  using stan::math::test::all_eigen_dense_exprs;
  using Eigen::Matrix;
  using Eigen::Array;
 all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_base, Matrix, -1, -1>;
 all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_base, Matrix, 1, -1>;
 all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_base, Matrix, -1, 1>;
 all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_base, Array, -1, -1>;
 all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_base, Array, 1, -1>;
 all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_base, Array, -1, 1>;
}
