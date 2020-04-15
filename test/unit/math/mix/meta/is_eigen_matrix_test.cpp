#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_eigen_matrix_hierarchy_tests) {
  using stan::is_eigen_matrix;
  using stan::math::test::all_eigen_dense;
  using Eigen::Matrix;
  using Eigen::Array;
  all_eigen_dense<false, true, true, true, true, is_eigen_matrix, Matrix, -1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_matrix, Matrix, 1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_matrix, Matrix, -1, 1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_matrix, Array, -1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_matrix, Array, 1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_matrix, Array, -1, 1>();
}

TEST(MathMetaPrim, is_eigen_matrix_sparse_tests) {
  using stan::is_eigen_matrix;
  using stan::math::test::all_eigen_sparse;
  using Eigen::Matrix;
  using Eigen::Array;
  all_eigen_sparse<false, false, false, false, is_eigen_matrix>();
}

TEST(MathMetaPrim, is_eigen_matrix_decomp_tests) {
  using stan::is_eigen_matrix;
  using stan::math::test::all_eigen_dense_decomp;
  using Eigen::Matrix;
  using Eigen::Array;
  all_eigen_dense_decomp<false, is_eigen_matrix>();

}

TEST(MathMetaPrim, is_eigen_matrix_expr_tests) {
  using stan::is_eigen_matrix;
  using stan::math::test::all_eigen_dense_exprs;
  using Eigen::Matrix;
  using Eigen::Array;
 all_eigen_dense_exprs<true, true, false, true, is_eigen_matrix, Matrix, -1, -1>;
 all_eigen_dense_exprs<true, true, false, true, is_eigen_matrix, Matrix, 1, -1>;
 all_eigen_dense_exprs<true, true, false, true, is_eigen_matrix, Matrix, -1, 1>;
 all_eigen_dense_exprs<true, true, false, true, is_eigen_matrix, Array, -1, -1>;
 all_eigen_dense_exprs<true, true, false, true, is_eigen_matrix, Array, 1, -1>;
 all_eigen_dense_exprs<true, true, false, true, is_eigen_matrix, Array, -1, 1>;
}
