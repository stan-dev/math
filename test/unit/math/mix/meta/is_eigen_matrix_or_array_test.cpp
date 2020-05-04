#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_eigen_matrix_or_array_hierarchy_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_matrix_or_array;
  using stan::math::test::all_eigen_dense;
  all_eigen_dense<false, false, true, true, true, is_eigen_matrix_or_array, Matrix,
                  -1, -1>();
  all_eigen_dense<false, false, true, true, true, is_eigen_matrix_or_array, Matrix,
                  1, -1>();
  all_eigen_dense<false, false, true, true, true, is_eigen_matrix_or_array, Matrix,
                  -1, 1>();
  all_eigen_dense<false, false, true, true, true, is_eigen_matrix_or_array, Array,
                  -1, -1>();
  all_eigen_dense<false, false, true, true, true, is_eigen_matrix_or_array, Array,
                  1, -1>();
  all_eigen_dense<false, false, true, true, true, is_eigen_matrix_or_array, Array,
                  -1, 1>();
}

TEST(MathMetaPrim, is_eigen_matrix_or_array_sparse_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_matrix_or_array;
  using stan::math::test::all_eigen_sparse;
  all_eigen_sparse<false, false, false, false, is_eigen_matrix_or_array>();
}

TEST(MathMetaPrim, is_eigen_matrix_or_array_decomp_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_matrix_or_array;
  using stan::math::test::all_eigen_dense_decomp;
  all_eigen_dense_decomp<false, is_eigen_matrix_or_array>();
}

TEST(MathMetaPrim, is_eigen_matrix_or_array_expr_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_matrix_or_array;
  using stan::math::test::all_eigen_dense_exprs;
  all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_or_array, Matrix,
                        -1, -1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_or_array, Matrix, 1,
                        -1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_or_array, Matrix,
                        -1, 1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_or_array, Array,
                        -1, -1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_or_array, Array,
                        1, -1>();
  all_eigen_dense_exprs<true, true, true, true, is_eigen_matrix_or_array, Array,
                        -1, 1>();
}
