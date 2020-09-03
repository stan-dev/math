#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_eigen_sparse_base_hierarchy_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_sparse_base;
  using stan::math::test::all_eigen_dense;
  all_eigen_dense<false, false, false, false, false, is_eigen_sparse_base,
                  Matrix, -1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_sparse_base,
                  Matrix, 1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_sparse_base,
                  Matrix, -1, 1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_sparse_base,
                  Array, -1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_sparse_base,
                  Array, 1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_sparse_base,
                  Array, -1, 1>();
}

TEST(MathMetaPrim, is_eigen_sparse_base_sparse_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_sparse_base;
  using stan::math::test::all_eigen_sparse;
  all_eigen_sparse<false, true, true, true, is_eigen_sparse_base>();
}

TEST(MathMetaPrim, is_eigen_sparse_base_decomp_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_sparse_base;
  using stan::math::test::all_eigen_dense_decomp;
  all_eigen_dense_decomp<false, is_eigen_sparse_base>();
}

TEST(MathMetaPrim, is_eigen_sparse_base_expr_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_sparse_base;
  using stan::math::test::all_eigen_dense_exprs;
  all_eigen_dense_exprs<false, false, false, false, is_eigen_sparse_base,
                        Matrix, -1, -1>();
  all_eigen_dense_exprs<false, false, false, false, is_eigen_sparse_base,
                        Matrix, 1, -1>();
  all_eigen_dense_exprs<false, false, false, false, is_eigen_sparse_base,
                        Matrix, -1, 1>();
  all_eigen_dense_exprs<false, false, false, false, is_eigen_sparse_base, Array,
                        -1, -1>();
  all_eigen_dense_exprs<false, false, false, false, is_eigen_sparse_base, Array,
                        1, -1>();
  all_eigen_dense_exprs<false, false, false, false, is_eigen_sparse_base, Array,
                        -1, 1>();
}
