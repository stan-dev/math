#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_eigen_dense_plain_type_hierarchy_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_dense_plain_type;
  using stan::math::test::all_eigen_dense;
  all_eigen_dense<false, false, false, false, false, is_eigen_dense_plain_type,
                  Matrix, -1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_dense_plain_type,
                  Matrix, 1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_dense_plain_type,
                  Matrix, -1, 1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_dense_plain_type,
                  Array, -1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_dense_plain_type,
                  Array, 1, -1>();
  all_eigen_dense<false, false, false, false, false, is_eigen_dense_plain_type,
                  Array, -1, 1>();
}

TEST(MathMetaPrim, is_eigen_dense_plain_type_sparse_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_dense_plain_type;
  using stan::math::test::all_eigen_sparse;
  all_eigen_sparse<false, false, false, false, is_eigen_dense_plain_type>();
}

TEST(MathMetaPrim, is_eigen_dense_plain_type_decomp_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_dense_plain_type;
  using stan::math::test::all_eigen_dense_decomp;
  all_eigen_dense_decomp<false, is_eigen_dense_plain_type>();
}

TEST(MathMetaPrim, is_eigen_dense_plain_type_expr_tests) {
  using Eigen::Array;
  using Eigen::Matrix;
  using stan::is_eigen_dense_plain_type;
  using stan::math::test::all_eigen_dense_exprs;
  all_eigen_dense_exprs<true, false, false, false, is_eigen_dense_plain_type,
                        Matrix, -1, -1>();
  all_eigen_dense_exprs<true, false, false, false, is_eigen_dense_plain_type,
                        Matrix, 1, -1>();
  all_eigen_dense_exprs<true, false, false, false, is_eigen_dense_plain_type,
                        Matrix, -1, 1>();
  all_eigen_dense_exprs<true, false, false, false, is_eigen_dense_plain_type,
                        Array, -1, -1>();
  all_eigen_dense_exprs<true, false, false, false, is_eigen_dense_plain_type,
                        Array, 1, -1>();
  all_eigen_dense_exprs<true, false, false, false, is_eigen_dense_plain_type,
                        Array, -1, 1>();
}
