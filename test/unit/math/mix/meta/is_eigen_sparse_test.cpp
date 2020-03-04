#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_eigen_sparse_tests) {
  using stan::is_eigen_sparse_matrix;
  using stan::math::test::all_eigen_dense_array;
  using stan::math::test::all_eigen_dense_decomp;
  using stan::math::test::all_eigen_dense_matrix;
  using stan::math::test::all_eigen_sparse;
  // Same for all compile time dimensions
  all_eigen_dense_matrix<false, false, false, false, false,
                         is_eigen_sparse_matrix, -1, -1>();
  all_eigen_dense_array<false, false, false, false, false,
                        is_eigen_sparse_matrix, -1, -1>();
  all_eigen_dense_decomp<false, is_eigen_sparse_matrix>();
  // Should this pass for eigenbase?
  all_eigen_sparse<false, true, true, true, is_eigen_sparse_matrix>();
}
