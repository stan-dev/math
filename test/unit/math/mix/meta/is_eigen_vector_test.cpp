#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_eigen_vector_tests) {
  using stan::is_eigen_vector;
  using stan::math::test::all_eigen_dense_array;
  using stan::math::test::all_eigen_dense_decomp;
  using stan::math::test::all_eigen_dense_matrix;
  using stan::math::test::all_eigen_sparse;
  all_eigen_dense_matrix<false, false, false, false, false, is_eigen_vector, -1,
                         -1>();
  all_eigen_dense_matrix<false, true, true, true, true, is_eigen_vector, 1,
                         -1>();
  all_eigen_dense_matrix<false, true, true, true, true, is_eigen_vector, -1,
                         1>();
  all_eigen_dense_array<false, false, false, false, false, is_eigen_vector, -1,
                        -1>();
  all_eigen_dense_array<false, true, true, true, true, is_eigen_vector, 1,
                        -1>();
  all_eigen_dense_array<false, true, true, true, true, is_eigen_vector, -1,
                        1>();
  all_eigen_dense_decomp<false, is_eigen_vector>();
  all_eigen_sparse<false, false, false, false, is_eigen_vector>();
}
