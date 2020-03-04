#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, eigen_dense_tests) {
  using stan::is_eigen_contiguous_map;
  using stan::math::test::all_eigen_dense_decomp;
  using stan::math::test::all_eigen_sparse;
  using stan::math::test::all_eigen_dense_matrix;
  using stan::math::test::all_eigen_dense_array;

  all_eigen_dense_matrix<false, false, false, false, true, is_eigen_contiguous_map, -1,
                              -1>();
  all_eigen_dense_matrix<false, false, false, false, true, is_eigen_contiguous_map, 1,
                              -1>();
  all_eigen_dense_matrix<false, false, false, false, true, is_eigen_contiguous_map, -1,
                              1>();
  all_eigen_dense_array<false, false, false, false, true, is_eigen_contiguous_map, -1,
                             -1>();
  all_eigen_dense_array<false, false, false, false, true, is_eigen_contiguous_map, 1,
                             -1>();
  all_eigen_dense_array<false, false, false, false, true, is_eigen_contiguous_map, -1,
                             1>();
  all_eigen_dense_decomp<false, is_eigen_contiguous_map>();
  all_eigen_sparse<false, false, false, false, is_eigen_contiguous_map>();
}
