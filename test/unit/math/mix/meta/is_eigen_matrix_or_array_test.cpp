#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

TEST(MathMetaPrim, eigen_matrix_or_array_tests) {
  using stan::is_eigen_matrix_or_array;
  using stan::math::test::all_eigen_dense_decomp;
  using stan::math::test::all_eigen_sparse;
  using stan::math::test::all_eigen_dense_matrix;
  using stan::math::test::all_eigen_dense_array;
  /**
   * Dense Flags:
   * eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v,
   *  eigen_dense_solver_v
   *
   * Sparse Flags:
   * eigen_base_v, eigen_sparse_compressed_v, eigen_sparse_matrix_v,
   * eigen_sparse_map_v
   */
  all_eigen_dense_matrix<false, true, true, true, true,
                              is_eigen_matrix_or_array, -1, -1>();
  all_eigen_dense_array<false, true, true, true, true,
                             is_eigen_matrix_or_array, -1, -1>();
  all_eigen_dense_decomp<true, is_eigen_matrix_or_array>();
  all_eigen_sparse<false, false, false, false, is_eigen_matrix_or_array>();
}
