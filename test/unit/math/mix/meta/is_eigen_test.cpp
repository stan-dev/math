#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>


TEST(MathMetaPrim, eigen_tests) {
  using stan::math::var;
  using stan::math::fvar;
  using stan::is_eigen;
  /**
   * Dense Flags:
   * eigen_base_v, eigen_dense_v, eigen_matrix_v, eigen_map_v, eigen_mat_exp_v,
   *  eigen_dense_solver_v
   *
   * Sparse Flags:
   * eigen_base_v, eigen_sparse_compressed_v, eigen_sparse_matrix_v,
   * eigen_sparse_map_v
   */
  test_eigen_dense<true, true, true, true, true, true, double, is_eigen>();
  test_eigen_sparse<true, true, true, true, double, is_eigen>();
  test_eigen_dense<true, true, true, true, true, true, var, is_eigen>();
  test_eigen_sparse<true, true, true, true, var, is_eigen>();
  test_eigen_dense<true, true, true, true, true, true, fvar<double>, is_eigen>();
  test_eigen_sparse<true, true, true, true, fvar<double>, is_eigen>();
  test_eigen_dense<true, true, true, true, true, true, fvar<var>, is_eigen>();
  test_eigen_sparse<true, true, true, true, fvar<var>, is_eigen>();
  test_eigen_dense<true, true, true, true, true, true, fvar<fvar<var>>, is_eigen>();
  test_eigen_sparse<true, true, true, true, fvar<fvar<var>>, is_eigen>();
}
