#include <iostream>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/fun/to_soa_sparse_matrix.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/util.hpp>
#include <vector>

TEST_F(AgradRev, to_soa_sparse_matrix_matrix_double) {
  using stan::math::to_soa_sparse_matrix;
  using stan::math::var;
  using stan::math::var_value;
  std::vector<int> v{0, 1, 2, 0, 1};
  std::vector<int> u{0, 1, 2, 3, 4, 5};
  Eigen::VectorXd w(5);
  int m = 5;
  int n = 5;
  w << 1, 2, 3, 4, 5;
  var_value<Eigen::SparseMatrix<double, Eigen::RowMajor>> w_mat_arena
      = to_soa_sparse_matrix<Eigen::RowMajor>(m, n, w, u, v);
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(w_mat_arena.val().valuePtr()[i], w(i));
  }
}

TEST_F(AgradRev, to_soa_sparse_matrix_matrix_var) {
  using stan::math::to_soa_sparse_matrix;
  using stan::math::var;
  using stan::math::var_value;
  std::vector<int> v{0, 1, 2, 0, 1};
  std::vector<int> u{0, 1, 2, 3, 4, 5};
  Eigen::VectorXd w(5);
  int m = 5;
  int n = 5;
  w << 1, 2, 3, 4, 5;
  Eigen::Matrix<var, -1, 1> w_var(w);
  var_value<Eigen::SparseMatrix<double, Eigen::RowMajor>> w_mat_arena
      = to_soa_sparse_matrix<Eigen::RowMajor>(m, n, w_var, u, v);
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(w_mat_arena.val().valuePtr()[i], w.val()(i));
  }
}

TEST_F(AgradRev, to_soa_sparse_matrix_var_matrix) {
  using stan::math::to_soa_sparse_matrix;
  using stan::math::var;
  using stan::math::var_value;
  std::vector<int> v{0, 1, 2, 0, 1};
  std::vector<int> u{0, 1, 2, 3, 4, 5};
  Eigen::VectorXd w(5);
  int m = 5;
  int n = 5;
  w << 1, 2, 3, 4, 5;
  var_value<Eigen::VectorXd> w_var(w);
  var_value<Eigen::SparseMatrix<double, Eigen::RowMajor>> w_mat_arena
      = to_soa_sparse_matrix<Eigen::RowMajor>(m, n, w_var, u, v);
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(w_mat_arena.val().valuePtr()[i], w.val()(i));
  }
  // Changing this value should change the adjoint of the sparse matrix
  for (int i = 0; i < 5; ++i) {
    w_mat_arena.adj().valuePtr()[i] = i;
  }
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(w_mat_arena.adj().valuePtr()[i], w_var.adj()(i));
  }
  for (int i = 0; i < 5; ++i) {
    w_mat_arena.adj().valuePtr()[i] = 0;
  }
  for (int i = 0; i < 5; ++i) {
    w_var.adj().coeffRef(i) = i;
  }
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(w_mat_arena.adj().valuePtr()[i], w_var.adj()(i));
  }
}
