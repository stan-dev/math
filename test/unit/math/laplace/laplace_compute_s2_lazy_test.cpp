#include <stan/math/mix.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <ostream>
#include <tuple>
#include <vector>

namespace {

// test matrix without zero entries
Eigen::MatrixXd dense_test_matrix(int n, double shift) {
  Eigen::MatrixXd m(n, n);
  for (int c = 0; c < n; ++c) {
    for (int r = 0; r < n; ++r) {
      m(r, c) = 1.0 / (shift + r + 2 * c) + (r == c ? 2.0 : 0.0);
    }
  }
  return m;
}

}  // namespace

// laplace_marginal_density() passes A to compute_s2() as an unevaluated
// expression holding a lazy product, because compute_s2() reads only the
// diagonal blocks of A: the result must match the one for the evaluated
// matrix, for each form of A that laplace_marginal_density() builds.
TEST(laplace_compute_s2_lazy, lazy_expression_matches_evaluated_matrix) {
  const int n = 6;
  const Eigen::MatrixXd K = dense_test_matrix(n, 1.0);
  const Eigen::MatrixXd C = dense_test_matrix(n, 2.5);
  const Eigen::MatrixXd D = dense_test_matrix(n, 4.0);
  Eigen::VectorXd theta(n);
  theta << 0.3, -0.2, 0.5, 1.0, -0.7, 0.1;
  const std::vector<int> y{1, 0, 2, 3, 0, 1, 2};
  const std::vector<int> y_index{1, 2, 3, 4, 5, 6, 3};
  const Eigen::VectorXd mean = Eigen::VectorXd::Constant(n, 0.4);
  std::ostream* msgs = nullptr;
  auto s2 = [&](const auto& A, int block_size) {
    return Eigen::VectorXd(stan::math::laplace_likelihood::compute_s2(
        stan::math::poisson_log_likelihood{}, theta, A, block_size,
        std::forward_as_tuple(y, y_index, mean), msgs));
  };
  for (int block_size : {1, 2, 3}) {
    // solver 1: A = K - C' C
    const Eigen::MatrixXd solver_1_dense = K - C.transpose() * C;
    EXPECT_MATRIX_FLOAT_EQ(s2(solver_1_dense, block_size),
                           s2(K - C.transpose().lazyProduct(C), block_size));
    // solver 2: A = C' C
    const Eigen::MatrixXd solver_2_dense = C.transpose() * C;
    EXPECT_MATRIX_FLOAT_EQ(s2(solver_2_dense, block_size),
                           s2(C.transpose().lazyProduct(C), block_size));
    // solver 3: A = K - C D
    const Eigen::MatrixXd solver_3_dense = K - C * D;
    EXPECT_MATRIX_FLOAT_EQ(s2(solver_3_dense, block_size),
                           s2(K - C.lazyProduct(D), block_size));
  }
}
