#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <ostream>
#include <tuple>
#include <vector>

namespace {

// symmetric test matrix without zero entries
Eigen::MatrixXd dense_test_matrix(int n) {
  Eigen::MatrixXd m(n, n);
  for (int c = 0; c < n; ++c) {
    for (int r = 0; r < n; ++r) {
      m(r, c) = 1.0 / (1.0 + r + c) + (r == c ? 2.0 : 0.0);
    }
  }
  return m;
}

}  // namespace

TEST(laplace_block_diagonal_entries, keeps_only_the_diagonal_blocks) {
  const int n = 6;
  const Eigen::MatrixXd full = dense_test_matrix(n);
  for (int block_size : {1, 2, 3}) {
    const Eigen::MatrixXd blocks = stan::math::internal::block_diagonal_entries(
        n, block_size,
        [&](Eigen::Index r, Eigen::Index c) { return full(r, c); });
    ASSERT_EQ(n, blocks.rows());
    ASSERT_EQ(n, blocks.cols());
    for (int c = 0; c < n; ++c) {
      for (int r = 0; r < n; ++r) {
        if (r / block_size == c / block_size) {
          EXPECT_EQ(full(r, c), blocks(r, c));
        } else {
          EXPECT_EQ(0.0, blocks(r, c));
        }
      }
    }
  }
}

// laplace_marginal_density() passes only the diagonal blocks of A to
// compute_s2(), because compute_s2() reads no other entry of A: the result
// is the same for the full matrix and for its diagonal blocks.
TEST(laplace_block_diagonal_entries,
     compute_s2_reads_only_the_diagonal_blocks) {
  const int n = 6;
  const Eigen::MatrixXd full = dense_test_matrix(n);
  Eigen::VectorXd theta(n);
  theta << 0.3, -0.2, 0.5, 1.0, -0.7, 0.1;
  const std::vector<int> y{1, 0, 2, 3, 0, 1, 2};
  const std::vector<int> y_index{1, 2, 3, 4, 5, 6, 3};
  const Eigen::VectorXd mean = Eigen::VectorXd::Constant(n, 0.4);
  std::ostream* msgs = nullptr;
  for (int block_size : {1, 2, 3}) {
    const Eigen::MatrixXd blocks = stan::math::internal::block_diagonal_entries(
        n, block_size,
        [&](Eigen::Index r, Eigen::Index c) { return full(r, c); });
    const Eigen::VectorXd s2_full = stan::math::laplace_likelihood::compute_s2(
        stan::math::poisson_log_likelihood{}, theta, full, block_size,
        std::forward_as_tuple(y, y_index, mean), msgs);
    const Eigen::VectorXd s2_blocks
        = stan::math::laplace_likelihood::compute_s2(
            stan::math::poisson_log_likelihood{}, theta, blocks, block_size,
            std::forward_as_tuple(y, y_index, mean), msgs);
    ASSERT_EQ(n, s2_full.size());
    for (int i = 0; i < n; ++i) {
      EXPECT_EQ(s2_full(i), s2_blocks(i)) << "block size " << block_size;
    }
  }
}
