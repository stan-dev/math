#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <vector>

namespace {

std::size_t tape_size() {
  const auto& stack = *stan::math::ChainableStack::instance_;
  return stack.var_stack_.size() + stack.var_nochain_stack_.size();
}

// Nodes put on the reverse-mode tape by one forward-mode (fvar<var>) sweep
// of a GLM log density, the sweep hessian_times_vector runs for every
// direction: the design-matrix products must not add one node per matrix
// coefficient, so the count has to be O(n), not O(n * k).
template <typename F>
std::size_t nodes_per_sweep(F&& f, int n, int k) {
  using stan::math::fvar;
  using stan::math::var;
  stan::math::nested_rev_autodiff nested;
  const Eigen::MatrixXd x = Eigen::MatrixXd::Random(n, k);
  Eigen::Matrix<fvar<var>, -1, 1> beta(k);
  for (int j = 0; j < k; ++j) {
    beta(j) = fvar<var>(var(0.01 * j), 1.0);
  }
  Eigen::Matrix<fvar<var>, -1, 1> alpha(n);
  for (int i = 0; i < n; ++i) {
    alpha(i) = fvar<var>(var(0.1), 1.0);
  }
  const std::size_t before = tape_size();
  fvar<var> lp = f(x, alpha, beta);
  EXPECT_TRUE(std::isfinite(lp.val_.val()));
  return tape_size() - before;
}

}  // namespace

TEST(mathMixProbGlm, forwardModeSweepAddsLinearNumberOfNodes) {
  const int n = 50, k = 30;
  std::vector<int> y_count(n, 3);
  std::vector<int> y_binary(n, 1);
  std::vector<int> y_trials(n, 5);
  Eigen::VectorXd y_real = Eigen::VectorXd::Constant(n, 0.5);
  // before this change a sweep added about 4 * n * k nodes for the products
  // plus ~20 * n elementwise ones; now only the elementwise ones remain
  const std::size_t bound = 40 * n;

  EXPECT_LT(nodes_per_sweep(
                [&](auto&& x, auto&& alpha, auto&& beta) {
                  return stan::math::neg_binomial_2_log_glm_lpmf(
                      y_count, x, alpha, beta, 5.0);
                },
                n, k),
            bound);
  EXPECT_LT(nodes_per_sweep(
                [&](auto&& x, auto&& alpha, auto&& beta) {
                  return stan::math::poisson_log_glm_lpmf(y_count, x, alpha,
                                                          beta);
                },
                n, k),
            bound);
  EXPECT_LT(nodes_per_sweep(
                [&](auto&& x, auto&& alpha, auto&& beta) {
                  return stan::math::bernoulli_logit_glm_lpmf(y_binary, x,
                                                              alpha, beta);
                },
                n, k),
            bound);
  EXPECT_LT(nodes_per_sweep(
                [&](auto&& x, auto&& alpha, auto&& beta) {
                  return stan::math::binomial_logit_glm_lpmf(y_count, y_trials,
                                                             x, alpha, beta);
                },
                n, k),
            bound);
  EXPECT_LT(nodes_per_sweep(
                [&](auto&& x, auto&& alpha, auto&& beta) {
                  return stan::math::normal_id_glm_lpdf(y_real, x, alpha, beta,
                                                        1.5);
                },
                n, k),
            bound);
}
