#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace {

using stan::math::bernoulli_logit_likelihood;
using stan::math::neg_binomial_2_log_likelihood;
using stan::math::neg_binomial_2_log_likelihood_summary;
using stan::math::poisson_log_likelihood;
using stan::math::laplace_likelihood::analytic_derivatives;
namespace lli = stan::math::laplace_likelihood::internal;

// The likelihood `LL` behind a lambda: no `analytic_derivatives`
// specialization, so the derivative functions use nested autodiff.
template <typename LL>
auto autodiff_version() {
  return [](auto&&... args) { return LL{}(args...); };
}

// relative tolerance for entries above 1 in absolute value,
// absolute tolerance below 1
void expect_near_rel(const Eigen::VectorXd& expected,
                     const Eigen::VectorXd& actual, double tol,
                     const char* what) {
  ASSERT_EQ(expected.size(), actual.size()) << what;
  for (Eigen::Index i = 0; i < expected.size(); ++i) {
    EXPECT_TRUE(std::isfinite(actual(i))) << what << " [" << i << "]";
    EXPECT_NEAR(expected(i), actual(i),
                tol * std::fmax(1.0, std::fabs(expected(i))))
        << what << " [" << i << "]";
  }
}

void expect_same_sparse(const Eigen::SparseMatrix<double>& expected,
                        const Eigen::SparseMatrix<double>& actual, double tol,
                        const char* what) {
  ASSERT_EQ(expected.rows(), actual.rows()) << what;
  ASSERT_EQ(expected.cols(), actual.cols()) << what;
  EXPECT_EQ(expected.nonZeros(), actual.nonZeros()) << what;
  const Eigen::MatrixXd e = expected;
  const Eigen::MatrixXd a = actual;
  for (Eigen::Index c = 0; c < e.cols(); ++c) {
    for (Eigen::Index r = 0; r < e.rows(); ++r) {
      EXPECT_NEAR(e(r, c), a(r, c), tol * std::fmax(1.0, std::fabs(e(r, c))))
          << what << " (" << r << ", " << c << ")";
    }
  }
}

// Each derivative function of laplace_likelihood::internal, for the
// analytic version (`LL{}`) and for the autodiff version, with
// hessian_block_size 1 and 2.
template <typename LL, typename... Args>
void expect_analytic_equals_autodiff(const Eigen::VectorXd& theta, double tol,
                                     const Args&... args) {
  static_assert(lli::use_analytic_derivatives_v<LL, const Eigen::VectorXd&,
                                                const Args&...>,
                "the analytic derivatives must be selected");
  auto ad = autodiff_version<LL>();
  std::ostream* msgs = nullptr;
  expect_near_rel(lli::theta_grad(ad, theta, msgs, args...),
                  lli::theta_grad(LL{}, theta, msgs, args...), tol,
                  "theta_grad");
  expect_near_rel(lli::diagonal_hessian(ad, theta, msgs, args...),
                  lli::diagonal_hessian(LL{}, theta, msgs, args...), tol,
                  "diagonal_hessian");
  for (int block_size : {1, 2}) {
    expect_same_sparse(
        lli::block_hessian(ad, theta, block_size, msgs, args...),
        lli::block_hessian(LL{}, theta, block_size, msgs, args...), tol,
        "block_hessian");
  }
  expect_near_rel(lli::third_diff(ad, theta, msgs, args...),
                  lli::third_diff(LL{}, theta, msgs, args...), tol,
                  "third_diff");
}

// 4 groups, 1 to 3 observations per group
const std::vector<int> y_counts{0, 3, 7, 2, 1, 4, 12, 9};
const std::vector<int> y_binary{0, 1, 1, 0, 1, 1, 0, 1};
const std::vector<int> y_index{1, 1, 1, 2, 3, 3, 4, 4};

// a likelihood with the formulas of poisson_log_likelihood that counts the
// calls of its analytic derivatives
struct counted_poisson_likelihood {
  template <typename... Args>
  auto operator()(Args&&... args) const {
    return poisson_log_likelihood()(std::forward<Args>(args)...);
  }
};
int n_analytic_calls = 0;

}  // namespace

namespace stan {
namespace math {
namespace laplace_likelihood {
template <>
struct analytic_derivatives<counted_poisson_likelihood> : std::true_type {
  template <typename... Args>
  static Eigen::VectorXd theta_grad(const Args&... args) {
    ++n_analytic_calls;
    return analytic_derivatives<poisson_log_likelihood>::theta_grad(args...);
  }
  template <typename... Args>
  static Eigen::VectorXd neg_hessian_diag(const Args&... args) {
    ++n_analytic_calls;
    return analytic_derivatives<poisson_log_likelihood>::neg_hessian_diag(
        args...);
  }
  template <typename... Args>
  static Eigen::VectorXd third_diff(const Args&... args) {
    ++n_analytic_calls;
    return analytic_derivatives<poisson_log_likelihood>::third_diff(args...);
  }
};
}  // namespace laplace_likelihood
}  // namespace math
}  // namespace stan

namespace {

struct se_covariance {
  template <typename T1, typename T2>
  auto operator()(const std::vector<double>& x, const T1& sigma, const T2& rho,
                  std::ostream*) const {
    return stan::math::add_diag(stan::math::gp_exp_quad_cov(x, sigma, rho),
                                1e-8);
  }
};
const std::vector<double> x_gp{0.0, 1.3, 2.1, 3.7};

}  // namespace

TEST(laplace_analytic_derivatives, trait) {
  EXPECT_TRUE(analytic_derivatives<poisson_log_likelihood>::value);
  EXPECT_TRUE(analytic_derivatives<bernoulli_logit_likelihood>::value);
  EXPECT_TRUE(analytic_derivatives<neg_binomial_2_log_likelihood>::value);
  EXPECT_TRUE(
      analytic_derivatives<neg_binomial_2_log_likelihood_summary>::value);
  auto ad = autodiff_version<poisson_log_likelihood>();
  EXPECT_FALSE(analytic_derivatives<decltype(ad)>::value);
  // a var argument selects the autodiff version
  EXPECT_FALSE((lli::use_analytic_derivatives_v<
                neg_binomial_2_log_likelihood, Eigen::VectorXd, stan::math::var,
                std::vector<int>, std::vector<int>, Eigen::VectorXd>));
}

TEST(laplace_analytic_derivatives, used_by_laplace_marginal) {
  const Eigen::VectorXd mean = Eigen::VectorXd::Constant(4, 0.5);
  std::ostream* msgs = nullptr;
  n_analytic_calls = 0;
  double lp = stan::math::laplace_marginal<false>(
      counted_poisson_likelihood{},
      std::forward_as_tuple(y_counts, y_index, mean), 1, se_covariance{},
      std::forward_as_tuple(x_gp, 1.0, 1.5), msgs);
  EXPECT_TRUE(std::isfinite(lp));
  EXPECT_GT(n_analytic_calls, 0);
}

TEST(laplace_analytic_derivatives, equal_to_autodiff) {
  Eigen::VectorXd theta(4);
  theta << 0.3, -1.2, 1.7, 0.05;
  const Eigen::VectorXd mean = Eigen::VectorXd::Constant(4, 0.4);
  const double tol = 1e-12;
  for (double shift : {0.0, -30.0, 30.0}) {
    const Eigen::VectorXd th = theta.array() + shift;
    expect_analytic_equals_autodiff<poisson_log_likelihood>(th, tol, y_counts,
                                                            y_index, mean);
    expect_analytic_equals_autodiff<bernoulli_logit_likelihood>(
        th, tol, y_binary, y_index, mean);
    expect_analytic_equals_autodiff<neg_binomial_2_log_likelihood>(
        th, tol, 2.5, y_counts, y_index, mean);
    Eigen::VectorXd eta(4);
    eta << 0.7, 2.5, 11.0, 1.3;
    expect_analytic_equals_autodiff<neg_binomial_2_log_likelihood>(
        th, tol, eta, y_counts, y_index, mean);
    const std::vector<int> n_per_group{3, 1, 2, 2};
    const std::vector<int> counts_per_group{10, 2, 5, 21};
    expect_analytic_equals_autodiff<neg_binomial_2_log_likelihood_summary>(
        th, tol, 2.5, y_counts, n_per_group, counts_per_group, mean);
  }
}

// The shifted cases of equal_to_autodiff compare entries of about 1e-13 with
// an absolute tolerance of 1e-12, so they do not see the cancellation in
// `1 - p`.  Reference values (mpmath, 50 digits, from the double values of
// `theta + mean` and `theta + mean - log(eta)`) check the stable forms.
TEST(laplace_analytic_derivatives, tail_reference_values) {
  Eigen::VectorXd theta(4);
  theta << 30.3, 28.8, 31.7, 30.05;
  const Eigen::VectorXd mean = Eigen::VectorXd::Constant(4, 0.4);
  std::ostream* msgs = nullptr;
  auto expect_rel = [](const std::vector<double>& ref,
                       const Eigen::VectorXd& actual, const char* what) {
    ASSERT_EQ(static_cast<Eigen::Index>(ref.size()), actual.size()) << what;
    for (int i = 0; i < 4; ++i) {
      EXPECT_NEAR(ref[i], actual(i), 1e-14 * std::fabs(ref[i]))
          << what << " [" << i << "]";
    }
  };
  expect_rel({1.3940574134239623e-13, 2.0825772910546344e-13,
              2.2918021714044091e-14, 1.1933367674388607e-13},
             lli::diagonal_hessian(bernoulli_logit_likelihood{}, theta, msgs,
                                   y_binary, y_index, mean),
             "bernoulli neg_hessian");
  expect_rel({1.3940574134238327e-13, 2.082577291053767e-13,
              2.2918021714043565e-14, 1.1933367674387183e-13},
             lli::third_diff(bernoulli_logit_likelihood{}, theta, msgs,
                             y_binary, y_index, mean),
             "bernoulli third_diff");
  expect_rel({2.0330003945763301e-12, 2.3428994524350021e-12,
              2.8647527142554155e-13, 3.8783444941756066e-12},
             lli::diagonal_hessian(neg_binomial_2_log_likelihood{}, theta, msgs,
                                   2.5, y_counts, y_index, mean),
             "neg_binomial_2 neg_hessian");
  expect_rel({2.0330003945758578e-12, 2.3428994524325625e-12,
              2.8647527142552514e-13, 3.8783444941744496e-12},
             lli::third_diff(neg_binomial_2_log_likelihood{}, theta, msgs, 2.5,
                             y_counts, y_index, mean),
             "neg_binomial_2 third_diff");
}

TEST(laplace_analytic_derivatives, finite_at_large_theta) {
  // exp(theta) overflows here; the inv_logit forms stay finite
  Eigen::VectorXd theta(4);
  theta << -800, 800, -750, 760;
  const Eigen::VectorXd mean = Eigen::VectorXd::Zero(4);
  expect_analytic_equals_autodiff<bernoulli_logit_likelihood>(
      theta, 1e-12, y_binary, y_index, mean);
  expect_analytic_equals_autodiff<neg_binomial_2_log_likelihood>(
      theta, 1e-12, 2.5, y_counts, y_index, mean);
}

TEST(laplace_analytic_derivatives, end_to_end_equal_to_autodiff) {
  using stan::math::var;
  const Eigen::VectorXd mean = Eigen::VectorXd::Constant(4, 0.5);
  std::ostream* msgs = nullptr;
  const double tol = 1e-7;
  for (int solver : {1, 2, 3}) {
    for (int block_size : {1, 2}) {
      auto ops = std::make_tuple(Eigen::VectorXd::Zero(4).eval(), 1e-10, 200,
                                 solver, 1000, false);
      // value and gradient w.r.t. sigma, rho (and eta) of `lpmf`
      auto run = [&](auto&& lpmf) {
        var sigma = 1.0, rho = 1.5, eta = 2.5;
        var lp = lpmf(sigma, rho, eta);
        stan::math::grad(lp.vi_);
        std::vector<double> out{lp.val(), sigma.adj(), rho.adj(), eta.adj()};
        stan::math::recover_memory();
        return out;
      };
      auto compare = [&](const std::vector<double>& a,
                         const std::vector<double>& b, const char* what) {
        for (size_t i = 0; i < a.size(); ++i) {
          EXPECT_NEAR(a[i], b[i], tol * std::fmax(1.0, std::fabs(a[i])))
              << what << " solver " << solver << " block size " << block_size
              << " output " << i;
        }
      };
      compare(run([&](auto s, auto r, auto) {
                return stan::math::laplace_marginal_tol<false>(
                    autodiff_version<poisson_log_likelihood>(),
                    std::forward_as_tuple(y_counts, y_index, mean), block_size,
                    se_covariance{}, std::forward_as_tuple(x_gp, s, r), ops,
                    msgs);
              }),
              run([&](auto s, auto r, auto) {
                return stan::math::laplace_marginal_tol_poisson_log_lpmf(
                    y_counts, y_index, mean, block_size, se_covariance{},
                    std::forward_as_tuple(x_gp, s, r), ops, msgs);
              }),
              "poisson");
      compare(run([&](auto s, auto r, auto) {
                return stan::math::laplace_marginal_tol<false>(
                    autodiff_version<bernoulli_logit_likelihood>(),
                    std::forward_as_tuple(y_binary, y_index, mean), block_size,
                    se_covariance{}, std::forward_as_tuple(x_gp, s, r), ops,
                    msgs);
              }),
              run([&](auto s, auto r, auto) {
                return stan::math::laplace_marginal_tol_bernoulli_logit_lpmf(
                    y_binary, y_index, mean, block_size, se_covariance{},
                    std::forward_as_tuple(x_gp, s, r), ops, msgs);
              }),
              "bernoulli");
      // eta as a var: the Newton iterations use the analytic derivatives,
      // the gradient w.r.t. eta the autodiff ones
      compare(run([&](auto s, auto r, auto e) {
                return stan::math::laplace_marginal_tol<false>(
                    autodiff_version<neg_binomial_2_log_likelihood>(),
                    std::forward_as_tuple(e, y_counts, y_index, mean),
                    block_size, se_covariance{},
                    std::forward_as_tuple(x_gp, s, r), ops, msgs);
              }),
              run([&](auto s, auto r, auto e) {
                return stan::math::laplace_marginal_tol_neg_binomial_2_log_lpmf(
                    y_counts, y_index, e, mean, block_size, se_covariance{},
                    std::forward_as_tuple(x_gp, s, r), ops, msgs);
              }),
              "neg_binomial_2 (eta var)");
    }
  }
}
