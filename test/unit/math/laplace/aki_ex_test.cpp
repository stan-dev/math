#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/laplace/roach_data/sigmaz.hpp>
#include <test/unit/math/laplace/roach_data/y.hpp>
#include <test/unit/math/laplace/csv_reader.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

struct poisson_re_log_ll_functor {
  template <typename T0, typename T1, typename T2>
  stan::return_type_t<stan::base_type_t<T0>, stan::base_type_t<T2>> operator()(
      const T0& theta_arg, const T1& y_arg, const T2& mu_arg,
      std::ostream* pstream) const {
    auto&& theta = stan::math::to_ref(theta_arg);
    auto&& mu = stan::math::to_ref(mu_arg);
    auto mu_theta = stan::math::eval(
        stan::math::add(stan::math::as_column_vector_or_scalar(mu), theta));
    return stan::math::poisson_log_lpmf<false>(y_arg, mu_theta);
  }
};

struct integrand_functor {
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  stan::return_type_t<T0, T1, stan::base_type_t<T2>, stan::base_type_t<T3>>
  operator()(const T0& theta, const T1& notused, const T2& phi, const T3& X_i,
             const T4& y_i, std::ostream* pstream) const {
    using local_scalar_t = stan::return_type_t<T0, T1, stan::base_type_t<T2>,
                                               stan::base_type_t<T3>>;
    // suppress unused var warning
    local_scalar_t sigma = phi[0];
    local_scalar_t mu = phi[1];
    local_scalar_t p = stan::math::exp(
        (stan::math::normal_lpdf<false>(theta, 0, sigma)
         + stan::math::poisson_log_lpmf<false>(y_i, (theta + mu))));
    return p;
  }
};
struct integrand_vec_functor {
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  stan::return_type_t<T0, T1, stan::base_type_t<T2>, stan::base_type_t<T3>>
  operator()(const T0& theta, const T1& notused, const T2& phi, const T3& X_i,
             const T4& y_i, std::ostream* pstream) const {
    // suppress unused var warning
    auto sigma = phi[0];
    auto mu
        = stan::math::as_column_vector_or_scalar(phi).tail(theta.size()).eval();
    auto p = stan::math::exp((stan::math::normal_lpdf<false>(theta, 0, sigma)
                              + stan::math::poisson_log_lpmf<false>(
                                  y_i, stan::math::add(theta, mu))));
    return p;
  }
};
struct cov_fun_functor {
  template <typename T0>
  Eigen::Matrix<stan::return_type_t<T0>, -1, -1> operator()(
      const T0& sigma, const int& N, std::ostream* pstream) const {
    return stan::math::diag_matrix(
        stan::math::rep_vector(stan::math::pow(sigma, 2), N));
  }
};

TEST(WriteArrayBodySimple, ExceededIteration) {
  stan::test::relative_tolerance rel_tol(5e-2);
  const double integrate_1d_reltol = 1e-8;
  auto mu_samples = stan::math::test::laplace::read_matrix_csv(
      "./test/unit/math/laplace/roach_data/mu_bad.csv");
  auto sigmaz_samples = stan::math::test::laplace::read_matrix_csv(
      "./test/unit/math/laplace/roach_data/sigma_bad.csv");
  auto y_samples_dbl = stan::math::test::laplace::read_matrix_csv(
      "./test/unit/math/laplace/roach_data/y_bad.csv");
  auto y_samples = y_samples_dbl.cast<int>();
  const int N = mu_samples.rows();
  std::ostream* pstream = nullptr;
  for (int i = 1; i <= N; ++i) {
    auto y = y_samples(i - 1, 0);
    auto mu = mu_samples(i - 1, 0);
    auto sigmaz = sigmaz_samples(i - 1, 0);
    double ll_laplace_val{0};
    try {
      ll_laplace_val = stan::math::laplace_marginal(
          poisson_re_log_ll_functor(), std::forward_as_tuple(y, mu),
          1, cov_fun_functor(), std::tuple<double, int>(sigmaz, 1), pstream);
    } catch (const std::domain_error& e) {
      // Log bad values to CSV files
      ADD_FAILURE() << "Laplace failed"
                    << "(i, y, mu, sigma)  = (" << i << ", " << y << ", " << mu
                    << ", " << sigmaz << ")"
                    << "\nerror: " << e.what();
      continue;
    }
    double piece{0};
    try {
      piece = stan::math::integrate_1d(
          integrand_functor(), stan::math::negative_infinity(),
          stan::math::positive_infinity(), std::vector<double>{sigmaz, mu},
          std::vector<double>{0}, std::vector<int>{y}, pstream,
          integrate_1d_reltol);
      std::string msg = std::string("for (i) = (") + std::to_string(i)
                        + "), laplace and integrated results should be close";
      expect_near_rel(msg, ll_laplace_val, std::log(piece), rel_tol,
                      "laplace_val", "integrated_val");
    } catch (const std::domain_error& e) {
      // NOTE: Failures for integration our fine since we are testing laplace.
      continue;
    }
  }
}

TEST(WriteArrayBodySimple, ExecutesBodyWithHardcodedData) {
  stan::test::relative_tolerance rel_tol(5e-1);
  const double integrate_1d_reltol = 1e-8;
  auto&& y = stan::math::test::roaches::y;
  auto&& sigmaz_samples = stan::math::test::roaches::sigmaz;
  auto mu_samples = stan::math::test::laplace::read_matrix_csv(
      "./test/unit/math/laplace/roach_data/mu.csv");
  const int num_samples = mu_samples.cols();
  const int N = mu_samples.rows();
  std::ostream* pstream = nullptr;
  for (int iter = 0; iter < num_samples; ++iter) {
    std::vector<double> ll_laplace_vec;
    double ll_integrate_1d = 0;
    double ll_laplace = 0;
    std::vector<double> ll_integrate_1d_vec;
    auto mu = mu_samples.col(iter);
    auto sigmaz = sigmaz_samples(0, iter);
    for (int i = 1; i <= N; ++i) {
      double ll_laplace_val{0};
      try {
        ll_laplace_val = stan::math::laplace_marginal(
            poisson_re_log_ll_functor(),
            std::forward_as_tuple(y[i - 1], mu[i - 1]), 1, cov_fun_functor(),
            std::tuple<double, int>(sigmaz, 1), pstream);
      } catch (const std::domain_error& e) {
        ADD_FAILURE() << "LAPLACE FAILURE: y and mu for i = " << i << ": ("
                      << y[i - 1] << ", " << mu[i - 1] << ")"
                      << "\nerror: " << e.what() << std::endl;
        continue;
      }
      double piece{0};
      try {
        piece = stan::math::integrate_1d(
            integrand_functor(), stan::math::negative_infinity(),
            stan::math::positive_infinity(),
            std::vector<double>{sigmaz, mu[i - 1]}, std::vector<double>{0},
            std::vector<int>{y[i - 1]}, pstream, integrate_1d_reltol);
        ll_laplace_vec.push_back(ll_laplace_val);
        ll_integrate_1d_vec.push_back(std::log(piece));
        ll_integrate_1d += std::log(piece);
        ll_laplace += ll_laplace_val;
        std::string msg = std::string("for (i) = (") + std::to_string(i)
                          + "), laplace and integrated results should be close";
        expect_near_rel(msg, ll_laplace_val, std::log(piece), rel_tol,
                        "laplace_val", "integrated_val");
      } catch (const std::domain_error& e) {
        // Note: Integration failures are fine since we are testing laplace.
        continue;
      }
    }
    auto ll_laplace_all = stan::math::laplace_marginal(
        poisson_re_log_ll_functor(), std::forward_as_tuple(y, mu),
        1, cov_fun_functor(), std::tuple<double, int>(sigmaz, N), pstream);
    stan::test::relative_tolerance sum_rel_tol(3e-2);
    expect_near_rel("sum laplace vs integrated sum", ll_laplace,
                    ll_integrate_1d, sum_rel_tol, "laplace_sum",
                    "integrated_sum");
    expect_near_rel("total laplace vs integrated sum", ll_laplace_all,
                    ll_integrate_1d, sum_rel_tol, "laplace_sum",
                    "integrated_sum");
    EXPECT_TRUE(std::isfinite(ll_laplace)) << "Laplace result should be finite";
    EXPECT_TRUE(std::isfinite(ll_integrate_1d))
        << "Integrated result should be finite";
  }
}
