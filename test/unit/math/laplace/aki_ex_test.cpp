#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/laplace/aki_synth_data/x1.hpp>

#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>


struct poisson_re_log_ll_functor {
  template <typename T0, typename T1, typename T2>
  stan::return_type_t<stan::base_type_t<T0>, stan::base_type_t<T2>>
  operator()(const T0& theta_arg, const T1& y_arg, const T2& mu_arg,
             std::ostream* pstream) const {
    using local_scalar_t = stan::return_type_t<stan::base_type_t<T0>,
                              stan::base_type_t<T2>>;
    const auto& theta = stan::math::to_ref(theta_arg);
    const auto& mu = stan::math::to_ref(mu_arg);
    return stan::math::poisson_log_lpmf<false>(y_arg, stan::math::add(stan::math::as_column_vector_or_scalar(mu), theta));
  }
};
struct integrand_functor {
  template <typename T0, typename T1, typename T2, typename T3,
            typename T4>
  stan::return_type_t<T0, T1, stan::base_type_t<T2>,
    stan::base_type_t<T3>>
  operator()(const T0& theta, const T1& notused, const T2& phi,
             const T3& X_i, const T4& y_i, std::ostream* pstream) const {
    using local_scalar_t = stan::return_type_t<T0, T1,
                              stan::base_type_t<T2>,
                              stan::base_type_t<T3>>;
    // suppress unused var warning
    static constexpr bool propto = true;
    local_scalar_t sigma = phi[0];
    local_scalar_t mu = phi[1];
    local_scalar_t p = stan::math::exp((stan::math::normal_lpdf<false>(theta, 0, sigma) +
          stan::math::poisson_log_lpmf<false>(y_i, (theta + mu))));
    return p;
  }
};
struct cov_fun_functor {
  template <typename T0>
  Eigen::Matrix<stan::return_type_t<T0>,-1,-1>
  operator()(const T0& sigma, const int& N, std::ostream* pstream) const {
  using local_scalar_t = stan::return_type_t<T0>;
  return stan::math::diag_matrix(
             stan::math::rep_vector(stan::math::pow(sigma, 2), N));
  }
};

TEST(WriteArrayBodySimple, ExecutesBodyWithHardcodedData) {
  const std::vector<int>    y{183, 91, 171};
  const std::vector<double> mu{0.6, -0.2, 0.6};
  const int  N = y.size();
  const double sigmaz = 2.0;
  const double integrate_1d_reltol = 1e-6;
  std::ostream* pstream = nullptr;

  auto ll_laplace = stan::math::laplace_marginal(
    poisson_re_log_ll_functor(),
    std::forward_as_tuple(y, mu),
    cov_fun_functor(),
    std::tuple<double, int>(sigmaz, N),
    pstream);
  std::vector<double> ll_laplace_vec;
  for (int i = 1; i <= N; ++i) {
    auto ll_laplace_val = stan::math::laplace_marginal(
      poisson_re_log_ll_functor(),
      std::forward_as_tuple(y[i - 1], mu[i - 1]),
      cov_fun_functor(),
      std::tuple<double, int>(sigmaz, 1),
      pstream);
    ll_laplace_vec.push_back(ll_laplace_val);
  }
  double ll_integrate_1d = 0;
  std::vector<double> ll_integrate_1d_vec;
  for (int i = 1; i <= N; ++i) {
    double piece = stan::math::integrate_1d(
      integrand_functor(),
      stan::math::negative_infinity(),
      stan::math::positive_infinity(),
      std::vector<double>{sigmaz, mu[i-1]},
      std::vector<double>{0},
      std::vector<int>{ y[i-1] },
      pstream,
      integrate_1d_reltol
    );
    ll_integrate_1d_vec.push_back(std::log(piece));
    ll_integrate_1d += std::log(piece);
  }
  std::cout << "Laplace result: " << ll_laplace << std::endl;
  std::cout << "Integrated result: " << ll_integrate_1d << std::endl;
  for (int i = 0; i < ll_integrate_1d_vec.size(); ++i) {
    std::cout << "results for y,mu[" << i << "]: ("
              << y[i] << ", " << mu[i] << "): \n"
              << "\tlaplace:    " << ll_laplace_vec[i] << "\n"
              << "\tintegrated: " << ll_integrate_1d_vec[i] << "\n"
              << "\tdifference: " << ll_laplace_vec[i] - ll_integrate_1d_vec[i] << std::endl;
  }

  // Assertions
  EXPECT_NEAR(ll_laplace, ll_integrate_1d, 5e-2)
      << "Laplace and integrated results should be close";
  EXPECT_TRUE(std::isfinite(ll_laplace)) << "Laplace result should be finite";
  EXPECT_TRUE(std::isfinite(ll_integrate_1d)) << "Integrated result should be finite";
}

