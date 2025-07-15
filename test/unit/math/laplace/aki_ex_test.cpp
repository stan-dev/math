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


struct poisson_re_log_ll_functor__ {
  template <typename T0__, typename T2__>
  stan::return_type_t<stan::base_type_t<T0__>, stan::base_type_t<T2__>>
  operator()(const T0__& theta_arg__, const std::vector<int>& y_arg__, const T2__& mu_arg__,
             std::ostream* pstream__) const {
    using local_scalar_t__ = stan::return_type_t<stan::base_type_t<T0__>,
                              stan::base_type_t<T2__>>;
    const auto& theta = stan::math::to_ref(theta_arg__);
    const auto& mu = stan::math::to_ref(mu_arg__);
    return stan::math::poisson_log_lpmf<false>(y_arg__, stan::math::add(stan::math::as_column_vector_or_scalar(mu), theta));
  }
};
struct integrand_functor__ {
  template <typename T0__, typename T1__, typename T2__, typename T3__,
            typename T4__>
  stan::return_type_t<T0__, T1__, stan::base_type_t<T2__>,
    stan::base_type_t<T3__>>
  operator()(const T0__& theta, const T1__& notused, const T2__& phi,
             const T3__& X_i, const T4__& y_i, std::ostream* pstream__) const {
    using local_scalar_t__ = stan::return_type_t<T0__, T1__,
                              stan::base_type_t<T2__>,
                              stan::base_type_t<T3__>>;
    // suppress unused var warning
    static constexpr bool propto__ = true;
    local_scalar_t__ sigma = phi[0];
    local_scalar_t__ mu = phi[1];
    local_scalar_t__ p = stan::math::exp((stan::math::normal_lpdf<false>(theta, 0, sigma) +
          stan::math::poisson_log_lpmf<false>(y_i, (theta + mu))));
    return p;
  }
};
struct cov_fun_functor__ {
  template <typename T0__>
  Eigen::Matrix<stan::return_type_t<T0__>,-1,-1>
  operator()(const T0__& sigma, const int& N, std::ostream* pstream__) const {
  using local_scalar_t__ = stan::return_type_t<T0__>;
  return stan::math::diag_matrix(
             stan::math::rep_vector(stan::math::pow(sigma, 2), N));
  }
};

TEST(WriteArrayBodySimple, ExecutesBodyWithHardcodedData) {
  const int  N = 1;
  const std::vector<int>    y{183};
  const std::vector<double> mu{0.5};
  const double sigmaz = 2.5;
  const double integrate_1d_reltol = 1e-6;
  std::ostream* pstream__ = nullptr;
  auto ll_laplace = stan::math::laplace_marginal(
      poisson_re_log_ll_functor__(),
      std::tuple<const std::vector<int>&, const std::vector<double>&>(y, mu),
      cov_fun_functor__(),
      std::tuple<double, int>(sigmaz, N),
      pstream__);

  double ll_integrate_1d = 0;
  for (int i = 1; i <= N; ++i) {
    double piece = stan::math::integrate_1d(
      integrand_functor__(),
      stan::math::negative_infinity(),
      stan::math::positive_infinity(),
      std::vector<double>{sigmaz, mu[i-1]},
      std::vector<double>{0},
      std::vector<int>{ y[i-1] },
      pstream__,
      integrate_1d_reltol
    );
    ll_integrate_1d += std::log(piece);
  }
  std::cout << "Laplace result: " << ll_laplace << std::endl;
  std::cout << "Integrated result: " << ll_integrate_1d << std::endl;
  // --- end inlined body ---

  // Assertions
  EXPECT_TRUE(std::isfinite(ll_laplace)) << "Laplace result should be finite";
  EXPECT_TRUE(std::isfinite(ll_integrate_1d)) << "Integrated result should be finite";
}
