#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/laplace/roach_data/y.hpp>
#include <test/unit/math/laplace/roach_data/sigmaz.hpp>
#include <test/unit/math/laplace/csv_reader.hpp>

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
    auto&& theta = stan::math::to_ref(theta_arg);
    auto&& mu = stan::math::to_ref(mu_arg);
    auto mu_theta = stan::math::eval(stan::math::add(stan::math::as_column_vector_or_scalar(mu), theta));
    return stan::math::poisson_log_lpmf<false>(y_arg, mu_theta);
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
struct integrand_vec_functor {
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
    auto sigma = phi[0];
    auto mu = stan::math::as_column_vector_or_scalar(phi).tail(theta.size()).eval();
    auto p = stan::math::exp((stan::math::normal_lpdf<false>(theta, 0, sigma) +
          stan::math::poisson_log_lpmf<false>(y_i, stan::math::add(theta, mu))));
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

TEST(WriteArrayBodySimple, ExceededIteration) {
  const double integrate_1d_reltol = 1e-8;
  auto mu_samples = stan::math::test::laplace::read_matrix_csv("./test/unit/math/laplace/roach_data/mu_bad.csv");
  auto sigmaz_samples = stan::math::test::laplace::read_matrix_csv("./test/unit/math/laplace/roach_data/sigma_bad.csv");
  auto y_samples_dbl = stan::math::test::laplace::read_matrix_csv("./test/unit/math/laplace/roach_data/y_bad.csv");
  auto y_samples = y_samples_dbl.cast<int>();
  std::cout << "y(" << y_samples.rows() << ", " << y_samples.cols() << ")\n";
  std::cout << "sigmaz(" << sigmaz_samples.rows() << ", " << sigmaz_samples.cols() << ")\n";
  std::cout << "mu(" << mu_samples.rows() << ", " << mu_samples.cols() << ")\n";
  const int num_samples = mu_samples.cols();
  const int  N = mu_samples.rows();
  std::ostream* pstream = nullptr;
  for (int i = 1; i <= N; ++i) {
    auto y = y_samples(i-1, 0);
    auto mu = mu_samples(i-1, 0);
    auto sigmaz = sigmaz_samples(i-1, 0);

    std::cout << "(i, y, mu, sigma)  = (" << i << ", "
              << y << ", " << mu << ", " << sigmaz << ")" << std::endl;
    double ll_laplace_val{0};
    try {
        ll_laplace_val = stan::math::laplace_marginal(
          poisson_re_log_ll_functor(),
          std::forward_as_tuple(y, mu),
          cov_fun_functor(),
          std::tuple<double, int>(sigmaz, 1),
          pstream);
    } catch (const std::domain_error& e) {
        // Log bad values to CSV files
        ADD_FAILURE() << "Laplace failed" << "(i, y, mu, sigma)  = (" << i << ", "
              << y << ", " << mu << ", " << sigmaz << ")" << "\nerror: " << e.what();
        continue;
    }
    double piece{0};
    try {
      piece = stan::math::integrate_1d(
        integrand_functor(),
        stan::math::negative_infinity(),
        stan::math::positive_infinity(),
        std::vector<double>{sigmaz, mu},
        std::vector<double>{0},
        std::vector<int>{ y },
        pstream,
        integrate_1d_reltol
      );
      EXPECT_NEAR(ll_laplace_val, std::log(piece), 8e-2) <<
        "for (i) = (" << i << "), laplace and integrated results should be close";
    } catch (const std::domain_error& e) {
      std::cout << "Integration Failed: y and mu for i = " << i << ": ("
                << y << ", " << mu << ")" << std::endl;
      std::cout << "Failed: " << e.what() << std::endl;
      continue;
    }
  }
}

TEST(WriteArrayBodySimple, ExecutesBodyWithHardcodedData) {
  const double integrate_1d_reltol = 1e-8;
  auto&& y = stan::math::test::roaches::y;
  auto&& sigmaz_samples = stan::math::test::roaches::sigmaz;
  auto mu_samples = stan::math::test::laplace::read_matrix_csv("./test/unit/math/laplace/roach_data/mu.csv");
  auto ll_samples = stan::math::test::laplace::read_matrix_csv("./test/unit/math/laplace/roach_data/integrated_likelihood.csv");
  std::cout << "y(" << y.size() << ")\n";
  std::cout << "sigmaz(" << sigmaz_samples.rows() << ", " << sigmaz_samples.cols() << ")\n";
  std::cout << "mu(" << mu_samples.rows() << ", " << mu_samples.cols() << ")\n";
  std::cout << "ll_samples(" << ll_samples.rows() << ", " << ll_samples.cols() << ")\n";
  const int num_samples = mu_samples.cols();
  const int  N = mu_samples.rows();
  std::ostream* pstream = nullptr;
  for (int iter = 0; iter < num_samples; ++iter) {
    std::vector<double> ll_laplace_vec;
    double ll_integrate_1d = 0;
    double ll_laplace = 0;
    std::vector<double> ll_integrate_1d_vec;
    auto mu = mu_samples.col(iter);
    auto sigmaz = sigmaz_samples(0, iter);
    for (int i = 1; i <= N; ++i) {
      //      std::cout << "y and mu for (i, iter) = (" << i << ", " << iter << "): ("
      //                << y[i - 1] << ", " << mu[i - 1] << ")" << std::endl;
      double ll_laplace_val{0};
      try {
          ll_laplace_val = stan::math::laplace_marginal(
            poisson_re_log_ll_functor(),
            std::forward_as_tuple(y[i - 1], mu[i - 1]),
            cov_fun_functor(),
            std::tuple<double, int>(sigmaz, 1),
            pstream);
      } catch (const std::domain_error& e) {
          // Log bad values to CSV files

/*
          std::ofstream y_bad("./test/unit/math/laplace/roach_data/y_bad.csv", std::ios::app);
          std::ofstream mu_bad("./test/unit/math/laplace/roach_data/mu_bad.csv", std::ios::app);
          std::ofstream sigma_bad("./test/unit/math/laplace/roach_data/sigma_bad.csv", std::ios::app);
          if (y_bad && mu_bad && sigma_bad) {
            y_bad << y[i - 1] << '\n';
            mu_bad << mu[i - 1] << '\n';
            sigma_bad << sigmaz << '\n';
          }
*/
          ADD_FAILURE() << "LAPLACE FAILURE: y and mu for i = " << i << ": ("
                    << y[i - 1] << ", " << mu[i - 1] << ")" << "\nerror: " << e.what() << std::endl;
          continue;
      }
      double piece{0};
      try {
        piece = stan::math::integrate_1d(
          integrand_functor(),
          stan::math::negative_infinity(),
          stan::math::positive_infinity(),
          std::vector<double>{sigmaz, mu[i-1]},
          std::vector<double>{0},
          std::vector<int>{ y[i-1] },
          pstream,
          integrate_1d_reltol
        );
        ll_laplace_vec.push_back(ll_laplace_val);
        ll_integrate_1d_vec.push_back(std::log(piece));
        ll_integrate_1d += std::log(piece);
        ll_laplace += ll_laplace_val;
        EXPECT_NEAR(ll_laplace_val, std::log(piece), 8e-2) <<
         "for (i, iter) = (" << i << ", " << iter << "), laplace and integrated results should be close";
      } catch (const std::domain_error& e) {
        std::cout << "Integration Failed: y and mu for i = " << i << ": ("
                  << y[i - 1] << ", " << mu[i - 1] << ")" << std::endl;
        std::cout << "Failed: " << e.what() << std::endl;
        continue;
      }
    }
    auto ll_laplace_all = stan::math::laplace_marginal(
      poisson_re_log_ll_functor(),
      std::forward_as_tuple(y, mu),
      cov_fun_functor(),
      std::tuple<double, int>(sigmaz, N),
      pstream);
    for (int i = 0; i < ll_integrate_1d_vec.size(); ++i) {
      break;
      std::cout << "results for y,mu[" << i << "]: ("
                << y[i] << ", " << mu[i] << "): \n"
                << "\tlaplace:    " << ll_laplace_vec[i] << "\n"
                << "\tintegrated: " << ll_integrate_1d_vec[i] << "\n"
                << "\tdifference: " << ll_laplace_vec[i] - ll_integrate_1d_vec[i] << std::endl;
    }
    // Assertions
    EXPECT_NEAR(ll_laplace, ll_integrate_1d, 2)
        << "For iter " << iter << ", Laplace and integrated results should be close";
    EXPECT_TRUE(std::isfinite(ll_laplace)) << "Laplace result should be finite";
    EXPECT_TRUE(std::isfinite(ll_integrate_1d)) << "Integrated result should be finite";
  }
}

