#include <limits>
#include <vector>
#include <gtest/gtest.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/tools/promotion.hpp>

#include <stan/math/rev.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/log1m.hpp>

template <typename T1, typename T2, typename T3, typename T4>
inline typename boost::math::tools::promote_args<T1, T2, T3, T4>::type
skew_de_ccdf_test(const T1& y, const T2& mu, const T3& sigma, const T4& tau) {
  using stan::math::log1m;
  using stan::math::log1m_exp;
  using std::exp;
  using std::log;

  if (y < mu) {
    return log1m_exp(log(tau) - 2 / sigma * (1 - tau) * (mu - y));
  } else {
    return log1m_exp(log1m((1 - tau) * exp(-2 / sigma * tau * (y - mu))));
  }
}

TEST(ProbDistributionsSkewedDoubleExponential,
     lccdf_computes_correct_gradients) {
  using stan::math::skew_double_exponential_lccdf;

  for (double ys : {-1.7, 0.2, 0.5, 0.9, 1.1, 3.2, 8.3}) {
    for (double mus : {-1.8, 0.1, 0.55, 0.89, 1.3, 4.2, 9.3}) {
      for (double sigmas : {0.1, 0.5, 1.1, 10.1}) {
        for (double taus : {0.01, 0.1, 0.5, 0.6, 0.9, 0.99}) {
          stan::math::var y = ys;
          stan::math::var mu = mus;
          stan::math::var sigma = sigmas;
          stan::math::var tau = taus;

          stan::math::var lp = skew_double_exponential_lccdf(y, mu, sigma, tau);
          std::vector<stan::math::var> theta;
          theta.push_back(y);
          theta.push_back(mu);
          theta.push_back(sigma);
          theta.push_back(tau);
          std::vector<double> grads;
          lp.grad(theta, grads);

          stan::math::var y_true = ys;
          stan::math::var mu_true = mus;
          stan::math::var sigma_true = sigmas;
          stan::math::var tau_true = taus;

          stan::math::var lp_test
              = skew_de_ccdf_test(y_true, mu_true, sigma_true, tau_true);
          std::vector<stan::math::var> theta_true;
          theta_true.push_back(y_true);
          theta_true.push_back(mu_true);
          theta_true.push_back(sigma_true);
          theta_true.push_back(tau_true);
          std::vector<double> grads_true;
          lp_test.grad(theta_true, grads_true);

          EXPECT_NEAR(grads_true[0], grads[0], 0.01);
          EXPECT_NEAR(grads_true[1], grads[1], 0.01);
          EXPECT_NEAR(grads_true[2], grads[2], 0.01);
          EXPECT_NEAR(grads_true[3], grads[3], 0.01);
        }
      }
    }
  }
}

TEST(ProbDistributionsSkewedDoubleExponential,
     lccdf_works_on_scalar_arguments) {
  using stan::math::skew_double_exponential_lccdf;

  for (double ys : {0.2, 0.9, 1.1, 3.2}) {
    for (double mus : {0.1, 1.3, 3.0}) {
      for (double sigmas : {0.1, 1.1, 3.2}) {
        for (double taus : {0.01, 0.1, 0.5, 0.9, 0.99}) {
          EXPECT_NEAR(skew_de_ccdf_test(ys, mus, sigmas, taus),
                      skew_double_exponential_lccdf(ys, mus, sigmas, taus),
                      0.001);
        }
      }
    }
  }
}

TEST(ProbDistributionsSkewedDoubleExponential,
     lccdf_works_on_vector_arguments) {
  using stan::math::skew_double_exponential_lccdf;

  std::vector<double> ys{0.2, 0.9, 1.1, 3.2};

  for (double mus : {0.1, 1.3, 3.0}) {
    for (double sigmas : {0.1, 1.1, 3.2}) {
      for (double taus : {0.01, 0.1, 0.5, 0.9, 0.99}) {
        double x = 0.0;
        for (double y : ys)
          x += skew_de_ccdf_test(y, mus, sigmas, taus);
        EXPECT_NEAR(x, skew_double_exponential_lccdf(ys, mus, sigmas, taus),
                    0.001);
      }
    }
  }
}

TEST(ProbDistributionsSkewedDoubleExponential,
     lccdf_works_on_vectorial_y_and_mu) {
  using stan::math::skew_double_exponential_lccdf;
  std::vector<double> ys{0.2, 0.9, 1.1};
  std::vector<double> mus{0.1, 1.3, 3.0};

  for (double sigmas : {0.1, 1.1, 3.2}) {
    for (double taus : {0.01, 0.1, 0.5, 0.9, 0.99}) {
      double x = 0.0;
      for (int i = 0; i < 3; i++)
        x += skew_de_ccdf_test(ys[i], mus[i], sigmas, taus);

      EXPECT_NEAR(x, skew_double_exponential_lccdf(ys, mus, sigmas, taus),
                  0.001);
    }
  }
}

TEST(ProbDistributionsSkewedDoubleExponential,
     lccdf_works_on_vectorial_y_and_sigma) {
  using stan::math::skew_double_exponential_lccdf;
  std::vector<double> ys{0.2, 0.9, 1.1};
  std::vector<double> sigmas{0.1, 1.1, 3.2};

  for (double mus : {0.1, 1.3, 3.0}) {
    for (double taus : {0.01, 0.1, 0.5, 0.9, 0.99}) {
      double x = 0.0;
      for (int i = 0; i < 3; i++)
        x += skew_de_ccdf_test(ys[i], mus, sigmas[i], taus);

      EXPECT_NEAR(x, skew_double_exponential_lccdf(ys, mus, sigmas, taus),
                  0.001);
    }
  }
}

TEST(ProbDistributionsSkewedDoubleExponential,
     lccdf_works_on_vectorial_y_and_tau) {
  using stan::math::skew_double_exponential_lccdf;
  std::vector<double> ys{0.2, 0.9, 1.1};
  std::vector<double> taus{0.1, 0.5, 0.9};

  for (double mus : {0.1, 1.3, 3.0}) {
    for (double sigmas : {0.1, 1.1, 3.2}) {
      double x = 0.0;
      for (int i = 0; i < 3; i++)
        x += skew_de_ccdf_test(ys[i], mus, sigmas, taus[i]);

      EXPECT_NEAR(x, skew_double_exponential_lccdf(ys, mus, sigmas, taus),
                  0.001);
    }
  }
}

TEST(ProbDistributionsSkewedDoubleExponential,
     lccdf_works_on_vectorial_mu_sigma_and_tau) {
  using stan::math::skew_double_exponential_lccdf;

  std::vector<double> mus{0.1, 1.3, 3.0};
  std::vector<double> sigmas{0.1, 1.1, 3.2};
  std::vector<double> taus{0.1, 0.5, 0.9};

  for (double ys : {0.1, 1.3, 3.0}) {
    double x = 0.0;
    for (int i = 0; i < 3; i++)
      x += skew_de_ccdf_test(ys, mus[i], sigmas[i], taus[i]);
    EXPECT_NEAR(x, skew_double_exponential_lccdf(ys, mus, sigmas, taus), 0.001);
  }
}

TEST(ProbDistributionsSkewedDoubleExponential, lccdf_check_errors) {
  using stan::math::skew_double_exponential_lccdf;
  static double inff = std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::skew_double_exponential_lccdf(1.0, 0.0, -1, 0.5),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lccdf(1.0, 0.0, 0.1, -0.5),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lccdf(inff, 0.0, 0.1, 1.5),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lccdf(1.0, inff, 0.1, 1.5),
               std::domain_error);
}

TEST(ProbDistributionsSkewedDoubleExponential, lccdf_check_inconsistent_size) {
  using stan::math::skew_double_exponential_lccdf;

  std::vector<double> mus{0.1, 1.3, 3.0};
  std::vector<double> sigmas{0.1, 1.1, 3.2, 1.0};
  std::vector<double> taus{0.1, 0.5, 0.9};
  EXPECT_THROW(
      stan::math::skew_double_exponential_lccdf(1.0, mus, sigmas, taus),
      std::invalid_argument);
}

TEST(ProbDistributionsSkewedDoubleExponential, cdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = 2;
  double sigma = 2.3;
  double tau = 0.1;

  EXPECT_FLOAT_EQ(
      (stan::math::skew_double_exponential_lccdf(y, mu, sigma, tau)),
      (stan::math::skew_double_exponential_ccdf_log(y, mu, sigma, tau)));
  EXPECT_FLOAT_EQ(
      (stan::math::skew_double_exponential_lccdf<double, double, double,
                                                 double>(y, mu, sigma, tau)),
      (stan::math::skew_double_exponential_ccdf_log<double, double, double,
                                                    double>(y, mu, sigma,
                                                            tau)));
}
