#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::matrix_cl;
using stan::math::var;
using stan::test::expect_near_rel;
using std::vector;

TEST(ProbDistributionsNegBinomial2LogGLM, error_checking) {
  int N = 3;
  int M = 2;

  vector<int> y{0, 1, 5};
  vector<int> y_size{0, 1, 5, 0};
  vector<int> y_value{1, 4, -23};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, Dynamic> x_size1(N - 1, M);
  x_size1 << -12, 46, -42, 24;
  Matrix<double, Dynamic, Dynamic> x_size2(N, M - 1);
  x_size2 << -12, 46, -42;
  Matrix<double, Dynamic, Dynamic> x_value(N, M);
  x_value << -12, 46, -42, 24, 25, -INFINITY;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> beta_size(M + 1, 1);
  beta_size << 0.3, 2, 0.4;
  Matrix<double, Dynamic, 1> beta_value(M, 1);
  beta_value << 0.3, INFINITY;
  Matrix<double, Dynamic, 1> alpha(N, 1);
  alpha << 0.3, -0.8, 1.8;
  Matrix<double, Dynamic, 1> alpha_size(N - 1, 1);
  alpha_size << 0.3, -0.8;
  Matrix<double, Dynamic, 1> alpha_value(N, 1);
  alpha_value << 0.3, -0.8, NAN;
  double phi1 = 1.2;
  Matrix<double, Dynamic, 1> phi2(N, 1);
  phi2 << 0.1, 0.1, 3.2;
  Matrix<double, Dynamic, 1> phi_size(N - 1, 1);
  phi_size << 0.3, 0.8;
  Matrix<double, Dynamic, 1> phi_value1(N, 1);
  phi_value1 << 0.3, 0.8, NAN;
  Matrix<double, Dynamic, 1> phi_value2(N, 1);
  phi_value2 << 0.3, -0.8, 3;

  matrix_cl<double> x_cl(x);
  matrix_cl<double> x_size1_cl(x_size1);
  matrix_cl<double> x_size2_cl(x_size2);
  matrix_cl<double> x_value_cl(x_value);
  matrix_cl<int> y_cl(y);
  matrix_cl<int> y_size_cl(y_size);
  matrix_cl<int> y_value_cl(y_value);
  matrix_cl<double> beta_cl(beta);
  matrix_cl<double> beta_size_cl(beta_size);
  matrix_cl<double> beta_value_cl(beta_value);
  matrix_cl<double> alpha_cl(alpha);
  matrix_cl<double> alpha_size_cl(alpha_size);
  matrix_cl<double> alpha_value_cl(alpha_value);
  matrix_cl<double> phi2_cl(phi2);
  matrix_cl<double> phi_size_cl(phi_size);
  matrix_cl<double> phi_value1_cl(phi_value1);
  matrix_cl<double> phi_value2_cl(phi_value2);

  EXPECT_NO_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                          beta_cl, phi1));
  EXPECT_NO_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                          beta_cl, phi2_cl));

  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_size_cl, x_cl,
                                                       alpha_cl, beta_cl, phi1),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_size1_cl,
                                                       alpha_cl, beta_cl, phi1),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_size2_cl,
                                                       alpha_cl, beta_cl, phi1),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(
                   y_cl, x_cl, alpha_size_cl, beta_cl, phi1),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                       beta_size_cl, phi1),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                       beta_size_cl, phi1),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                       beta_cl, phi_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_value_cl, x_cl,
                                                       alpha_cl, beta_cl, phi1),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_value_cl,
                                                       alpha_cl, beta_cl, phi1),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(
                   y_cl, x_cl, alpha_value_cl, beta_cl, phi1),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                       beta_value_cl, phi1),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                       beta_cl, phi_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                       beta_cl, phi_value2_cl),
               std::domain_error);
}

auto neg_binomial_2_log_glm_lpmf_functor
    = [](const auto& y, const auto& x, const auto& alpha, const auto& beta,
         const auto& phi) {
        return stan::math::neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi);
      };
auto neg_binomial_2_log_glm_lpmf_functor_propto
    = [](const auto& y, const auto& x, const auto& alpha, const auto& beta,
         const auto& phi) {
        return stan::math::neg_binomial_2_log_glm_lpmf<true>(y, x, alpha, beta,
                                                             phi);
      };

TEST(ProbDistributionsNegBinomial2LogGLM, opencl_matches_cpu_small_simple) {
  int N = 3;
  int M = 2;

  vector<int> y{0, 1, 5};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  double alpha = 0.3;
  double phi = 13.2;

  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor, y, x, alpha, beta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor_propto, y, x, alpha, beta, phi);
}

TEST(ProbDistributionsNegBinomial2LogGLM, opencl_broadcast_y) {
  int N = 3;
  int M = 2;

  int y_scal = 1;
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  double alpha = 0.3;
  double phi = 13.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_log_glm_lpmf_functor, y_scal, x, alpha, beta, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_log_glm_lpmf_functor_propto, y_scal, x, alpha, beta, phi);
}

TEST(ProbDistributionsNegBinomial2LogGLM, opencl_matches_cpu_zero_instances) {
  int N = 0;
  int M = 2;

  vector<int> y{};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  double alpha = 0.3;
  double phi = 13.2;

  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor, y, x, alpha, beta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor_propto, y, x, alpha, beta, phi);
}

TEST(ProbDistributionsNegBinomial2LogGLM, opencl_matches_cpu_zero_attributes) {
  int N = 3;
  int M = 0;

  vector<int> y{0, 1, 5};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, 1> beta(M, 1);
  double alpha = 0.3;
  double phi = 13.2;

  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor, y, x, alpha, beta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor_propto, y, x, alpha, beta, phi);
}

TEST(ProbDistributionsNegBinomial2LogGLM,
     opencl_matches_cpu_small_vector_alpha_phi) {
  int N = 3;
  int M = 2;

  vector<int> y{0, 1, 5};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> alpha(N, 1);
  alpha << 0.3, -0.8, 1.8;
  Matrix<double, Dynamic, 1> phi(N, 1);
  phi << 0.1, 0.5, 1.2;

  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor, y, x, alpha, beta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor_propto, y, x, alpha, beta, phi);
}

TEST(ProbDistributionsNegBinomial2LogGLM, opencl_matches_cpu_big) {
  int N = 153;
  int M = 71;

  vector<int> y(N);
  for (int i = 0; i < N; i++) {
    y[i] = Array<int, Dynamic, 1>::Random(1, 1).abs()(0);
  }
  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Random(N, M);
  Matrix<double, Dynamic, 1> beta = Matrix<double, Dynamic, 1>::Random(M, 1);
  Matrix<double, Dynamic, 1> alpha = Matrix<double, Dynamic, 1>::Random(N, 1);
  Matrix<double, Dynamic, 1> phi
      = Array<double, Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor, y, x, alpha, beta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_glm_lpmf_functor_propto, y, x, alpha, beta, phi);
}

TEST(ProbDistributionsNegBinomial2LogGLM, opencl_large_shapes_reference) {
  // The tests above compare OpenCL against the CPU and cannot see an error
  // that both share. These are absolute references from mpmath at 90 digits,
  // checked against 140 digits, as in
  // test/unit/math/rev/prob/large_shapes_test.cpp. The
  // value formed from lgamma(phi), lgamma(y + phi) and phi log(phi), and
  // the phi partial formed from digamma(y + phi) - digamma(phi), keep no
  // correct digits for phi near 1e15. One attribute; the phi partial is
  // compared as phi * d/dphi. The last row is a small phi where the plain
  // differences are correct.
  struct TestValue {
    vector<int> y;
    vector<double> x;
    double alpha;
    double beta;
    double phi;
    double value;
    double d_alpha;
    double d_beta;
    double d_phi;
  };
  const vector<double> x5{0.5, -0.25, 1.25, 0.0, 2.0};
  const vector<TestValue> test_values = {
      {{0},
       {0.5},
       1.0,
       0.75,
       0x1.c6bf526340000p+49,
       -3.9550767229205693,
       -3.9550767229205615,
       -1.9775383614602807,
       -7.8213159420940446e-30},
      {{5},
       {0.5},
       1.0,
       0.75,
       0x1.c6bf526340000p+49,
       -1.8675684657026251,
       1.0449232770794187,
       5.2246163853970937e-1,
       1.9540676725087929e-30},
      {{0},
       {0.5},
       -0x1.999999999999ap-3,
       0.75,
       0x1.37807ed5e8000p+50,
       -1.1912462166123576,
       -1.1912462166123571,
       -5.9562310830617854e-1,
       -3.7803493755481261e-31},
      {{2, 4, 3, 5, 6},
       x5,
       1.0,
       0.75,
       0x1.c6bf526340000p+49,
       -1.3267966555421410e+1,
       -8.0507631204932393,
       -1.8705862362560038e+1,
       -2.2918189242185528e-29},
      {{0, 1, 5, 3, 57},
       x5,
       1.0,
       0.75,
       0x1.c6bf526340000p+49,
       -5.5025862739499810e+1,
       3.7949236879506146e+1,
       8.5544137637438705e+1,
       -9.8183556706826074e-28},
      {{2, 4, 3, 5, 6},
       x5,
       1.0,
       0.75,
       0x1.d1a94a2000000p+39,
       -1.3267966555398514e+1,
       -8.0507631203930684,
       -1.8705862362370543e+1,
       -2.2918189241713975e-23},
      {{0},
       {0.5},
       1.0,
       0.75,
       0x1.9000000000000p+6,
       -3.8788665246722319,
       -3.8046018026251338,
       -1.9023009013125669,
       -7.4264722047098149e-4},
      {{0, 1, 5, 3, 57},
       x5,
       1.0,
       0.75,
       0x1.9000000000000p+6,
       -4.7240731087544224e+1,
       3.3378922648644250e+1,
       7.6036039699167587e+1,
       -6.2120571346832957e-2},
  };
  auto expect_reference = [](const TestValue& t, const var& lp,
                             double alpha_adj, double beta_adj, double phi_adj,
                             const char* signature) {
    const std::string where = std::string(signature)
                              + ": y[0] = " + std::to_string(t.y[0])
                              + ", N = " + std::to_string(t.y.size())
                              + ", phi = " + std::to_string(t.phi);
    EXPECT_NEAR(lp.val(), t.value, 1e-12 * std::max(1.0, std::fabs(t.value)))
        << where;
    EXPECT_NEAR(alpha_adj, t.d_alpha,
                1e-11 * std::max(1.0, std::fabs(t.d_alpha)))
        << where;
    EXPECT_NEAR(beta_adj, t.d_beta, 1e-11 * std::max(1.0, std::fabs(t.d_beta)))
        << where;
    const double gp = t.phi * t.d_phi;
    EXPECT_NEAR(t.phi * phi_adj, gp, 1e-11 * std::max(1.0, std::fabs(gp)))
        << where;
  };
  for (const auto& t : test_values) {
    const int N = t.y.size();
    Matrix<double, Dynamic, Dynamic> x(N, 1);
    for (int i = 0; i < N; ++i) {
      x(i, 0) = t.x[i];
    }
    matrix_cl<int> y_cl(t.y);
    matrix_cl<double> x_cl(x);

    // y vector, alpha and phi scalars: the kernel sums the phi partial
    {
      var alpha = t.alpha;
      Matrix<var, Dynamic, 1> beta(1);
      beta << t.beta;
      var phi = t.phi;
      auto beta_cl = stan::math::to_matrix_cl(beta);
      var lp = stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha,
                                                       beta_cl, phi);
      lp.grad();
      expect_reference(t, lp, alpha.adj(), beta(0).adj(), phi.adj(),
                       "scalar alpha and phi");
      stan::math::recover_memory();
    }

    // alpha and phi vectors: the kernel computes the value term and the phi
    // partial of each instance
    {
      Matrix<var, Dynamic, 1> alpha(N);
      Matrix<var, Dynamic, 1> beta(1);
      beta << t.beta;
      Matrix<var, Dynamic, 1> phi(N);
      for (int i = 0; i < N; ++i) {
        alpha(i) = t.alpha;
        phi(i) = t.phi;
      }
      auto alpha_cl = stan::math::to_matrix_cl(alpha);
      auto beta_cl = stan::math::to_matrix_cl(beta);
      auto phi_cl = stan::math::to_matrix_cl(phi);
      var lp = stan::math::neg_binomial_2_log_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                       beta_cl, phi_cl);
      lp.grad();
      double alpha_adj = 0;
      double phi_adj = 0;
      for (int i = 0; i < N; ++i) {
        alpha_adj += alpha(i).adj();
        phi_adj += phi(i).adj();
      }
      expect_reference(t, lp, alpha_adj, beta(0).adj(), phi_adj,
                       "vector alpha and phi");
      stan::math::recover_memory();
    }

    // y, alpha and phi scalars: the value term is computed on the host
    if (N == 1) {
      var alpha = t.alpha;
      Matrix<var, Dynamic, 1> beta(1);
      beta << t.beta;
      var phi = t.phi;
      auto beta_cl = stan::math::to_matrix_cl(beta);
      var lp = stan::math::neg_binomial_2_log_glm_lpmf(t.y[0], x_cl, alpha,
                                                       beta_cl, phi);
      lp.grad();
      expect_reference(t, lp, alpha.adj(), beta(0).adj(), phi.adj(),
                       "scalar y, alpha and phi");
      stan::math::recover_memory();
    }
  }
}
#endif
