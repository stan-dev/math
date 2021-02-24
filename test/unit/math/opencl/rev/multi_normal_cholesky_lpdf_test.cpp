#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsMultiNormalCholesky, error_checking) {
  int N = 3;
  int M = 2;

  Eigen::MatrixXd y1(N, M);
  y1 << -0.1, 0.3, 1.4, 0.7, 5.5, 0;
  Eigen::VectorXd y2(N);
  y2 << 0.1, 0.4, 3.5;
  Eigen::MatrixXd y_size1(N - 1, M);
  y_size1 << 0.4, 0.5, 0.6, 0.5;
  Eigen::MatrixXd y_size2(N, M + 1);
  y_size2 << 0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4;
  Eigen::VectorXd y_value(N);
  y_value << NAN, 0.5, 0.6;

  Eigen::MatrixXd mu1(N, M);
  mu1 << 0.3, 0.8, -1.0, 0.12, 0.5, -12.3;
  Eigen::VectorXd mu2(N);
  mu2 << 0.3, -0.8, 1.0;
  Eigen::MatrixXd mu_size1(N - 1, M);
  mu_size1 << 0.3, 0.8, 0.3, 0.8;
  Eigen::MatrixXd mu_size2(N, M + 1);
  mu_size2 << 0.3, 0.8, 1.0, 0.12, 0.5, 12.3, 1, 2, 3;
  Eigen::VectorXd mu_value(N);
  mu_value << 0.3, INFINITY, 0.5;

  Eigen::MatrixXd L(N, N);
  L << 1, 0, 0, 2, 3, 0, -INFINITY, 5, -6;
  Eigen::MatrixXd L_size1(N, N - 1);
  L_size1 << 1, 0, 0, 2, 3, 0;
  Eigen::MatrixXd L_size2(N - 1, N);
  L_size2 << 1, 0, 0, 2, 3, 0;

  stan::math::matrix_cl<double> y1_cl(y1);
  stan::math::matrix_cl<double> y2_cl(y2);
  stan::math::matrix_cl<double> y_size1_cl(y_size1);
  stan::math::matrix_cl<double> y_size2_cl(y_size2);
  stan::math::matrix_cl<double> y_value_cl(y_value);

  stan::math::matrix_cl<double> mu1_cl(mu1);
  stan::math::matrix_cl<double> mu2_cl(mu2);
  stan::math::matrix_cl<double> mu_size1_cl(mu_size1);
  stan::math::matrix_cl<double> mu_size2_cl(mu_size2);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);

  stan::math::matrix_cl<double> L_cl(L);
  stan::math::matrix_cl<double> L_size1_cl(L_size1);
  stan::math::matrix_cl<double> L_size2_cl(L_size2);

  EXPECT_NO_THROW(stan::math::multi_normal_cholesky_lpdf(y1_cl, mu1_cl, L_cl));
  EXPECT_NO_THROW(stan::math::multi_normal_cholesky_lpdf(y1_cl, mu2_cl, L_cl));
  EXPECT_NO_THROW(stan::math::multi_normal_cholesky_lpdf(y2_cl, mu1_cl, L_cl));

  EXPECT_THROW(stan::math::multi_normal_cholesky_lpdf(y_size1_cl, mu1_cl, L_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_cholesky_lpdf(y_size2_cl, mu1_cl, L_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_cholesky_lpdf(y1_cl, mu_size1_cl, L_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_cholesky_lpdf(y1_cl, mu_size2_cl, L_cl),
               std::invalid_argument);
  EXPECT_THROW(
      stan::math::multi_normal_cholesky_lpdf(y1_cl, mu1_cl, L_size1_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::multi_normal_cholesky_lpdf(y1_cl, mu1_cl, L_size2_cl),
      std::invalid_argument);

  EXPECT_THROW(stan::math::multi_normal_cholesky_lpdf(y_value_cl, mu1_cl, L_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::multi_normal_cholesky_lpdf(y1_cl, mu_value_cl, L_cl),
               std::domain_error);
}

auto multi_normal_cholesky_lpdf_functor
    = [](const auto& y, const auto& mu, const auto& L) {
        return stan::math::multi_normal_cholesky_lpdf(y, mu, L);
      };
auto multi_normal_cholesky_lpdf_functor_propto
    = [](const auto& y, const auto& mu, const auto& L) {
        return stan::math::multi_normal_cholesky_lpdf<true>(y, mu, L);
      };

TEST(ProbDistributionsMultiNormalCholesky, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y1(N);
  y1 << 0.5, 0.4, 0.1;
  Eigen::RowVectorXd y2(N);
  y2 << 0.6, 0.2, 0.2;
  std::vector<Eigen::VectorXd> y3{y1, y2};
  std::vector<Eigen::RowVectorXd> y4{y2, y1};
  Eigen::VectorXd mu1(N);
  mu1 << 0.5, 0.1, 12.3;
  Eigen::RowVectorXd mu2(N);
  mu2 << 2.1, 3.4, 2.3;
  std::vector<Eigen::VectorXd> mu3{mu1, mu2};
  std::vector<Eigen::RowVectorXd> mu4{mu2, mu1};
  Eigen::MatrixXd L(N, N);
  L << 1, 0, 0, 2, 3, 0, -4, 5, -6;

  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y1, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y1, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y1, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y1, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y1, mu4, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y1, mu4, L);

  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y3, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y3, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y3, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y3, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y3, mu4, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y3, mu4, L);

  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y4, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y4, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y4, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y4, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y4, mu4, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y4, mu4, L);
}

TEST(ProbDistributionsMultiNormalCholesky, opencl_matches_cpu_big) {
  int N = 73;
  int M = 11;
  Eigen::VectorXd y1;
  Eigen::VectorXd mu1;
  std::vector<Eigen::VectorXd> y3;
  std::vector<Eigen::VectorXd> mu3;
  std::vector<Eigen::RowVectorXd> y4;
  std::vector<Eigen::RowVectorXd> mu4;
  Eigen::MatrixXd L
      = Eigen::MatrixXd::Random(N, N).triangularView<Eigen::Lower>();

  for (int i = 0; i < M; i++) {
    y1 = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
    mu1 = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
    y3.push_back(y1);
    mu3.push_back(mu1);
    y4.push_back(y1);
    mu4.push_back(mu1);
  }
  Eigen::RowVectorXd y2 = y1;
  Eigen::RowVectorXd mu2 = mu1;

  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y1, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y1, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y1, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y1, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y1, mu4, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y1, mu4, L);

  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y3, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y3, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y3, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y3, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y3, mu4, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y3, mu4, L);

  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y4, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y4, mu1, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y4, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y4, mu3, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor, y4, mu4, L);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multi_normal_cholesky_lpdf_functor_propto, y4, mu4, L);
}

#endif
