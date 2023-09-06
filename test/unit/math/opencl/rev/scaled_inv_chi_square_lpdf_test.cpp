#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsScaledInvChiSquare, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  Eigen::VectorXd nu(N);
  nu << 0.3, 0.8, 1.0;
  Eigen::VectorXd nu_size(N - 1);
  nu_size << 0.3, 0.8;
  Eigen::VectorXd nu_value1(N);
  nu_value1 << 0.3, 0, 0.5;
  Eigen::VectorXd nu_value2(N);
  nu_value2 << 0.3, INFINITY, 0.5;

  Eigen::VectorXd s(N);
  s << 0.3, 0.8, 1.0;
  Eigen::VectorXd s_size(N - 1);
  s_size << 0.3, 0.8;
  Eigen::VectorXd s_value1(N);
  s_value1 << 0.3, 0, 0.5;
  Eigen::VectorXd s_value2(N);
  s_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> nu_cl(nu);
  stan::math::matrix_cl<double> nu_size_cl(nu_size);
  stan::math::matrix_cl<double> nu_value1_cl(nu_value1);
  stan::math::matrix_cl<double> nu_value2_cl(nu_value2);
  stan::math::matrix_cl<double> s_cl(s);
  stan::math::matrix_cl<double> s_size_cl(s_size);
  stan::math::matrix_cl<double> s_value1_cl(s_value1);
  stan::math::matrix_cl<double> s_value2_cl(s_value2);

  EXPECT_NO_THROW(stan::math::scaled_inv_chi_square_lpdf(y_cl, nu_cl, s_cl));

  EXPECT_THROW(stan::math::scaled_inv_chi_square_lpdf(y_size_cl, nu_cl, s_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::scaled_inv_chi_square_lpdf(y_cl, nu_size_cl, s_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::scaled_inv_chi_square_lpdf(y_cl, nu_cl, s_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::scaled_inv_chi_square_lpdf(y_value_cl, nu_cl, s_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::scaled_inv_chi_square_lpdf(y_cl, nu_value1_cl, s_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::scaled_inv_chi_square_lpdf(y_cl, nu_value2_cl, s_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::scaled_inv_chi_square_lpdf(y_cl, nu_cl, s_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::scaled_inv_chi_square_lpdf(y_cl, nu_cl, s_value2_cl),
               std::domain_error);
}

auto scaled_inv_chi_square_lpdf_functor
    = [](const auto& y, const auto& nu, const auto& s) {
        return stan::math::scaled_inv_chi_square_lpdf(y, nu, s);
      };
auto scaled_inv_chi_square_lpdf_functor_propto
    = [](const auto& y, const auto& nu, const auto& s) {
        return stan::math::scaled_inv_chi_square_lpdf<true>(y, nu, s);
      };

TEST(ProbDistributionsScaledInvChiSquare, opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd nu(N);
  nu << 0.3, 0.8, 1.0;
  Eigen::VectorXd s(N);
  s << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor_propto, y, nu, s);
  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor, y, nu, s);
  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor_propto, y.transpose().eval(),
      nu.transpose().eval(), s.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor, y.transpose().eval(),
      nu.transpose().eval(), s.transpose().eval());
}

TEST(ProbDistributionsScaledInvChiSquare, opencl_matches_cpu_small_y_negative) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, -0.8, 1.0;
  Eigen::VectorXd nu(N);
  nu << 0.3, 0.8, 1.0;
  Eigen::VectorXd s(N);
  s << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor, y, nu, s);
  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor_propto, y, nu, s);
}

TEST(ProbDistributionsScaledInvChiSquare, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd nu(N);
  nu << 0.5, 1.2, 1.0;
  Eigen::VectorXd s(N);
  s << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      scaled_inv_chi_square_lpdf_functor, y_scal, nu, s);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      scaled_inv_chi_square_lpdf_functor_propto, y_scal, nu, s);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      scaled_inv_chi_square_lpdf_functor, y_scal, nu.transpose().eval(), s);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      scaled_inv_chi_square_lpdf_functor_propto, y_scal, nu,
      s.transpose().eval());
}

TEST(ProbDistributionsScaledInvChiSquare, opencl_broadcast_nu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double nu_scal = 12.3;
  Eigen::VectorXd s(N);
  s << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      scaled_inv_chi_square_lpdf_functor, y, nu_scal, s);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      scaled_inv_chi_square_lpdf_functor_propto, y, nu_scal, s);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      scaled_inv_chi_square_lpdf_functor, y.transpose().eval(), nu_scal, s);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      scaled_inv_chi_square_lpdf_functor_propto, y, nu_scal,
      s.transpose().eval());
}

TEST(ProbDistributionsScaledInvChiSquare, opencl_broadcast_s) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd nu(N);
  nu << 0.3, 0.8, 1.0;
  double s_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      scaled_inv_chi_square_lpdf_functor, y, nu, s_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      scaled_inv_chi_square_lpdf_functor_propto, y, nu, s_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      scaled_inv_chi_square_lpdf_functor, y.transpose().eval(), nu, s_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      scaled_inv_chi_square_lpdf_functor_propto, y, nu.transpose().eval(),
      s_scal);
}

TEST(ProbDistributionsScaledInvChiSquare, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> nu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> s
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor, y, nu, s);
  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor_propto, y, nu, s);
  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor, y.transpose().eval(),
      nu.transpose().eval(), s.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      scaled_inv_chi_square_lpdf_functor_propto, y.transpose().eval(),
      nu.transpose().eval(), s.transpose().eval());
}

#endif
