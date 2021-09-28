#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsChiSquare, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.0, 1.2;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, -0.8, 0.5;
  Eigen::VectorXd nu(N);
  nu << 0.3, 0.8, 1.2;
  Eigen::VectorXd nu_size(N - 1);
  nu_size << 0.3, 0.8;
  Eigen::VectorXd nu_value(N);
  nu_value << 0.3, 0, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> nu_cl(nu);
  stan::math::matrix_cl<double> nu_size_cl(nu_size);
  stan::math::matrix_cl<double> nu_value_cl(nu_value);

  EXPECT_NO_THROW(stan::math::chi_square_lpdf(y_cl, nu_cl));

  EXPECT_THROW(stan::math::chi_square_lpdf(y_size_cl, nu_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::chi_square_lpdf(y_cl, nu_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::chi_square_lpdf(y_value_cl, nu_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::chi_square_lpdf(y_cl, nu_value_cl),
               std::domain_error);
}

auto chi_square_lpdf_functor = [](const auto& y, const auto& nu) {
  return stan::math::chi_square_lpdf(y, nu);
};
auto chi_square_lpdf_functor_propto = [](const auto& y, const auto& nu) {
  return stan::math::chi_square_lpdf<true>(y, nu);
};

TEST(ProbDistributionsChiSquare, opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.0, 1.2;
  Eigen::VectorXd nu(N);
  nu << 0.2, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(chi_square_lpdf_functor, y, nu);
  stan::math::test::compare_cpu_opencl_prim_rev(chi_square_lpdf_functor_propto,
                                                y, nu);
  stan::math::test::compare_cpu_opencl_prim_rev(chi_square_lpdf_functor,
                                                y.transpose().eval(), nu);
  stan::math::test::compare_cpu_opencl_prim_rev(chi_square_lpdf_functor_propto,
                                                y, nu.transpose().eval());
}

TEST(ProbDistributionsChiSquare, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 1.1;
  Eigen::VectorXd nu(N);
  nu << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      chi_square_lpdf_functor, y_scal, nu);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      chi_square_lpdf_functor_propto, y_scal, nu);
}

TEST(ProbDistributionsChiSquare, opencl_broadcast_nu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.0, 1.2;
  double nu_scal = 0.4;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      chi_square_lpdf_functor, y, nu_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      chi_square_lpdf_functor_propto, y, nu_scal);
}

TEST(ProbDistributionsChiSquare, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> nu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(chi_square_lpdf_functor, y, nu);
  stan::math::test::compare_cpu_opencl_prim_rev(chi_square_lpdf_functor_propto,
                                                y, nu);
  stan::math::test::compare_cpu_opencl_prim_rev(
      chi_square_lpdf_functor, y.transpose().eval(), nu.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(chi_square_lpdf_functor_propto,
                                                y.transpose().eval(),
                                                nu.transpose().eval());
}

#endif
