#ifdef STAN_OPENCL
#include <test/unit/math/opencl/util.hpp>
#include <stan/math/opencl/prim.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::matrix_cl;
using stan::math::var;
using std::vector;

TEST(ProbDistributionsOrderedLogisitc, error_checking) {
  int N = 3;
  int C = 5;

  vector<int> y{1, 3, 2};
  vector<int> y_size{1, 3, 1, 2};
  vector<int> y_value1{0, 1, 2};
  vector<int> y_value2{1, 2, 6};
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 2, -0.3;
  Eigen::VectorXd lambda_size(N + 1);
  lambda_size << 0.3, 2, 0.4, 0.3;
  Eigen::VectorXd lambda_value(N);
  lambda_value << 0.3, 2, INFINITY;
  Eigen::VectorXd cuts1(C - 1);
  cuts1 << -0.3, 0.8, 1.8, 3.0;
  Eigen::VectorXd cuts2(C - 1);
  cuts2 << -0.6, 1.8, 2.4, 3.2;
  std::vector<Eigen::VectorXd> cuts{cuts1, cuts2, cuts1};
  Eigen::VectorXd cuts_value1(C - 1);
  cuts_value1 << -0.3, -0.8, 3.0, INFINITY;
  Eigen::VectorXd cuts_value2(C - 1);
  cuts_value2 << 0.3, -0.8, 0, 4.5;

  matrix_cl<int> y_cl(y);
  matrix_cl<int> y_size_cl(y_size);
  matrix_cl<int> y_value1_cl(y_value1);
  matrix_cl<int> y_value2_cl(y_value2);
  matrix_cl<double> lambda_cl(lambda);
  matrix_cl<double> lambda_size_cl(lambda_size);
  matrix_cl<double> lambda_value_cl(lambda_value);
  matrix_cl<double> cuts1_cl(cuts1);
  matrix_cl<double> cuts_cl(cuts);
  matrix_cl<double> cuts_value1_cl(cuts_value1);
  matrix_cl<double> cuts_value2_cl(cuts_value2);

  EXPECT_NO_THROW(stan::math::ordered_logistic_lpmf(y_cl, lambda_cl, cuts_cl));
  EXPECT_NO_THROW(stan::math::ordered_logistic_lpmf(y_cl, lambda_cl, cuts1_cl));

  EXPECT_THROW(stan::math::ordered_logistic_lpmf(y_size_cl, lambda_cl, cuts_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::ordered_logistic_lpmf(y_cl, lambda_size_cl, cuts_cl),
               std::invalid_argument);

  EXPECT_THROW(
      stan::math::ordered_logistic_lpmf(y_value1_cl, lambda_cl, cuts_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_lpmf(y_value2_cl, lambda_cl, cuts_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_lpmf(y_cl, lambda_value_cl, cuts_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_lpmf(y_cl, lambda_cl, cuts_value1_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_lpmf(y_cl, lambda_cl, cuts_value2_cl),
      std::domain_error);
}

auto ordered_logistic_lpmf_functor
    = [](const auto& y, const auto& lambda, const auto& cuts) {
        return stan::math::ordered_logistic_lpmf(y, lambda, cuts);
      };
auto ordered_logistic_lpmf_functor_propto
    = [](const auto& y, const auto& lambda, const auto& cuts) {
        return stan::math::ordered_logistic_lpmf<true>(y, lambda, cuts);
      };

TEST(ProbDistributionsOrderedLogisitc, opencl_matches_cpu_small_simple) {
  int N = 3;
  int C = 5;

  vector<int> y{2, 1, 5};
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 2, 1.1;
  Eigen::VectorXd cuts1(C - 1);
  cuts1 << -0.4, 0.1, 0.3, 4.5;
  Eigen::VectorXd cuts2(C - 1);
  cuts2 << -0.8, 0.2, 0.3, 4.9;
  std::vector<Eigen::VectorXd> cuts{cuts1, cuts2, cuts1};

  stan::math::test::compare_cpu_opencl_prim_rev(ordered_logistic_lpmf_functor,
                                                y, lambda, cuts1);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_lpmf_functor_propto, y, lambda, cuts1);
  stan::math::test::compare_cpu_opencl_prim_rev(ordered_logistic_lpmf_functor,
                                                y, lambda, cuts);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_lpmf_functor_propto, y, lambda, cuts);
}

TEST(ProbDistributionsOrderedLogisitc, opencl_matches_cpu_zero_instances) {
  int N = 0;
  int C = 5;

  vector<int> y{};
  Eigen::VectorXd lambda(N);
  Eigen::VectorXd cuts1(C - 1);
  cuts1 << -0.4, 0.1, 0.3, 4.5;
  std::vector<Eigen::VectorXd> cuts{};

  stan::math::test::compare_cpu_opencl_prim_rev(ordered_logistic_lpmf_functor,
                                                y, lambda, cuts1);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_lpmf_functor_propto, y, lambda, cuts1);
  stan::math::test::compare_cpu_opencl_prim_rev(ordered_logistic_lpmf_functor,
                                                y, lambda, cuts);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_lpmf_functor_propto, y, lambda, cuts);
}

TEST(ProbDistributionsOrderedLogisitc, opencl_matches_cpu_single_class) {
  int N = 3;
  int C = 1;

  vector<int> y{1, 1, 1};
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 2, 0.4;
  Eigen::VectorXd cuts1(C - 1);
  std::vector<Eigen::VectorXd> cuts{cuts1, cuts1, cuts1};

  stan::math::test::compare_cpu_opencl_prim_rev(ordered_logistic_lpmf_functor,
                                                y, lambda, cuts1);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_lpmf_functor_propto, y, lambda, cuts1);
  stan::math::test::compare_cpu_opencl_prim_rev(ordered_logistic_lpmf_functor,
                                                y, lambda, cuts);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_lpmf_functor_propto, y, lambda, cuts);
}

TEST(ProbDistributionsOrderedLogisitc, opencl_matches_cpu_big) {
  int N = 153;
  int C = 43;

  vector<int> y(N);
  for (int i = 0; i < N; i++) {
    y[i] = Array<int, Dynamic, 1>::Random(1, 1).abs()(0) % C + 1;
  }
  Eigen::VectorXd lambda = Eigen::VectorXd::Random(N, 1);

  Eigen::VectorXd cuts1;
  std::vector<Eigen::VectorXd> cuts;
  for (int i = 0; i < N; i++) {
    cuts1 = Array<double, Dynamic, 1>::Random(C - 1, 1).abs() + 1e-5;
    for (int i = 1; i < C - 1; i++) {
      cuts1[i] += cuts1[i - 1];
    }
    cuts.push_back(cuts1);
  }

  stan::math::test::compare_cpu_opencl_prim_rev(ordered_logistic_lpmf_functor,
                                                y, lambda, cuts1);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_lpmf_functor_propto, y, lambda, cuts1);
  stan::math::test::compare_cpu_opencl_prim_rev(ordered_logistic_lpmf_functor,
                                                y, lambda, cuts);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_lpmf_functor_propto, y, lambda, cuts);
}

#endif
