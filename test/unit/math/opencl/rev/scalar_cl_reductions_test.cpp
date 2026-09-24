#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <type_traits>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using stan::math::matrix_cl;
using stan::math::var;
using stan::math::var_value;
using stan::math::opencl::ScalarCl;
using stan::math::opencl::to_host;

TEST(ScalarClReductions, prim_return_types) {
  matrix_cl<double> m_cl(MatrixXd::Random(4, 3));
  matrix_cl<double> v_cl(VectorXd::Random(5));
  matrix_cl<int> i_cl(Eigen::MatrixXi::Ones(4, 3));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::sum(m_cl)), ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::sum(m_cl * 2.0)),
                            ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::sum(i_cl)), int>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::dot_product(v_cl, v_cl)),
                            ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::dot_self(v_cl)),
                            ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::log_sum_exp(m_cl)),
                            ScalarCl<double>>::value));
  EXPECT_TRUE((
      std::is_same<decltype(stan::math::prod(m_cl)), ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::trace(m_cl)),
                            ScalarCl<double>>::value));
  EXPECT_TRUE((
      std::is_same<decltype(stan::math::mean(m_cl)), ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::variance(m_cl)),
                            ScalarCl<double>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::sd(m_cl)), ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::squared_distance(v_cl, v_cl)),
                            ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::distance(v_cl, v_cl)),
                            ScalarCl<double>>::value));
}

TEST(ScalarClReductions, rev_return_types) {
  var_value<matrix_cl<double>> m{matrix_cl<double>(MatrixXd::Random(4, 3))};
  var_value<matrix_cl<double>> v{matrix_cl<double>(VectorXd::Random(5))};
  matrix_cl<double> v_cl(VectorXd::Random(5));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::sum(m)), ScalarCl<var>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::dot_product(v, v_cl)),
                            ScalarCl<var>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::dot_self(v)), ScalarCl<var>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::log_sum_exp(m)),
                            ScalarCl<var>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::prod(m)), ScalarCl<var>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::trace(m)), ScalarCl<var>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::mean(m)), ScalarCl<var>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::variance(m)), ScalarCl<var>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::sd(m)), ScalarCl<var>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::squared_distance(v, v_cl)),
                            ScalarCl<var>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::distance(v, v_cl)),
                            ScalarCl<var>>::value));
  stan::math::recover_memory();
}

TEST(ScalarClReductions, sum_sizes) {
  for (auto dims : std::vector<std::pair<int, int>>{
           {0, 0}, {1, 1}, {5, 1}, {1, 5}, {3, 2000}, {2000, 3}, {700, 700}}) {
    MatrixXd m = MatrixXd::Random(dims.first, dims.second);
    matrix_cl<double> m_cl(m);
    EXPECT_NEAR(to_host(stan::math::sum(m_cl)), m.sum(), 1e-9)
        << dims.first << "x" << dims.second;
    EXPECT_NEAR(to_host(stan::math::sum(m_cl * 2.0 - 1.0)),
                (m.array() * 2.0 - 1.0).sum(), 1e-9)
        << dims.first << "x" << dims.second;
  }
}

TEST(ScalarClReductions, dot_product_checks_sizes) {
  matrix_cl<double> a_cl(VectorXd::Random(3));
  matrix_cl<double> b_cl(VectorXd::Random(4));
  EXPECT_THROW(stan::math::dot_product(a_cl, b_cl), std::invalid_argument);
}

TEST(ScalarClReductions, log_sum_exp_non_finite) {
  double inf = std::numeric_limits<double>::infinity();
  VectorXd a(3);
  a << 1.0, inf, 2.0;
  matrix_cl<double> a_cl(a);
  EXPECT_EQ(to_host(stan::math::log_sum_exp(a_cl)), inf);
  VectorXd b = VectorXd::Constant(4, -inf);
  matrix_cl<double> b_cl(b);
  EXPECT_EQ(to_host(stan::math::log_sum_exp(b_cl)), -inf);
  matrix_cl<double> empty_cl(0, 1);
  EXPECT_EQ(to_host(stan::math::log_sum_exp(empty_cl)), -inf);
  VectorXd c = VectorXd::Random(1000) * 50;
  matrix_cl<double> c_cl(c);
  EXPECT_NEAR(to_host(stan::math::log_sum_exp(c_cl)),
              stan::math::log_sum_exp(c), 1e-10);
  EXPECT_NEAR(to_host(stan::math::log_sum_exp(c_cl * 2.0)),
              stan::math::log_sum_exp(VectorXd(c * 2.0)), 1e-10);
}

TEST_F(AgradRev, ScalarClReductions_resident_chain_gradient) {
  // sum(x) stays on the device and feeds a matrix expression through
  // log_softmax and softmax, whose reverse passes use device reductions.
  VectorXd x_val = VectorXd::Random(50);
  stan::math::vector_v x_cpu = x_val;
  var cpu_res = stan::math::sum(stan::math::log_softmax(x_cpu))
                + stan::math::dot_self(stan::math::softmax(x_cpu))
                + stan::math::sd(x_cpu);
  cpu_res.grad();
  VectorXd cpu_grad = x_cpu.adj();
  stan::math::recover_memory();

  var_value<matrix_cl<double>> x{matrix_cl<double>(x_val)};
  ScalarCl<var> total = stan::math::sum(stan::math::log_softmax(x))
                        + stan::math::dot_self(stan::math::softmax(x))
                        + stan::math::sd(x);
  var res = to_host(total);
  EXPECT_NEAR(res.val(), cpu_res.val(), 1e-10);
  res.grad();
  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(x.adj()), cpu_grad, 1e-10);
}
#endif
