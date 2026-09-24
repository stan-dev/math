#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using Eigen::VectorXd;
using stan::math::matrix_cl;
using stan::math::var;
using stan::math::var_value;
using stan::math::vector_v;
using stan::math::opencl::ScalarCl;
using stan::math::opencl::to_host;

namespace {
/**
 * Compares a constraint with a device var log density accumulator against
 * the CPU version with a CPU var accumulator.
 * @param cpu_f calls the CPU constraint with (x, lp)
 * @param cl_f calls the OpenCL constraint with (x, lp)
 */
template <typename F_cpu, typename F_cl>
void expect_scalar_cl_lp(const F_cpu& cpu_f, const F_cl& cl_f) {
  VectorXd x_val(5);
  x_val << -1.5, -0.3, 0.0, 0.7, 2.1;

  vector_v x_cpu = x_val;
  var lp_cpu = 0.5;
  vector_v res_cpu = cpu_f(x_cpu, lp_cpu);
  var f_cpu = stan::math::sum(res_cpu) + 2.0 * lp_cpu;
  f_cpu.grad();
  VectorXd x_grad = x_cpu.adj();
  double f_val = f_cpu.val();
  stan::math::recover_memory();

  var_value<matrix_cl<double>> x{matrix_cl<double>(x_val)};
  ScalarCl<var> lp(0.5);
  var_value<matrix_cl<double>> res = cl_f(x, lp);
  ScalarCl<var> f = stan::math::sum(res) + 2.0 * lp;
  var f_host = to_host(f);
  EXPECT_NEAR(f_host.val(), f_val, 1e-10);
  f_host.grad();
  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(x.adj()), x_grad, 1e-10);
  stan::math::recover_memory();
}
}  // namespace

TEST_F(AgradRev, ScalarClConstraints_device_lp) {
  expect_scalar_cl_lp(
      [](const auto& x, auto& lp) {
        return stan::math::lb_constrain(x, 1.5, lp);
      },
      [](const auto& x, auto& lp) {
        return stan::math::lb_constrain(x, ScalarCl<double>(1.5), lp);
      });
  expect_scalar_cl_lp(
      [](const auto& x, auto& lp) {
        return stan::math::ub_constrain(x, -0.5, lp);
      },
      [](const auto& x, auto& lp) {
        return stan::math::ub_constrain(x, -0.5, lp);
      });
  expect_scalar_cl_lp(
      [](const auto& x, auto& lp) {
        return stan::math::lub_constrain(x, -1.0, 3.0, lp);
      },
      [](const auto& x, auto& lp) {
        return stan::math::lub_constrain(x, ScalarCl<double>(-1.0), 3.0, lp);
      });
  expect_scalar_cl_lp(
      [](const auto& x, auto& lp) {
        return stan::math::offset_multiplier_constrain(x, 0.3, 1.7, lp);
      },
      [](const auto& x, auto& lp) {
        return stan::math::offset_multiplier_constrain(
            x, 0.3, ScalarCl<double>(1.7), lp);
      });
}

TEST_F(AgradRev, ScalarClConstraints_prim_device_lp) {
  VectorXd x_val(4);
  x_val << -1.0, 0.0, 0.5, 2.0;
  matrix_cl<double> x(x_val);
  ScalarCl<double> lp(0.25);
  double lp_host = 0.25;
  matrix_cl<double> res
      = stan::math::lub_constrain(x, 0.0, ScalarCl<double>(2.0), lp);
  VectorXd res_cpu = stan::math::lub_constrain(x_val, 0.0, 2.0, lp_host);
  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(res), res_cpu, 1e-12);
  EXPECT_NEAR(to_host(lp), lp_host, 1e-12);
}
#endif
