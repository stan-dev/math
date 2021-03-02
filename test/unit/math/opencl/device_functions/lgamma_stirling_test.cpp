#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/kernels/device_functions/lgamma_stirling.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/expect_near_rel.hpp>
#include <string>

static const std::string test__lgamma_stirling_kernel_code
    = STRINGIFY(__kernel void test(__global double *B, __global double *A) {
        const int i = get_global_id(0);
        B[i] = lgamma_stirling(A[i]);
      });

const stan::math::opencl_kernels::kernel_cl<
    stan::math::opencl_kernels::out_buffer,
    stan::math::opencl_kernels::in_buffer>
    lgamma_stirling(
        "test", {stan::math::opencl_kernels::lgamma_stirling_device_function,
                 test__lgamma_stirling_kernel_code});

TEST(MathMatrixCL, lgamma_stirling) {
  Eigen::VectorXd a = Eigen::VectorXd::Random(1000).array() * 30;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl(1000, 1);
  lgamma_stirling(cl::NDRange(1000), res_cl, a_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  EXPECT_NEAR_REL(res, a.unaryExpr([](double x) {
    return stan::math::lgamma_stirling(x);
  }));
}

TEST(MathMatrixCL, lgamma_stirling_edge_cases) {
  Eigen::VectorXd a(3);
  a << NAN, -1, 1.0E50;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl(3, 1);
  lgamma_stirling(cl::NDRange(3), res_cl, a_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  EXPECT_NEAR_REL(res, a.unaryExpr([](double x) {
    return stan::math::lgamma_stirling(x);
  }));
}

#endif
