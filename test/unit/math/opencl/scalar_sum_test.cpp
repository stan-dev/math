#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <stan/math/opencl/kernels/scalar_sum.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <string>
#include <vector>

using Eigen::MatrixXd;
using stan::math::matrix_cl;
using stan::math::opencl::ScalarCl;
using stan::math::opencl::to_host;
using stan::math::opencl::internal::sum_into;

namespace {
/**
 * Builds a program from the given sources with `-cl-std=CL1.2` and the same
 * macro definitions the library uses.
 */
void build_cl12(const std::vector<std::string>& sources) {
  auto& ctx = stan::math::opencl_context;
  std::string opts = "-cl-std=CL1.2";
  for (auto&& opt : ctx.base_opts()) {
    opts += " -D" + opt.first + "=" + std::to_string(opt.second);
  }
  cl::Program program(ctx.context(), sources);
  try {
    program.build({ctx.device()}, opts.c_str());
  } catch (const cl::Error& e) {
    FAIL() << "OpenCL 1.2 build failed: "
           << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(ctx.device()[0]);
  }
}
}  // namespace

TEST(ScalarClKernels, scalar_sum_builds_as_opencl_1_2) {
  build_cl12({stan::math::opencl_kernels::scalar_sum_kernel_code});
}

TEST(ScalarClKernels, scalar_expression_builds_as_opencl_1_2) {
  matrix_cl<double> m_cl(MatrixXd::Ones(3, 2));
  matrix_cl<double> res_cl(3, 2);
  ScalarCl<double> s(2.0);
  auto expr = m_cl * s + s;
  build_cl12({stan::math::view_kernel_helpers,
              expr.get_kernel_source_for_evaluating_into(res_cl)});
  ScalarCl<double> t;
  auto scalar_expr
      = stan::math::scalar_result_<decltype(stan::math::as_operation_cl(s)
                                            * 3.0)>(
          stan::math::as_operation_cl(s) * 3.0);
  build_cl12({stan::math::view_kernel_helpers,
              scalar_expr.get_kernel_source_for_evaluating_into(t.matrix())});
}

TEST(ScalarClKernels, handwritten_kernel_takes_scalar_cl) {
  MatrixXd m = MatrixXd::Random(5, 3);
  matrix_cl<double> m_cl(m);
  ScalarCl<double> out(100.0);
  stan::math::opencl_kernels::scalar_sum(cl::NDRange(64), cl::NDRange(64), out,
                                         m_cl, m_cl.size(), 0.5, 0);
  EXPECT_NEAR(to_host(out), m.sum() + 0.5, 1e-12);
  stan::math::opencl_kernels::scalar_sum(cl::NDRange(64), cl::NDRange(64), out,
                                         m_cl, m_cl.size(), 0.0, 1);
  EXPECT_NEAR(to_host(out), 2 * m.sum() + 0.5, 1e-12);
}

TEST(ScalarClKernels, sum_into_sizes) {
  for (auto dims : std::vector<std::pair<int, int>>{{0, 0},
                                                    {0, 5},
                                                    {1, 1},
                                                    {7, 1},
                                                    {1, 7},
                                                    {3, 1000},
                                                    {1000, 3},
                                                    {100, 100},
                                                    {1000, 1000}}) {
    MatrixXd m = MatrixXd::Random(dims.first, dims.second);
    matrix_cl<double> m_cl(m);
    ScalarCl<double> res(-7.0);
    sum_into(res, m_cl, false);
    EXPECT_NEAR(to_host(res), m.sum(), 1e-9)
        << dims.first << "x" << dims.second;
  }
}

TEST(ScalarClKernels, sum_into_accumulate_and_offset) {
  MatrixXd m = MatrixXd::Random(40, 30);
  matrix_cl<double> m_cl(m);
  ScalarCl<double> res(1.5);
  sum_into(res, m_cl, true);
  EXPECT_NEAR(to_host(res), 1.5 + m.sum(), 1e-10);
  sum_into(res, m_cl, true, 2.0);
  EXPECT_NEAR(to_host(res), 1.5 + 2 * m.sum() + 2.0, 1e-10);
  sum_into(res, m_cl, false, -1.0);
  EXPECT_NEAR(to_host(res), m.sum() - 1.0, 1e-10);

  matrix_cl<double> empty_cl(0, 3);
  sum_into(res, empty_cl, true, 4.0);
  EXPECT_NEAR(to_host(res), m.sum() + 3.0, 1e-10);
  sum_into(res, empty_cl, false, 4.0);
  EXPECT_EQ(to_host(res), 4.0);
}

TEST(ScalarClKernels, sum_into_expression) {
  MatrixXd m = MatrixXd::Random(200, 13);
  matrix_cl<double> m_cl(m);
  ScalarCl<double> s(3.0);
  ScalarCl<double> res;
  sum_into(res, m_cl * s + 1.0, false);
  EXPECT_NEAR(to_host(res), (m.array() * 3.0 + 1.0).sum(), 1e-9);
}

TEST(ScalarClKernels, generated_handwritten_generated_chain) {
  MatrixXd m = MatrixXd::Random(50, 20);
  matrix_cl<double> m_cl(m);
  ScalarCl<double> s(2.0);
  for (int i = 0; i < 10; ++i) {
    // generated kernel writes s
    s = s * 0.5 + 1.0;
    // handwritten kernel reads m and accumulates into s
    sum_into(s, m_cl, true);
  }
  // generated kernel reads s
  matrix_cl<double> res_cl = m_cl * s;
  double expected = 2.0;
  for (int i = 0; i < 10; ++i) {
    expected = expected * 0.5 + 1.0 + m.sum();
  }
  EXPECT_NEAR(to_host(s), expected, 1e-9);
  EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(res_cl), m * expected, 1e-9);
}
#endif
