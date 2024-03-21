#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto add_functor
    = [](const auto& a, const auto& b) { return stan::math::add(a, b).eval(); };

TEST(OpenCLPrim, add_v_small_zero) {
  stan::math::vector_d d1(3), d2(3);
  d1 << 1, 2, 3;
  d2 << 3, 2, 1;
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d1, d2);

  stan::math::vector_d d0(0);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d0, d0);

  double d3 = 3.0;
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d1, d3);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d3, d1);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d3, d0);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d0, d3);
}

TEST(OpenCLPrim, add_rv_small_zero) {
  stan::math::row_vector_d d1(3), d2(3);
  d1 << 1, 2, 3;
  d2 << 3, 2, 1;
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d1, d2);

  stan::math::vector_d d0(0);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d0, d0);

  double d3 = 3.0;
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d1, d3);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d3, d1);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d3, d0);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d0, d3);
}

TEST(OpenCLPrim, add_m_small_zero) {
  stan::math::matrix_d d1(3, 3), d2(3, 3);
  d1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  d2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d1, d2);

  stan::math::matrix_d d0(0, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d0, d0);

  double d3 = 3.0;
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d1, d3);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d3, d1);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d3, d0);
  stan::math::test::compare_cpu_opencl_prim_rev(add_functor, d0, d3);
}

TEST(OpenCLPrim, add_rev_exceptions) {
  using stan::math::matrix_cl;
  stan::math::vector_d vd1(2), vd2(3);
  stan::math::var_value<matrix_cl<double>> vd11 = stan::math::to_matrix_cl(vd1);
  stan::math::var_value<matrix_cl<double>> vd22 = stan::math::to_matrix_cl(vd2);
  EXPECT_THROW(stan::math::add(vd11, vd22), std::invalid_argument);

  stan::math::row_vector_d rvd1(2), rvd2(3);
  stan::math::var_value<matrix_cl<double>> rvd11
      = stan::math::to_matrix_cl(rvd1);
  stan::math::var_value<matrix_cl<double>> rvd22
      = stan::math::to_matrix_cl(rvd2);
  EXPECT_THROW(stan::math::add(rvd11, rvd22), std::invalid_argument);

  stan::math::matrix_d md1(2, 2), md2(3, 3);
  stan::math::var_value<matrix_cl<double>> md11 = stan::math::to_matrix_cl(md1);
  stan::math::var_value<matrix_cl<double>> md22 = stan::math::to_matrix_cl(md2);
  EXPECT_THROW(stan::math::add(md11, md22), std::invalid_argument);
}

TEST(OpenCLPrim, add_aliasing) {
  stan::math::matrix_d d1(3, 3);
  d1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  using stan::math::matrix_cl;
  using stan::math::var;
  using stan::math::var_value;
  using varmat_cl = var_value<matrix_cl<double>>;
  varmat_cl d11 = stan::math::to_matrix_cl(d1);
  // Add the same matrix as the left and right hand side
  var res = stan::math::sum(stan::math::add(d11, d11));
  res.grad();
  // Get back adjoints
  Eigen::MatrixXd grad_res = stan::math::from_matrix_cl(d11.adj());
  stan::math::recover_memory();
  Eigen::Matrix<var, -1, -1> d_host = d1;
  // Same op as above but on the host
  var res_host = stan::math::sum(stan::math::add(d_host, d_host));
  res_host.grad();
  Eigen::MatrixXd grad_res_host = d_host.adj();
  std::cout << "OpenCL Adjoints: " << std::endl;
  std::cout << grad_res << std::endl;
  std::cout << "CPU Adjoints: " << std::endl;
  std::cout << grad_res_host << std::endl;
}

#endif
