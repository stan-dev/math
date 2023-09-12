#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>



TEST(OpenCLPrim, add_aliasing) {
  stan::math::matrix_d d1(3, 3);
  d1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  using stan::math::matrix_cl;
  using stan::math::var_value;
  using stan::math::var;
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
