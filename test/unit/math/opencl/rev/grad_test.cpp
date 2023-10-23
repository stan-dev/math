#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(OpenCLGradTest, exceptions) {
  using stan::math::matrix_cl;
  using stan::math::to_matrix_cl;
  using stan::math::var;
  using stan::math::var_value;
  using varmat_cl = var_value<matrix_cl<double>>;
  Eigen::VectorXd a(6);
  a << 1, 2, 3, 4, 5, 6;
  varmat_cl a_cl(to_matrix_cl(a));
  varmat_cl b_cl(a_cl.block(0, 0, 3, 1));
  varmat_cl c_cl(a_cl.block(3, 0, 3, 1));
  Eigen::VectorXd ret_grads_cl(6);
  using stan::math::add;
  using stan::math::subtract;
  var ret = stan::math::sum(add(b_cl, b_cl));
  stan::math::grad(ret, a_cl, ret_grads_cl);
  std::cout << "opencl grads: \n" << ret_grads_cl << std::endl;
  stan::math::recover_memory();
  Eigen::Matrix<var, -1, 1> a_host = a;
  Eigen::Matrix<var, -1, 1> b_host = a_host.segment(0, 3);
  Eigen::Matrix<var, -1, 1> c_host = a_host.segment(3, 3);
  var ret_host = stan::math::sum(add(b_host, b_host));
  Eigen::VectorXd ret_grads_host(6);
  stan::math::grad(ret_host, a_host, ret_grads_host);
  std::cout << "host grads: \n" << ret_grads_host << std::endl;
  stan::math::recover_memory();
}

#endif
