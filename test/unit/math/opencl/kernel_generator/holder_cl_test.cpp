#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;



template <typename T>
auto f2(T&& a) {
  return stan::math::make_holder([](const auto& a) { return a + a; },
                                    std::forward<T>(a));
}

TEST(KernelGenerator, make_holder_cl_lvalue_test) {
  MatrixXd m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  matrix_cl<double> m_cl(m);
  matrix_cl<double> res_cl = f2(m_cl);

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m + m;
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, make_holder_cl_rvalue_test) {
  MatrixXd m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  matrix_cl<double> m_cl(m);
  matrix_cl<double> res_cl = f2(std::move(m_cl));

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m + m;
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

#endif
