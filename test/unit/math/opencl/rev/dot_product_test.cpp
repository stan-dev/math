#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto dot_product_functor = [](const auto& a, const auto& b) {
  return stan::math::dot_product(a, b);
};

TEST(OpenCLDotProduct, errors) {
  Eigen::VectorXd a(3);
  Eigen::VectorXd b(8);
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);
  EXPECT_THROW(stan::math::dot_product(a_cl, b_cl), std::invalid_argument);
  EXPECT_THROW(stan::math::dot_product(a_cl, b_cl), std::invalid_argument);
}

TEST(OpenCLDotProduct, prim_rev_small_vector) {
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  Eigen::VectorXd b(8);
  b << 1, 2, 3, 4, 5, 6, 7, 8;
  stan::math::test::compare_cpu_opencl_prim_rev(dot_product_functor, a, b);
}

TEST(OpenCLDotProduct, prim_rev_small_row_vector) {
  Eigen::RowVectorXd a(8);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  Eigen::RowVectorXd b(8);
  b << 1, 2, 3, 4, 5, 6, 7, 8;
  stan::math::test::compare_cpu_opencl_prim_rev(dot_product_functor, a, b);
}

TEST(OpenCLDotProduct, prim_rev_size_0) {
  int N = 0;

  Eigen::VectorXd a(N);
  stan::math::test::compare_cpu_opencl_prim_rev(dot_product_functor, a, a);
}

TEST(OpenCLDotProduct, prim_rev_values_large) {
  int N = 71;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(dot_product_functor, a, b);
}

#endif
