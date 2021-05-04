#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(OpenCLAdjointResults, errors) {
  stan::math::matrix_cl<double> a(3, 2);
  stan::math::matrix_cl<double> b(6, 1);
  stan::math::var_value<stan::math::matrix_cl<double>> c(
      stan::math::matrix_cl<double>(3, 2));
  stan::math::var d = 0;
  stan::math::var_value<stan::math::matrix_cl<double>> e(
      stan::math::matrix_cl<double>(2, 3));

  EXPECT_THROW(
      stan::math::adjoint_results(c, d) += stan::math::expressions(b, a),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::adjoint_results(c, e) += stan::math::expressions(a, e.val()),
      std::invalid_argument);
  EXPECT_NO_THROW(stan::math::adjoint_results() += stan::math::expressions());
}

TEST(OpenCLAdjointResults, matrix_expressions) {
  Eigen::MatrixXd a(3, 2);
  a << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::var_value<stan::math::matrix_cl<double>> b(
      stan::math::matrix_cl<double>(3, 2));
  b.adj() = a_cl * 0.5;
  stan::math::var c = 0;
  c.adj() = 2;
  double d = 5;

  stan::math::adjoint_results(b, c, d, a_cl)
      += stan::math::expressions(a_cl + 1, a_cl - 1, a_cl + 2, a_cl - 2);

  EXPECT_MATRIX_EQ(stan::math::from_matrix_cl(b.adj()),
                   a.array() * 0.5 + a.array() + 1);
  EXPECT_EQ(c.adj(), 2 + stan::math::sum(a.array() - 1));
  EXPECT_MATRIX_EQ(stan::math::from_matrix_cl(a_cl), a);
  EXPECT_EQ(d, 5);
}

TEST(OpenCLAdjointResults, vector_expressions) {
  Eigen::VectorXd a(6);
  a << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::var_value<stan::math::matrix_cl<double>> b(
      stan::math::matrix_cl<double>(6, 1));
  b.adj() = a_cl * 0.5;
  stan::math::var c = 0;
  c.adj() = 2;
  double d = 5;

  stan::math::adjoint_results(b, c, d, a_cl)
      += stan::math::expressions(a_cl + 1, a_cl - 1, a_cl + 2, a_cl - 2);

  EXPECT_MATRIX_EQ(stan::math::from_matrix_cl(b.adj()),
                   a.array() * 0.5 + a.array() + 1);
  EXPECT_EQ(c.adj(), 2 + stan::math::sum(a.array() - 1));
  EXPECT_MATRIX_EQ(stan::math::from_matrix_cl(a_cl), a);
  EXPECT_EQ(d, 5);
}

#endif
