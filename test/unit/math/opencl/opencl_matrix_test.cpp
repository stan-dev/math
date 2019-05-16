#ifdef STAN_OPENCL
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/fwd/mat/fun/typedefs.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixCL, matrix_cl_double_creation) {
  stan::math::vector_d d1;
  stan::math::matrix_d d2;
  stan::math::matrix_d d3;

  d1.resize(3);
  d2.resize(2, 3);
  EXPECT_NO_THROW(stan::math::matrix_cl<double> A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_cl<double> d11(d1));
  EXPECT_NO_THROW(stan::math::matrix_cl<double> d22(d2));
  EXPECT_NO_THROW(stan::math::matrix_cl<double> d33(d3));
}

TEST(MathMatrixCL, matrix_cl_var_creation) {
  stan::math::vector_v d1;
  stan::math::matrix_v d2;
  stan::math::matrix_v d3;

  d1.resize(3);
  d2.resize(2, 3);
  d3.resize(3, 3);
  d1 << 1, 2, 3;
  d2 << 1, 2, 3, 1, 2, 3;
  d3 << 1, 2, 3, 1, 2, 3, 1, 2, 3;

  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::var> A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::var> d22(d2));
  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::var> d33(d3));
  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::var> d11(d1));
}

TEST(MathMatrixCL, matrix_cl_fvar_creation) {
  stan::math::vector_fd d1;
  stan::math::matrix_fd d2;
  stan::math::matrix_fd d3;

  d1.resize(3);
  d2.resize(2, 3);
  d3.resize(3, 3);
  d1 << 1, 2, 3;
  d2 << 1, 2, 3, 1, 2, 3;
  d3 << 1, 2, 3, 1, 2, 3, 1, 2, 3;

  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::fvar<double>> A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::fvar<double>> d22(d2));
  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::fvar<double>> d33(d3));
  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::fvar<double>> d11(d1));
}

/* This doesn't work yet.
TEST(MathMatrixCL, matrix_cl_fvar_var_creation) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::matrix_cl;

  typedef Eigen::Matrix<fvar<var>, Eigen::Dynamic, Eigen::Dynamic> matrix_fv;
  typedef Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> vector_fv;
  vector_fv d1;
  matrix_fv d2;
  matrix_fv d3;

  d1.resize(3);
  d2.resize(2, 3);
  d3.resize(3, 3);
  d1 << 1, 2, 3;
  d2 << 1, 2, 3, 1, 2, 3;
  d3 << 1, 2, 3, 1, 2, 3, 1, 2, 3;

  EXPECT_NO_THROW(matrix_cl<fvar<var>> A(1, 1));
  EXPECT_NO_THROW(matrix_cl<fvar<var>> d22(d2));
  EXPECT_NO_THROW(matrix_cl<fvar<var>> d33(d3));
  EXPECT_NO_THROW(matrix_cl<fvar<var>> d11(d1));
}

TEST(MathMatrixCL, matrix_cl_fvar_fvar_creation) {
  stan::math::vector_ffd d1;
  stan::math::matrix_ffd d2;
  stan::math::matrix_ffd d3;

  d1.resize(3);
  d2.resize(2, 3);
  d3.resize(3, 3);
  d1 << 1, 2, 3;
  d2 << 1, 2, 3, 1, 2, 3;
  d3 << 1, 2, 3, 1, 2, 3, 1, 2, 3;
  using stan::math::fvar;
  using stan::math::matrix_cl;
  EXPECT_NO_THROW(matrix_cl<fvar<fvar<double>>> A(1, 1));
//  EXPECT_NO_THROW(matrix_cl<fvar<fvar<double>>> d22(d2));
  EXPECT_NO_THROW(matrix_cl<fvar<fvar<double>>> d33(d3));
  EXPECT_NO_THROW(matrix_cl<fvar<fvar<double>>> d11(d1));
}
*/
#endif
