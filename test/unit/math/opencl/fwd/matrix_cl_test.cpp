#ifdef STAN_OPENCL
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/fwd/mat/fun/typedefs.hpp>
#include <stan/math/mix/mat.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/fwd/matrix_cl.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixCL, matrix_cl_fvar_creation) {
  stan::math::vector_fd vec_1;
  stan::math::matrix_fd mat_1;
  stan::math::matrix_fd mat_2;

  vec_1.resize(3);
  mat_1.resize(2, 3);
  mat_2.resize(3, 3);
  vec_1 << 1, 2, 3;
  mat_1 << 1, 2, 3, 1, 2, 3;
  mat_2 << 1, 2, 3, 1, 2, 3, 1, 2, 3;

  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::fvar<double>> A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::fvar<double>> d22(mat_1));
  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::fvar<double>> d33(mat_2));
  EXPECT_NO_THROW(stan::math::matrix_cl<stan::math::fvar<double>> d11(vec_1));
}


TEST(MathMatrixCL, matrix_cl_fvar_var_creation) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::matrix_cl;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  typedef Matrix<fvar<var>, Dynamic, Dynamic> matrix_fv;
  typedef Matrix<fvar<var>, Dynamic, 1> vector_fv;
  vector_fv vec_1(10);
  matrix_fv mat_1(5, 2);
  matrix_fv mat_2(10, 1);

  auto filler = 1.2;
  for (size_t i = 0; i < 10; ++i) {
    vec_1[i] = filler;
    vec_1[i].d_ = filler;
    mat_1(i) = filler;
    mat_1(i).d_ = filler;
    mat_2(i) = filler;
    mat_2(i).d_ = filler;
    filler++;
  }

  EXPECT_NO_THROW(matrix_cl<fvar<var>> A(1, 1));
  EXPECT_NO_THROW(matrix_cl<fvar<var>> d22(mat_1));
  EXPECT_NO_THROW(matrix_cl<fvar<var>> d33(mat_2));
  EXPECT_NO_THROW(matrix_cl<fvar<var>> d11(vec_1));
}

TEST(MathMatrixCL, matrix_cl_fvar_fvar_double_creation) {
  using stan::math::fvar;
  using stan::math::matrix_cl;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  typedef Matrix<fvar<fvar<double>>, Dynamic, Dynamic> matrix_fv;
  typedef Matrix<fvar<fvar<double>>, Dynamic, 1> vector_fv;
  typedef matrix_cl<fvar<fvar<double>>> matrix_ffv_cl;
  vector_fv vec_1(10);
  matrix_fv mat_1(5, 2);
  matrix_fv mat_2(10, 1);

  auto filler = 1.2;
  for (size_t i = 0; i < 10; ++i) {
    vec_1[i] = filler;
    vec_1[i].d_.val_ = filler;
    vec_1[i].d_.d_ = filler;
    mat_1(i) = filler;
    mat_1(i).d_.val_ = filler;
    mat_1(i).d_.d_ = filler;
    mat_2(i) = filler;
    mat_2(i).d_.val_ = filler;
    mat_2(i).d_.d_ = filler;
    filler++;
  }

  EXPECT_NO_THROW(matrix_ffv_cl A(1, 1));
  EXPECT_NO_THROW(matrix_ffv_cl d11(vec_1));
  EXPECT_NO_THROW(matrix_ffv_cl d22(mat_1));
  EXPECT_NO_THROW(matrix_ffv_cl d33(mat_2));
}

TEST(MathMatrixCL, matrix_cl_fvar_fvar_var_creation) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::matrix_cl;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  typedef Matrix<fvar<fvar<var>>, Dynamic, Dynamic> matrix_fv;
  typedef Matrix<fvar<fvar<var>>, Dynamic, 1> vector_fv;
  typedef matrix_cl<fvar<fvar<var>>> matrix_ffv_cl;
  vector_fv vec_1(10);
  matrix_fv mat_1(5, 2);
  matrix_fv mat_2(10, 1);

  auto filler = 1.2;
  for (size_t i = 0; i < 10; ++i) {
    vec_1[i] = filler;
    vec_1[i].d_ = filler;
    vec_1[i].d_.d_ = filler;
    mat_1(i) = filler;
    mat_1(i).d_ = filler;
    mat_1(i).d_.d_ = filler;
    mat_2(i) = filler;
    mat_2(i).d_ = filler;
    mat_2(i).d_.d_ = filler;
    filler++;
  }

  EXPECT_NO_THROW(matrix_ffv_cl A(1, 1));
  EXPECT_NO_THROW(matrix_ffv_cl d11(vec_1));
  EXPECT_NO_THROW(matrix_ffv_cl d22(mat_1));
  EXPECT_NO_THROW(matrix_ffv_cl d33(mat_2));
}

#endif
