#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/opencl.hpp>
#include <gtest/gtest.h>
#include <limits>

using stan::math::check_nan;

// ---------- check_nan: matrix tests ----------
TEST(ErrorHandlingScalarCL, check_nan_Matrix) {
  const char* function = "check_nan";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;
  using stan::math::matrix_cl;
  x.resize(3);
  x << -1, 0, 1;
  matrix_cl<double> xx(x);

  ASSERT_NO_THROW(check_nan(function, "xx", xx))
      << "check_nan should be true with finite xx";
}

TEST(ErrorHandlingScalarCL, check_nan_Matrix_quit_nan) {
  const char* function = "check_nan";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;
  using stan::math::matrix_cl;

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::quiet_NaN();
  matrix_cl<double> xx(x);
  EXPECT_THROW(check_nan(function, "xx", xx), std::domain_error)
      << "check_nan should throw exception on NaN";
}

TEST(ErrorHandlingScalarCL, check_nan_positions) {
  const char* function = "check_nan";
  double nan = std::numeric_limits<double>::quiet_NaN();
  using stan::math::matrix_cl;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x_mat(3);
  x_mat << nan, 0, 1;
  matrix_cl<double> xx_mat1(x_mat);
  EXPECT_THROW(check_nan(function, "xx_mat1", xx_mat1), std::domain_error);

  x_mat << 1, nan, 1;
  matrix_cl<double> xx_mat2(x_mat);
  EXPECT_THROW(check_nan(function, "xx_mat2", xx_mat2), std::domain_error);

  x_mat << 1, 0, nan;
  matrix_cl<double> xx_mat3(x_mat);
  EXPECT_THROW(check_nan(function, "xx_mat3", xx_mat3), std::domain_error);
}

TEST(ErrorHandlingScalarCL, check_rv_v_symmetric_cl) {
  const char* function = "check_symmetric";

  stan::math::row_vector_d rv;
  stan::math::vector_d v;
  rv.resize(3);
  v.resize(3);
  stan::math::matrix_cl<double> rvv(rv);
  stan::math::matrix_cl<double> vv(v);
  EXPECT_THROW(check_symmetric(function, "rvv_non_symm", rvv),
               std::invalid_argument);
  EXPECT_THROW(check_symmetric(function, "v_non_symm", vv),
               std::invalid_argument);
}

TEST(ErrorHandlingScalarCL, check_m_symmetric) {
  const char* function = "check_symmetric";

  stan::math::matrix_d m_ok(3, 3);
  stan::math::matrix_d m_fail(3, 3);
  m_fail << 1, 2, 3, 2, 4, -5, 3, 5, 6;
  stan::math::matrix_cl<double> mm_fail(m_fail);
  m_ok << 1, 2, 3, 2, 4, 5, 3, 5, 6;
  stan::math::matrix_cl<double> mm_ok(m_ok);
  EXPECT_THROW(check_symmetric(function, "m_non_symm", mm_fail),
               std::domain_error);
  EXPECT_NO_THROW(check_symmetric(function, "m_symm_mat1", mm_ok));
}

TEST(ErrorHandlingScalarCL, check_m_diagonal_zeros) {
  const char* function = "check_diagonal_zeros";

  stan::math::matrix_d m_ok(3, 3);
  stan::math::matrix_d m_fail(3, 3);
  m_ok << 1, 2, 3, 2, 4, -5, 3, 5, 6;
  m_fail << 1, 2, 3, 2, 0, -5, 3, 5, 6;
  stan::math::matrix_cl<double> mm_ok(m_ok);
  stan::math::matrix_cl<double> mm_fail(m_fail);
  EXPECT_NO_THROW(check_diagonal_zeros(function, "mm_ok", mm_ok));
  EXPECT_THROW(check_diagonal_zeros(function, "mm_fail", mm_fail),
               std::domain_error);
}

TEST(ErrorHandlingScalarCL, check_m_mat_size_one) {
  const char* function = "check_mat_size_one";

  stan::math::matrix_d m_ok(1, 1);
  stan::math::matrix_d m_fail0;
  stan::math::matrix_d m_fail33(3, 3);
  stan::math::matrix_d m_fail31(3, 1);
  stan::math::matrix_cl<double> mm_ok(m_ok);
  stan::math::matrix_cl<double> mm_fail0(m_fail0);
  stan::math::matrix_cl<double> mm_fail33(m_fail33);
  stan::math::matrix_cl<double> mm_fail31(m_fail31);
  EXPECT_NO_THROW(check_mat_size_one(function, "mm_ok", mm_ok));
  EXPECT_THROW(check_mat_size_one(function, "mm_fail0", mm_fail0),
               std::invalid_argument);
  EXPECT_THROW(check_mat_size_one(function, "mm_fail33", mm_fail33),
               std::invalid_argument);
  EXPECT_THROW(check_mat_size_one(function, "mm_fail31", mm_fail31),
               std::invalid_argument);
}

TEST(ErrorHandlingScalarCL, check_m_mat_not_size_one) {
  const char* function = "check_mat_not_size_one";

  stan::math::matrix_d m_ok31(3, 1);
  stan::math::matrix_d m_ok0;
  stan::math::matrix_d m_ok33(3, 3);
  stan::math::matrix_d m_fail(1, 1);
  stan::math::matrix_cl<double> mm_ok31(m_ok31);
  stan::math::matrix_cl<double> mm_ok0(m_ok0);
  stan::math::matrix_cl<double> mm_ok33(m_ok33);
  stan::math::matrix_cl<double> mm_fail(m_fail);
  EXPECT_NO_THROW(check_mat_not_size_one(function, "mm_ok31", mm_ok31));
  EXPECT_NO_THROW(check_mat_not_size_one(function, "mm_ok0", mm_ok0));
  EXPECT_NO_THROW(check_mat_not_size_one(function, "mm_ok33", mm_ok33));
  EXPECT_THROW(check_mat_not_size_one(function, "mm_fail", mm_fail),
               std::invalid_argument);
}

TEST(ErrorHandlingScalarCL, check_m_invalid_matrix_view) {
  const char* function = "check_invalid_matrix_view";

  stan::math::matrix_cl<double> m(5, 5);
  m.view(stan::math::matrix_cl_view::Lower);
  EXPECT_THROW(check_invalid_matrix_view(function, "m", m,
                                         stan::math::matrix_cl_view::Lower),
               std::invalid_argument);
  m.view(stan::math::matrix_cl_view::Upper);
  EXPECT_THROW(check_invalid_matrix_view(function, "m", m,
                                         stan::math::matrix_cl_view::Upper),
               std::invalid_argument);
  EXPECT_NO_THROW(check_invalid_matrix_view(function, "m", m,
                                            stan::math::matrix_cl_view::Lower));
}
#endif
