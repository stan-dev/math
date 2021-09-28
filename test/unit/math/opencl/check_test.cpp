#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

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
#endif
