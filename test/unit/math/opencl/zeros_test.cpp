#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixCL, zero_m_exception_pass) {
  stan::math::matrix_cl<double> m(1, 1);

  EXPECT_NO_THROW(m.zeros<stan::math::matrix_cl_view::Entire>());
  EXPECT_NO_THROW(m.zeros<stan::math::matrix_cl_view::Diagonal>());
  EXPECT_NO_THROW(m.zeros<stan::math::matrix_cl_view::Lower>());
  EXPECT_NO_THROW(m.zeros<stan::math::matrix_cl_view::Upper>());

  EXPECT_THROW(m.zeros_strict_tri<stan::math::matrix_cl_view::Entire>(),
               std::invalid_argument);
  EXPECT_THROW(m.zeros_strict_tri<stan::math::matrix_cl_view::Diagonal>(),
               std::invalid_argument);
  EXPECT_NO_THROW(m.zeros_strict_tri<stan::math::matrix_cl_view::Lower>());
  EXPECT_NO_THROW(m.zeros_strict_tri<stan::math::matrix_cl_view::Upper>());

  stan::math::matrix_cl<double> m0;
  EXPECT_NO_THROW(m0.zeros<stan::math::matrix_cl_view::Entire>());
  EXPECT_NO_THROW(m0.zeros<stan::math::matrix_cl_view::Diagonal>());
  EXPECT_NO_THROW(m0.zeros<stan::math::matrix_cl_view::Lower>());
  EXPECT_NO_THROW(m0.zeros<stan::math::matrix_cl_view::Upper>());

  EXPECT_THROW(m0.zeros_strict_tri<stan::math::matrix_cl_view::Entire>(),
               std::invalid_argument);
  EXPECT_THROW(m0.zeros_strict_tri<stan::math::matrix_cl_view::Diagonal>(),
               std::invalid_argument);
  EXPECT_NO_THROW(m0.zeros_strict_tri<stan::math::matrix_cl_view::Lower>());
  EXPECT_NO_THROW(m0.zeros_strict_tri<stan::math::matrix_cl_view::Upper>());
}

TEST(MathMatrixCL, zeros_strict_m_value_check) {
  stan::math::matrix_d m0(2, 2);
  stan::math::matrix_d m0_dst(2, 2);
  m0 << 2, 2, 2, 2;
  stan::math::matrix_cl<double> m(m0);
  stan::math::matrix_cl<double> m_upper(m0);
  stan::math::matrix_cl<double> m_lower(m0);

  EXPECT_THROW(m.zeros_strict_tri<stan::math::matrix_cl_view::Entire>(),
               std::invalid_argument);
  EXPECT_THROW(m.zeros_strict_tri<stan::math::matrix_cl_view::Diagonal>(),
               std::invalid_argument);
  EXPECT_NO_THROW(
      m_lower.zeros_strict_tri<stan::math::matrix_cl_view::Lower>());
  EXPECT_NO_THROW(
      m_upper.zeros_strict_tri<stan::math::matrix_cl_view::Upper>());

  m0_dst = stan::math::from_matrix_cl(m_lower);
  EXPECT_EQ(2, m0_dst(0, 0));
  EXPECT_EQ(2, m0_dst(0, 1));
  EXPECT_EQ(0, m0_dst(1, 0));
  EXPECT_EQ(2, m0_dst(1, 1));

  m0_dst = stan::math::from_matrix_cl(m_upper);
  EXPECT_EQ(2, m0_dst(0, 0));
  EXPECT_EQ(0, m0_dst(0, 1));
  EXPECT_EQ(2, m0_dst(1, 0));
  EXPECT_EQ(2, m0_dst(1, 1));
}

TEST(MathMatrixCL, zeros_m_value_check) {
  stan::math::matrix_d m0(2, 2);
  stan::math::matrix_d m0_dst(2, 2);
  m0 << 2, 2, 2, 2;
  stan::math::matrix_cl<double> m(m0);
  stan::math::matrix_cl<double> m_diagonal(m0);
  stan::math::matrix_cl<double> m_upper(m0);
  stan::math::matrix_cl<double> m_lower(m0);

  EXPECT_NO_THROW(m.zeros<stan::math::matrix_cl_view::Entire>());
  EXPECT_NO_THROW(m_diagonal.zeros<stan::math::matrix_cl_view::Diagonal>());
  EXPECT_NO_THROW(m_lower.zeros<stan::math::matrix_cl_view::Lower>());
  EXPECT_NO_THROW(m_upper.zeros<stan::math::matrix_cl_view::Upper>());

  m0_dst = stan::math::from_matrix_cl(m);
  EXPECT_EQ(0, m0_dst(0, 0));
  EXPECT_EQ(0, m0_dst(0, 1));
  EXPECT_EQ(0, m0_dst(1, 0));
  EXPECT_EQ(0, m0_dst(1, 1));

  m0_dst = stan::math::from_matrix_cl(m_diagonal);
  EXPECT_EQ(0, m0_dst(0, 0));
  EXPECT_EQ(2, m0_dst(0, 1));
  EXPECT_EQ(2, m0_dst(1, 0));
  EXPECT_EQ(0, m0_dst(1, 1));

  m0_dst = stan::math::from_matrix_cl(m_lower);
  EXPECT_EQ(0, m0_dst(0, 0));
  EXPECT_EQ(2, m0_dst(0, 1));
  EXPECT_EQ(0, m0_dst(1, 0));
  EXPECT_EQ(0, m0_dst(1, 1));

  m0_dst = stan::math::from_matrix_cl(m_upper);
  EXPECT_EQ(0, m0_dst(0, 0));
  EXPECT_EQ(0, m0_dst(0, 1));
  EXPECT_EQ(2, m0_dst(1, 0));
  EXPECT_EQ(0, m0_dst(1, 1));
}

#endif
