#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/rev/opencl.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixRevCL, triangular_transpose_m_exception_pass) {
  using stan::math::TriangularMapCL;
  using stan::math::matrix_cl;
  using stan::math::var;
  matrix_cl<var> m1(1, 1);
  matrix_cl<var> m0;
  matrix_cl<var> m2(5, 3);
  matrix_cl<var> m3(3, 4);

  EXPECT_NO_THROW(m0.triangular_transpose<TriangularMapCL::LowerToUpper>());
  EXPECT_NO_THROW(m0.triangular_transpose<TriangularMapCL::UpperToLower>());
  EXPECT_NO_THROW(m1.triangular_transpose<TriangularMapCL::LowerToUpper>());
  EXPECT_NO_THROW(m1.triangular_transpose<TriangularMapCL::UpperToLower>());
  EXPECT_THROW(m2.triangular_transpose<TriangularMapCL::LowerToUpper>(),
               std::invalid_argument);
  EXPECT_THROW(m2.triangular_transpose<TriangularMapCL::UpperToLower>(),
               std::invalid_argument);
  EXPECT_THROW(m3.triangular_transpose<TriangularMapCL::LowerToUpper>(),
               std::invalid_argument);
  EXPECT_THROW(m3.triangular_transpose<TriangularMapCL::UpperToLower>(),
               std::invalid_argument);
}

TEST(MathMatrixRevCL, triangular_transpose_m_pass) {
  using stan::math::TriangularMapCL;
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::var;
  matrix_v m0(2, 2);
  matrix_d m0_dst(2, 2);
  m0 << 1, 2, 3, 4;

  matrix_cl<var> m00(m0);
  matrix_cl<var> m11(m0);

  EXPECT_NO_THROW(m00.triangular_transpose<TriangularMapCL::LowerToUpper>());
  EXPECT_NO_THROW(m0_dst = from_matrix_cl(m00.val()));
  EXPECT_EQ(1, m0_dst(0, 0));
  EXPECT_EQ(3, m0_dst(0, 1));
  EXPECT_EQ(3, m0_dst(1, 0));
  EXPECT_EQ(4, m0_dst(1, 1));

  EXPECT_NO_THROW(m11.triangular_transpose<TriangularMapCL::UpperToLower>());
  EXPECT_NO_THROW(m0_dst = from_matrix_cl(m11.val()));
  EXPECT_EQ(1, m0_dst(0, 0));
  EXPECT_EQ(2, m0_dst(0, 1));
  EXPECT_EQ(2, m0_dst(1, 0));
  EXPECT_EQ(4, m0_dst(1, 1));
}

#endif
