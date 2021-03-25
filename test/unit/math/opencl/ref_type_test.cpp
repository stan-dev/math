#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, ref_type_matrix_cl) {
  using stan::ref_type_t;
  using stan::math::matrix_cl;
  Eigen::MatrixXd a_eig(3, 3);
  a_eig << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  matrix_cl<double> a(a_eig);
  matrix_cl<double> a2 = a;
  ref_type_t<matrix_cl<double>> a_ref1 = a;
  ref_type_t<matrix_cl<double>&> a_ref2 = a;
  ref_type_t<matrix_cl<double>&&> a_ref3 = std::move(a2);

  EXPECT_MATRIX_EQ(from_matrix_cl(a_ref1), from_matrix_cl(a));
  EXPECT_MATRIX_EQ(from_matrix_cl(a_ref2), from_matrix_cl(a));
  EXPECT_MATRIX_EQ(from_matrix_cl(a_ref3), from_matrix_cl(a));

  EXPECT_TRUE((std::is_same<decltype(a), ref_type_t<decltype(a)&&>>::value));
  EXPECT_TRUE(std::is_lvalue_reference<ref_type_t<matrix_cl<double>>>::value);
  EXPECT_TRUE(std::is_lvalue_reference<ref_type_t<matrix_cl<double>&>>::value);
  EXPECT_FALSE(std::is_reference<ref_type_t<matrix_cl<double>&&>>::value);
}

TEST(MathMetaPrim, ref_type_kg_light_expression) {
  using stan::ref_type_t;
  using stan::math::matrix_cl;
  Eigen::MatrixXd m_eig(3, 3);
  m_eig << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  matrix_cl<double> m(m_eig);
  auto a = stan::math::block_zero_based(m, 0, 1, 2, 2);
  ref_type_t<decltype(a)> a_ref1 = a;
  ref_type_t<decltype(a)&> a_ref2 = a;
  ref_type_t<decltype(a)&&> a_ref3
      = stan::math::block_zero_based(m, 0, 1, 2, 2);

  matrix_cl<double> a_eval = a;
  EXPECT_MATRIX_EQ(from_matrix_cl(a_ref1), from_matrix_cl(a_eval));
  EXPECT_MATRIX_EQ(from_matrix_cl(a_ref2), from_matrix_cl(a_eval));
  EXPECT_MATRIX_EQ(from_matrix_cl(a_ref3), from_matrix_cl(a_eval));

  EXPECT_FALSE((stan::is_matrix_cl<ref_type_t<decltype(a)>>::value));
  EXPECT_FALSE((stan::is_matrix_cl<ref_type_t<decltype(a)&>>::value));
  EXPECT_FALSE((stan::is_matrix_cl<ref_type_t<decltype(a)&&>>::value));
}

TEST(MathMetaPrim, ref_type_kg_heavy_expression) {
  using stan::ref_type_t;
  using stan::math::matrix_cl;
  Eigen::MatrixXd m_eig(3, 3);
  m_eig << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  matrix_cl<double> m(m_eig);
  auto a = m * 3;
  ref_type_t<decltype(a)> a_ref1 = a;
  ref_type_t<decltype(a)&> a_ref2 = a;
  ref_type_t<decltype(a)&&> a_ref3 = m * 3;

  matrix_cl<double> a_eval = a;
  EXPECT_MATRIX_EQ(from_matrix_cl(a_ref1), from_matrix_cl(a_eval));
  EXPECT_MATRIX_EQ(from_matrix_cl(a_ref2), from_matrix_cl(a_eval));
  EXPECT_MATRIX_EQ(from_matrix_cl(a_ref3), from_matrix_cl(a_eval));

  EXPECT_TRUE((std::is_same<stan::math::matrix_cl<double>,
                            ref_type_t<decltype(a)>>::value));
  EXPECT_TRUE((std::is_same<stan::math::matrix_cl<double>,
                            ref_type_t<decltype(a)&>>::value));
  EXPECT_TRUE((std::is_same<stan::math::matrix_cl<double>,
                            ref_type_t<decltype(a)&&>>::value));
}

#endif
