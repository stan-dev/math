#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>
boost::random::mt19937 rng;

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrix, vector_row_vector) {
  stan::math::vector_d v(3);
  stan::math::row_vector_d rv(3);
  stan::math::matrix_cl<double> v_cl(v);
  stan::math::matrix_cl<double> rv_cl(rv);
  stan::math::matrix_cl<double> m_cl(1, 1);
  EXPECT_NO_THROW(m_cl = rv_cl * v_cl);
}

TEST(MathMatrix, one_dim_zero_matrix) {
  stan::math::matrix_d m0(5, 0);
  stan::math::matrix_d m1(0, 3);

  stan::math::matrix_cl<double> m0_cl(m0);
  stan::math::matrix_cl<double> m1_cl(m1);
  EXPECT_NO_THROW(m0_cl * m1_cl);

  EXPECT_NO_THROW(m0_cl * 2.0);
  EXPECT_NO_THROW(2.0 * m0_cl);

  EXPECT_NO_THROW(m1_cl * 2.0);
  EXPECT_NO_THROW(2.0 * m1_cl);
}

TEST(MathMatrix, zero_result_matrix) {
  stan::math::matrix_d m0(0, 5);
  stan::math::matrix_d m1(5, 0);

  stan::math::matrix_cl<double> m0_cl(m0);
  stan::math::matrix_cl<double> m1_cl(m1);
  EXPECT_NO_THROW(m0_cl * m1_cl);
}

TEST(MathMatrix, zero_size_input_matrix) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d m1(0, 0);

  stan::math::matrix_cl<double> m0_cl(m0);
  stan::math::matrix_cl<double> m1_cl(m1);
  EXPECT_NO_THROW(m0_cl * m1_cl);

  EXPECT_NO_THROW(m0_cl * 2.0);
  EXPECT_NO_THROW(2.0 * m0_cl);
}

TEST(MathMatrix, non_matching_dim_excpetion) {
  stan::math::matrix_d m0(5, 3);
  stan::math::matrix_d m1(2, 6);

  stan::math::matrix_cl<double> m0_cl(m0);
  stan::math::matrix_cl<double> m1_cl(m1);
  EXPECT_THROW(m0_cl * m1_cl, std::invalid_argument);
}

TEST(MathMatrix, multiply_scalar) {
  auto v = stan::math::vector_d::Random(25).eval();
  stan::math::vector_d v_cl_res(25);
  auto rv = stan::math::row_vector_d::Random(25).eval();
  stan::math::row_vector_d rv_cl_res(25);
  auto m = stan::math::matrix_d::Random(5, 5).eval();
  stan::math::matrix_d m_cl_res(5, 5);

  stan::math::matrix_cl<double> v_cl(v);
  v_cl = v_cl * 2.0;
  v_cl_res = stan::math::from_matrix_cl(v_cl);

  stan::math::matrix_cl<double> rv_cl(rv);
  rv_cl = rv_cl * 2.0;
  rv_cl_res = stan::math::from_matrix_cl(rv_cl);

  stan::math::matrix_cl<double> m_cl(m);
  m_cl = m_cl * 2.0;
  m_cl_res = stan::math::from_matrix_cl(m_cl);

  v = v * 2.0;
  rv = rv * 2.0;
  m = m * 2.0;

  EXPECT_MATRIX_NEAR(v, v_cl_res, 1e-10);
  EXPECT_MATRIX_NEAR(rv, rv_cl_res, 1e-10);
  EXPECT_MATRIX_NEAR(m, m_cl_res, 1e-10);
}

TEST(MathMatrix, row_vector_vector) {
  auto v = stan::math::vector_d::Random(5).eval();
  auto rv = stan::math::row_vector_d::Random(5).eval();
  stan::math::matrix_d m0(1, 1);
  stan::math::matrix_d m0_cl_res(1, 1);
  stan::math::matrix_d m1(5, 5);
  stan::math::matrix_d m1_cl_res(5, 5);

  m0 = rv * v;
  m1 = v * rv;

  stan::math::matrix_cl<double> v_cl(v);
  stan::math::matrix_cl<double> rv_cl(rv);
  stan::math::matrix_cl<double> m0_cl(1, 1);
  stan::math::matrix_cl<double> m1_cl(5, 5);

  m0_cl = rv_cl * v_cl;
  m1_cl = v_cl * rv_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  m1_cl_res = stan::math::from_matrix_cl(m1_cl);

  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
  EXPECT_MATRIX_NEAR(m1, m1_cl_res, 1e-10);
}

TEST(MathMatrix, matrix_vector_small) {
  auto m = stan::math::matrix_d::Random(4, 5).eval();
  auto v = stan::math::vector_d::Random(5).eval();
  stan::math::matrix_d m0(4, 1);
  stan::math::matrix_d m0_cl_res(4, 1);

  m0 = m * v;

  stan::math::matrix_cl<double> v_cl(v);
  stan::math::matrix_cl<double> m_cl(m);
  stan::math::matrix_cl<double> m0_cl(4, 1);

  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);

  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
}

TEST(MathMatrix, matrix_vector_big) {
  auto m = stan::math::matrix_d::Random(400, 600).eval();
  auto v = stan::math::vector_d::Random(600).eval();
  stan::math::matrix_d m0(400, 1);
  stan::math::matrix_d m0_cl_res(400, 1);

  m0 = m * v;

  stan::math::matrix_cl<double> v_cl(v);
  stan::math::matrix_cl<double> m_cl(m);
  stan::math::matrix_cl<double> m0_cl(400, 1);

  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);

  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
}

TEST(MathMatrix, row_vector_matrix_small) {
  auto m = stan::math::matrix_d::Random(5, 4).eval();
  auto rv = stan::math::row_vector_d::Random(5).eval();
  stan::math::matrix_d m0(1, 4);
  stan::math::matrix_d m0_cl_res(1, 4);

  m0 = rv * m;

  stan::math::matrix_cl<double> m_cl(m);
  stan::math::matrix_cl<double> rv_cl(rv);
  stan::math::matrix_cl<double> m0_cl(1, 4);

  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);

  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
}

TEST(MathMatrix, row_vector_matrix_big) {
  auto m = stan::math::matrix_d::Random(600, 400).eval();
  auto rv = stan::math::row_vector_d::Random(600).eval();
  stan::math::matrix_d m0(1, 400);
  stan::math::matrix_d m0_cl_res(1, 400);

  m0 = rv * m;

  stan::math::matrix_cl<double> m_cl(m);
  stan::math::matrix_cl<double> rv_cl(rv);
  stan::math::matrix_cl<double> m0_cl(1, 400);

  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);

  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
}

TEST(MathMatrix, matrix_vector_tri_small) {
  auto m = stan::math::matrix_d::Random(4, 5).eval();
  auto v = stan::math::vector_d::Random(5).eval();
  stan::math::matrix_d m0(4, 1);
  stan::math::matrix_d m0_cl_res(4, 1);

  stan::math::matrix_cl<double> v_cl(v);
  stan::math::matrix_cl<double> m_cl(m);
  stan::math::matrix_cl<double> m0_cl(4, 1);

  m0 = m * v.triangularView<Eigen::Lower>();
  m_cl.view(stan::math::matrix_cl_view::Entire);
  v_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = m.triangularView<Eigen::Lower>() * v;
  m_cl.view(stan::math::matrix_cl_view::Lower);
  v_cl.view(stan::math::matrix_cl_view::Entire);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = m * v.triangularView<Eigen::Upper>();
  m_cl.view(stan::math::matrix_cl_view::Entire);
  v_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = m.triangularView<Eigen::Upper>() * v;
  m_cl.view(stan::math::matrix_cl_view::Upper);
  v_cl.view(stan::math::matrix_cl_view::Entire);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  // the cast is needed because operator* for for two triangular views does not
  // exist
  m0 = static_cast<Eigen::MatrixXd>(m.triangularView<Eigen::Lower>())
       * v.triangularView<Eigen::Lower>();
  m_cl.view(stan::math::matrix_cl_view::Lower);
  v_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(m.triangularView<Eigen::Upper>())
       * v.triangularView<Eigen::Upper>();
  m_cl.view(stan::math::matrix_cl_view::Upper);
  v_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(m.triangularView<Eigen::Lower>())
       * v.triangularView<Eigen::Upper>();
  m_cl.view(stan::math::matrix_cl_view::Lower);
  v_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(m.triangularView<Eigen::Upper>())
       * v.triangularView<Eigen::Lower>();
  m_cl.view(stan::math::matrix_cl_view::Upper);
  v_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
}

TEST(MathMatrix, matrix_vector_tri_big) {
  auto m = stan::math::matrix_d::Random(400, 600).eval();
  auto v = stan::math::vector_d::Random(600).eval();
  stan::math::matrix_d m0(400, 1);
  stan::math::matrix_d m0_cl_res(400, 1);

  stan::math::matrix_cl<double> v_cl(v);
  stan::math::matrix_cl<double> m_cl(m);
  stan::math::matrix_cl<double> m0_cl(400, 1);

  m0 = m * v.triangularView<Eigen::Lower>();
  m_cl.view(stan::math::matrix_cl_view::Entire);
  v_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = m.triangularView<Eigen::Lower>() * v;
  m_cl.view(stan::math::matrix_cl_view::Lower);
  v_cl.view(stan::math::matrix_cl_view::Entire);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = m * v.triangularView<Eigen::Upper>();
  m_cl.view(stan::math::matrix_cl_view::Entire);
  v_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = m.triangularView<Eigen::Upper>() * v;
  m_cl.view(stan::math::matrix_cl_view::Upper);
  v_cl.view(stan::math::matrix_cl_view::Entire);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  // the cast is needed because operator* for for two triangular views does not
  // exist
  m0 = static_cast<Eigen::MatrixXd>(m.triangularView<Eigen::Lower>())
       * v.triangularView<Eigen::Lower>();
  m_cl.view(stan::math::matrix_cl_view::Lower);
  v_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(m.triangularView<Eigen::Upper>())
       * v.triangularView<Eigen::Upper>();
  m_cl.view(stan::math::matrix_cl_view::Upper);
  v_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(m.triangularView<Eigen::Lower>())
       * v.triangularView<Eigen::Upper>();
  m_cl.view(stan::math::matrix_cl_view::Lower);
  v_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(m.triangularView<Eigen::Upper>())
       * v.triangularView<Eigen::Lower>();
  m_cl.view(stan::math::matrix_cl_view::Upper);
  v_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = m_cl * v_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
}

TEST(MathMatrix, row_vector_matrix_tri_small) {
  auto m = stan::math::matrix_d::Random(5, 4).eval();
  auto rv = stan::math::row_vector_d::Random(5).eval();
  stan::math::matrix_d m0(1, 4);
  stan::math::matrix_d m0_cl_res(1, 4);

  stan::math::matrix_cl<double> m_cl(m);
  stan::math::matrix_cl<double> rv_cl(rv);
  stan::math::matrix_cl<double> m0_cl(1, 4);

  m0 = rv * m;
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = rv * m.triangularView<Eigen::Lower>();
  rv_cl.view(stan::math::matrix_cl_view::Entire);
  m_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = rv.triangularView<Eigen::Lower>() * m;
  rv_cl.view(stan::math::matrix_cl_view::Lower);
  m_cl.view(stan::math::matrix_cl_view::Entire);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = rv * m.triangularView<Eigen::Upper>();
  rv_cl.view(stan::math::matrix_cl_view::Entire);
  m_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = rv.triangularView<Eigen::Upper>() * m;
  rv_cl.view(stan::math::matrix_cl_view::Upper);
  m_cl.view(stan::math::matrix_cl_view::Entire);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  // the cast is needed because operator* for for two triangular views does not
  // exist
  m0 = static_cast<Eigen::MatrixXd>(rv.triangularView<Eigen::Lower>())
       * m.triangularView<Eigen::Lower>();
  rv_cl.view(stan::math::matrix_cl_view::Lower);
  m_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(rv.triangularView<Eigen::Upper>())
       * m.triangularView<Eigen::Upper>();
  rv_cl.view(stan::math::matrix_cl_view::Upper);
  m_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(rv.triangularView<Eigen::Lower>())
       * m.triangularView<Eigen::Upper>();
  rv_cl.view(stan::math::matrix_cl_view::Lower);
  m_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(rv.triangularView<Eigen::Upper>())
       * m.triangularView<Eigen::Lower>();
  rv_cl.view(stan::math::matrix_cl_view::Upper);
  m_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
}

TEST(MathMatrix, row_vector_matrix_tri_big) {
  auto m = stan::math::matrix_d::Random(600, 400).eval();
  auto rv = stan::math::row_vector_d::Random(600).eval();
  stan::math::matrix_d m0(1, 400);
  stan::math::matrix_d m0_cl_res(1, 400);

  stan::math::matrix_cl<double> m_cl(m);
  stan::math::matrix_cl<double> rv_cl(rv);
  stan::math::matrix_cl<double> m0_cl(1, 400);

  m0 = rv * m;
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = rv * m.triangularView<Eigen::Lower>();
  rv_cl.view(stan::math::matrix_cl_view::Entire);
  m_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = rv.triangularView<Eigen::Lower>() * m;
  rv_cl.view(stan::math::matrix_cl_view::Lower);
  m_cl.view(stan::math::matrix_cl_view::Entire);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = rv * m.triangularView<Eigen::Upper>();
  rv_cl.view(stan::math::matrix_cl_view::Entire);
  m_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = rv.triangularView<Eigen::Upper>() * m;
  rv_cl.view(stan::math::matrix_cl_view::Upper);
  m_cl.view(stan::math::matrix_cl_view::Entire);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  // the cast is needed because operator* for for two triangular views does not
  // exist
  m0 = static_cast<Eigen::MatrixXd>(rv.triangularView<Eigen::Lower>())
       * m.triangularView<Eigen::Lower>();
  rv_cl.view(stan::math::matrix_cl_view::Lower);
  m_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(rv.triangularView<Eigen::Upper>())
       * m.triangularView<Eigen::Upper>();
  rv_cl.view(stan::math::matrix_cl_view::Upper);
  m_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(rv.triangularView<Eigen::Lower>())
       * m.triangularView<Eigen::Upper>();
  rv_cl.view(stan::math::matrix_cl_view::Lower);
  m_cl.view(stan::math::matrix_cl_view::Upper);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);

  m0 = static_cast<Eigen::MatrixXd>(rv.triangularView<Eigen::Upper>())
       * m.triangularView<Eigen::Lower>();
  rv_cl.view(stan::math::matrix_cl_view::Upper);
  m_cl.view(stan::math::matrix_cl_view::Lower);
  m0_cl = rv_cl * m_cl;
  m0_cl_res = stan::math::from_matrix_cl(m0_cl);
  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
}

TEST(MathMatrix, multiply_small) {
  using stan::math::multiply;
  auto m1 = stan::math::matrix_d::Random(3, 3).eval();
  auto m2 = stan::math::matrix_d::Random(3, 3).eval();
  stan::math::matrix_d m3_cl_res(3, 3);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, multiply_big) {
  using stan::math::multiply;
  int size = 512;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, lower_tri_rect_multiply_small) {
  auto m1 = stan::math::matrix_d::Random(3, 3).eval();
  auto m2 = stan::math::matrix_d::Random(3, 3).eval();
  stan::math::matrix_d m3_cl_res(3, 3);

  stan::math::matrix_cl<double> m11(m1, stan::math::matrix_cl_view::Lower);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1.triangularView<Eigen::Lower>() * m2).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, lower_tri_rect_multiply_big) {
  int size = 512;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1, stan::math::matrix_cl_view::Lower);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1.triangularView<Eigen::Lower>() * m2).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, lower_tri_rect_multiply_big_rect) {
  int size = 321;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size * 3).eval();
  stan::math::matrix_d m3_cl_res(size, size * 3);

  stan::math::matrix_cl<double> m11(m1, stan::math::matrix_cl_view::Lower);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1.triangularView<Eigen::Lower>() * m2).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, upper_tri_rect_multiply_small) {
  auto m1 = stan::math::matrix_d::Random(3, 3).eval();
  auto m2 = stan::math::matrix_d::Random(3, 3).eval();
  stan::math::matrix_d m3_cl_res(3, 3);

  stan::math::matrix_cl<double> m11(m1, stan::math::matrix_cl_view::Upper);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1.triangularView<Eigen::Upper>() * m2).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, upper_tri_rect_multiply_big) {
  int size = 472;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1, stan::math::matrix_cl_view::Upper);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1.triangularView<Eigen::Upper>() * m2).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, upper_tri_rect_multiply_big_rect) {
  int size = 463;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size * 3).eval();
  stan::math::matrix_d m3_cl_res(size, size * 3);

  stan::math::matrix_cl<double> m11(m1, stan::math::matrix_cl_view::Upper);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1.triangularView<Eigen::Upper>() * m2).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, rect_lower_tri_multiply_small) {
  auto m1 = stan::math::matrix_d::Random(3, 3).eval();
  auto m2 = stan::math::matrix_d::Random(3, 3).eval();
  stan::math::matrix_d m3_cl_res(3, 3);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2, stan::math::matrix_cl_view::Lower);

  auto m3 = (m1 * m2.triangularView<Eigen::Lower>()).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, rect_lower_tri_multiply_big) {
  int size = 451;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2, stan::math::matrix_cl_view::Lower);

  auto m3 = (m1 * m2.triangularView<Eigen::Lower>()).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, rect_lower_tri_multiply_big_rect) {
  int size = 444;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size * 3).eval();
  stan::math::matrix_d m3_cl_res(size, size * 3);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2, stan::math::matrix_cl_view::Lower);

  auto m3 = (m1 * m2.triangularView<Eigen::Lower>()).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, rect_upper_tri_multiply_small) {
  auto m1 = stan::math::matrix_d::Random(3, 3).eval();
  auto m2 = stan::math::matrix_d::Random(3, 3).eval();
  stan::math::matrix_d m3_cl_res(3, 3);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2, stan::math::matrix_cl_view::Upper);

  auto m3 = (m1 * m2.triangularView<Eigen::Upper>()).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, rect_upper_tri_multiply_big) {
  int size = 468;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2, stan::math::matrix_cl_view::Upper);

  auto m3 = (m1 * m2.triangularView<Eigen::Upper>()).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, rect_upper_tri_multiply_big_rect) {
  int size = 345;
  auto m1 = stan::math::matrix_d::Random(size * 3, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size).eval();
  stan::math::matrix_d m3_cl_res(size * 3, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2, stan::math::matrix_cl_view::Upper);

  auto m3 = (m1 * m2.triangularView<Eigen::Upper>()).eval();

  auto m33 = m11 * m22;

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, multiply_big_split_4) {
  using stan::math::multiply;
  int size = 512;
  auto m1 = stan::math::matrix_d::Random(size, size * 2).eval();
  auto m2 = stan::math::matrix_d::Random(size * 2, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = stan::math::opencl::multiply(m11, m22);

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, multiply_big_split_11) {
  using stan::math::multiply;
  int size = 433;
  auto m1 = stan::math::matrix_d::Random(size, size * 11).eval();
  auto m2 = stan::math::matrix_d::Random(size * 11, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = stan::math::opencl::multiply(m11, m22);

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrix, multiply_small_split_big) {
  using stan::math::multiply;
  int size = 32;
  auto m1 = stan::math::matrix_d::Random(size, size * 10).eval();
  auto m2 = stan::math::matrix_d::Random(size * 10, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = stan::math::opencl::multiply(m11, m22);

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}
#endif
