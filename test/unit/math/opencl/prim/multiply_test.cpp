#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto multiply_scalar_matrix_functor
    = [](const auto& a) { return stan::math::multiply(2.0, a); };
auto multiply_scalar_matrix_operator_functor
    = [](const auto& a) { return 2.0 * a; };
auto multiply_matrix_scalar_functor
    = [](const auto& a) { return stan::math::multiply(a, 2.0); };
auto multiply_matrix_scalar_operator_functor
    = [](const auto& a) { return a * 2.0; };
auto multiply_matrix_matrix_functor
    = [](const auto& a, const auto& b) { return stan::math::multiply(a, b); };
auto multiply_matrix_matrix_operator_functor
    = [](const auto& a, const auto& b) { return a * b; };

TEST(MathMatrixOpenCLPrim, one_dim_zero_matrix) {
  stan::math::matrix_d m0(5, 0);
  stan::math::matrix_d m1(0, 3);

  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_functor, m0);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_operator_functor, m0);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_functor, m0);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_operator_functor, m0);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_functor, m1);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_operator_functor, m1);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_functor, m1);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_operator_functor, m1);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_functor, m0, m1);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_operator_functor, m0, m1);
}

TEST(MathMatrixOpenCLPrim, zero_result_matrix) {
  stan::math::matrix_d m0(0, 5);
  stan::math::matrix_d m1(5, 0);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_functor, m0, m1);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_operator_functor, m0, m1);
}

TEST(MathMatrixOpenCLPrim, zero_size_input_matrix) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d m1(0, 0);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_functor, m0);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_operator_functor, m0);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_functor, m0);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_operator_functor, m0);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_functor, m0, m1);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_operator_functor, m0, m1);
}

TEST(MathMatrixOpenCLPrim, non_matching_dim_excpetion) {
  stan::math::matrix_d m0(5, 3);
  stan::math::matrix_d m1(2, 6);

  stan::math::matrix_cl<double> m0_cl(m0);
  stan::math::matrix_cl<double> m1_cl(m1);
  EXPECT_THROW(m0_cl * m1_cl, std::invalid_argument);
  EXPECT_THROW(stan::math::multiply(m0_cl, m1_cl), std::invalid_argument);
}

TEST(MathMatrixOpenCLPrim, multiply_scalar) {
  auto v = stan::math::vector_d::Random(25).eval();
  auto rv = stan::math::row_vector_d::Random(25).eval();
  auto m = stan::math::matrix_d::Random(5, 5).eval();

  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_functor, v);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_operator_functor, v);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_functor, v);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_operator_functor, v);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_functor, rv);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_operator_functor, rv);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_functor, rv);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_operator_functor, rv);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_functor, m);
  stan::math::test::compare_cpu_opencl_prim(multiply_scalar_matrix_operator_functor, m);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_functor, m);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_scalar_operator_functor, m);
  
}

TEST(MathMatrixOpenCLPrim, row_vector_vector) {
  stan::math::matrix_d v = stan::math::matrix_d::Random(5,1);
  stan::math::matrix_d rv = stan::math::matrix_d::Random(1,5);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_functor, rv, v);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_operator_functor, rv, v);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_functor, v, rv);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_operator_functor, v, rv);
}

TEST(MathMatrixOpenCLPrim, matrix_vector) {
  auto m = stan::math::matrix_d::Random(4, 5).eval();
  auto v = stan::math::vector_d::Random(5).eval();
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_functor, m, v);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_operator_functor, m, v);
  
  auto m1 = stan::math::matrix_d::Random(400, 600).eval();
  auto v1 = stan::math::vector_d::Random(600).eval();
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_functor, m1, v1);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_operator_functor, m1, v1);
}

TEST(MathMatrixOpenCLPrim, row_vector_matrix_small) {
  auto m = stan::math::matrix_d::Random(5, 4).eval();
  auto rv = stan::math::row_vector_d::Random(5).eval();
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_functor, rv, m);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_operator_functor, rv, m);
  auto m1 = stan::math::matrix_d::Random(600, 400).eval();
  auto rv1 = stan::math::row_vector_d::Random(600).eval();
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_functor, rv1, m1);
  stan::math::test::compare_cpu_opencl_prim(multiply_matrix_matrix_operator_functor, rv1, m1);
}

TEST(MathMatrixOpenCLPrim, matrix_vector_tri_small) {
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

TEST(MathMatrixOpenCLPrim, matrix_vector_tri_big) {
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

TEST(MathMatrixOpenCLPrim, row_vector_matrix_tri_small) {
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

TEST(MathMatrixOpenCLPrim, row_vector_matrix_tri_big) {
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

TEST(MathMatrixOpenCLPrim, multiply_small) {
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

TEST(MathMatrixOpenCLPrim, multiply_big) {
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

TEST(MathMatrixOpenCLPrim, lower_tri_rect_multiply_small) {
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

TEST(MathMatrixOpenCLPrim, lower_tri_rect_multiply_big) {
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

TEST(MathMatrixOpenCLPrim, lower_tri_rect_multiply_big_rect) {
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

TEST(MathMatrixOpenCLPrim, upper_tri_rect_multiply_small) {
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

TEST(MathMatrixOpenCLPrim, upper_tri_rect_multiply_big) {
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

TEST(MathMatrixOpenCLPrim, upper_tri_rect_multiply_big_rect) {
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

TEST(MathMatrixOpenCLPrim, rect_lower_tri_multiply_small) {
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

TEST(MathMatrixOpenCLPrim, rect_lower_tri_multiply_big) {
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

TEST(MathMatrixOpenCLPrim, rect_lower_tri_multiply_big_rect) {
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

TEST(MathMatrixOpenCLPrim, rect_upper_tri_multiply_small) {
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

TEST(MathMatrixOpenCLPrim, rect_upper_tri_multiply_big) {
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

TEST(MathMatrixOpenCLPrim, rect_upper_tri_multiply_big_rect) {
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

TEST(MathMatrixOpenCLPrim, multiply_big_split_4) {
  using stan::math::multiply;
  int size = 512;
  auto m1 = stan::math::matrix_d::Random(size, size * 2).eval();
  auto m2 = stan::math::matrix_d::Random(size * 2, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = stan::math::multiply(m11, m22);

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrixOpenCLPrim, multiply_big_split_11) {
  using stan::math::multiply;
  int size = 433;
  auto m1 = stan::math::matrix_d::Random(size, size * 11).eval();
  auto m2 = stan::math::matrix_d::Random(size * 11, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = stan::math::multiply(m11, m22);

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(MathMatrixOpenCLPrim, multiply_small_split_big) {
  using stan::math::multiply;
  int size = 32;
  auto m1 = stan::math::matrix_d::Random(size, size * 10).eval();
  auto m2 = stan::math::matrix_d::Random(size * 10, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl<double> m11(m1);
  stan::math::matrix_cl<double> m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = stan::math::multiply(m11, m22);

  m3_cl_res = stan::math::from_matrix_cl(m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}
#endif
