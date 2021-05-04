#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(OpenCLPrim, add_exceptions) {
  stan::math::vector_d vd1(2), vd2(3);
  stan::math::matrix_cl<double> vd11(vd1);
  stan::math::matrix_cl<double> vd22(vd2);
  EXPECT_THROW(stan::math::add(vd11, vd22), std::invalid_argument);

  stan::math::row_vector_d rvd1(2), rvd2(3);
  stan::math::matrix_cl<double> rvd11(rvd1);
  stan::math::matrix_cl<double> rvd22(rvd2);
  EXPECT_THROW(stan::math::add(rvd11, rvd22), std::invalid_argument);

  stan::math::matrix_d md1(2, 2), md2(3, 3);
  stan::math::matrix_cl<double> md11(md1);
  stan::math::matrix_cl<double> md22(md2);
  EXPECT_THROW(stan::math::add(md11, md22), std::invalid_argument);
}

TEST(OpenCLPrim, add_tri_value_check) {
  Eigen::MatrixXd a(3, 3);
  a << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::MatrixXd b = Eigen::MatrixXd::Ones(3, 3) * -3;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);
  stan::math::matrix_cl<double> c_cl(3, 3);
  stan::math::matrix_cl<double> c_cl_fun(3, 3);
  Eigen::MatrixXd c(3, 3);

  a_cl.view(stan::math::matrix_cl_view::Lower);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  c_cl = a_cl + b_cl;
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Lower);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR((Eigen::MatrixXd(a.triangularView<Eigen::Lower>())
                      + Eigen::MatrixXd(b.triangularView<Eigen::Lower>())),
                     c, 1E-8);

  c_cl_fun = add(a_cl, b_cl);
  EXPECT_EQ(c_cl_fun.view(), stan::math::matrix_cl_view::Lower);
  c = stan::math::from_matrix_cl(c_cl_fun);
  EXPECT_MATRIX_NEAR((Eigen::MatrixXd(a.triangularView<Eigen::Lower>())
                      + Eigen::MatrixXd(b.triangularView<Eigen::Lower>())),
                     c, 1E-8);

  a_cl.view(stan::math::matrix_cl_view::Lower);
  b_cl.view(stan::math::matrix_cl_view::Upper);
  c_cl = a_cl + b_cl;
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR((Eigen::MatrixXd(a.triangularView<Eigen::Lower>())
                      + Eigen::MatrixXd(b.triangularView<Eigen::Upper>())),
                     c, 1E-8);

  c_cl_fun = add(a_cl, b_cl);
  EXPECT_EQ(c_cl_fun.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl_fun);
  EXPECT_MATRIX_NEAR((Eigen::MatrixXd(a.triangularView<Eigen::Lower>())
                      + Eigen::MatrixXd(b.triangularView<Eigen::Upper>())),
                     c, 1E-8);

  a_cl.view(stan::math::matrix_cl_view::Upper);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  c_cl = a_cl + b_cl;
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR((Eigen::MatrixXd(a.triangularView<Eigen::Upper>())
                      + Eigen::MatrixXd(b.triangularView<Eigen::Lower>())),
                     c, 1E-8);

  c_cl_fun = add(a_cl, b_cl);
  EXPECT_EQ(c_cl_fun.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl_fun);
  EXPECT_MATRIX_NEAR((Eigen::MatrixXd(a.triangularView<Eigen::Upper>())
                      + Eigen::MatrixXd(b.triangularView<Eigen::Lower>())),
                     c, 1E-8);

  a_cl.view(stan::math::matrix_cl_view::Entire);
  b_cl.view(stan::math::matrix_cl_view::Lower);
  c_cl = a_cl + b_cl;
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR((a + Eigen::MatrixXd(b.triangularView<Eigen::Lower>())), c,
                     1E-8);

  c_cl_fun = add(a_cl, b_cl);
  EXPECT_EQ(c_cl_fun.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl_fun);
  EXPECT_MATRIX_NEAR((a + Eigen::MatrixXd(b.triangularView<Eigen::Lower>())), c,
                     1E-8);
}

TEST(OpenCLPrim, add_tri_scalar_value_check) {
  Eigen::MatrixXd a(3, 3);
  a << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> c_cl(3, 3);
  Eigen::MatrixXd c(3, 3);

  using stan::math::add;
  a_cl.view(stan::math::matrix_cl_view::Lower);
  c_cl = add(a_cl, 1.5);
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR(
      add(Eigen::MatrixXd(a.triangularView<Eigen::Lower>()), 1.5), c, 1E-8);

  a_cl.view(stan::math::matrix_cl_view::Lower);
  c_cl = add(1.5, a_cl);
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR(
      add(1.5, Eigen::MatrixXd(a.triangularView<Eigen::Lower>())), c, 1E-8);

  a_cl.view(stan::math::matrix_cl_view::Upper);
  c_cl = add(a_cl, 1.5);
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR(
      add(Eigen::MatrixXd(a.triangularView<Eigen::Upper>()), 1.5), c, 1E-8);

  a_cl.view(stan::math::matrix_cl_view::Upper);
  c_cl = add(1.5, a_cl);
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR(
      add(1.5, Eigen::MatrixXd(a.triangularView<Eigen::Upper>())), c, 1E-8);

  a_cl.view(stan::math::matrix_cl_view::Entire);
  c_cl = add(a_cl, 1.5);
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR(add(a, 1.5), c, 1E-8);

  a_cl.view(stan::math::matrix_cl_view::Entire);
  c_cl = add(1.5, a_cl);
  EXPECT_EQ(c_cl.view(), stan::math::matrix_cl_view::Entire);
  c = stan::math::from_matrix_cl(c_cl);
  EXPECT_MATRIX_NEAR(add(1.5, a), c, 1E-8);
}

TEST(OpenCLPrim, add_batch) {
  // used to represent 5 matrices of size 10x10
  const int batch_size = 11;
  const int size = 13;
  stan::math::matrix_d a(size, size * batch_size);
  stan::math::matrix_d a_res(size, size);
  for (int k = 0; k < batch_size; k++) {
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) {
        a(i, k * size + j) = k;
      }
  }
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> a_cl_res(size, size);
  stan::math::opencl_kernels::add_batch(cl::NDRange(size, size), a_cl_res, a_cl,
                                        size, size, batch_size);
  a_res = stan::math::from_matrix_cl(a_cl_res);
  for (int k = 0; k < batch_size; k++) {
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) {
        a(i, j) += a(i, k * size + j);
      }
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      EXPECT_EQ(a(i, j), a_res(i, j));
    }
  }
}
#endif
