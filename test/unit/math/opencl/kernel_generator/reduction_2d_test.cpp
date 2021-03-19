#ifdef STAN_OPENCL

#include <stan/math.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <string>

using Eigen::Dynamic;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using stan::math::matrix_cl;

TEST(KernelGenerator, sum_2d_test) {
  std::string kernel_filename = "sum_2d.cl";
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  auto tmp = stan::math::sum_2d(m_cl);

  matrix_cl<double> res_cl;
  std::string kernel_src = tmp.get_kernel_source_for_evaluating_into(res_cl);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res_cl = tmp;

  MatrixXd raw_res = stan::math::from_matrix_cl(res_cl);
  EXPECT_GE(m.rows(), raw_res.rows());
  EXPECT_GE(m.cols(), raw_res.cols());
  double res = raw_res.sum();
  double correct = m.sum();
  EXPECT_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, prod_2d_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res_cl = stan::math::prod_2d(m_cl);
  MatrixXd raw_res = stan::math::from_matrix_cl(res_cl);
  EXPECT_GE(m.rows(), raw_res.rows());
  EXPECT_GE(m.cols(), raw_res.cols());
  double res = raw_res.prod();
  double correct = m.prod();
  EXPECT_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, min_2d_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res_cl = stan::math::min_2d(m_cl);
  MatrixXd raw_res = stan::math::from_matrix_cl(res_cl);
  EXPECT_GE(m.rows(), raw_res.rows());
  EXPECT_GE(m.cols(), raw_res.cols());
  double res = raw_res.minCoeff();
  double correct = m.minCoeff();
  EXPECT_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, max_2d_test) {
  MatrixXd m(2, 3);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res_cl = stan::math::max_2d(m_cl);
  MatrixXd raw_res = stan::math::from_matrix_cl(res_cl);
  EXPECT_GE(m.rows(), raw_res.rows());
  EXPECT_GE(m.cols(), raw_res.cols());
  double res = raw_res.maxCoeff();
  double correct = m.maxCoeff();
  EXPECT_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, sum_2d_triangular) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m, stan::math::matrix_cl_view::Lower);
  MatrixXd m1 = m;
  m1.triangularView<Eigen::StrictlyUpper>() = MatrixXd::Constant(3, 2, 0);

  matrix_cl<double> res1_cl = stan::math::sum_2d(m_cl);
  MatrixXd raw_res1 = stan::math::from_matrix_cl(res1_cl);
  EXPECT_GE(m1.rows(), raw_res1.rows());
  EXPECT_GE(m1.cols(), raw_res1.cols());
  double res1 = raw_res1.sum();
  double correct1 = m1.sum();
  EXPECT_NEAR(correct1, res1, 1e-9);

  m_cl.view(stan::math::matrix_cl_view::Upper);
  MatrixXd m2 = m;
  m2.triangularView<Eigen::StrictlyLower>() = MatrixXd::Constant(3, 2, 0);

  matrix_cl<double> res2_cl = stan::math::sum_2d(m_cl);
  MatrixXd raw_res2 = stan::math::from_matrix_cl(res2_cl);
  EXPECT_GE(m2.rows(), raw_res2.rows());
  EXPECT_GE(m2.cols(), raw_res2.cols());
  double res2 = raw_res2.sum();
  double correct2 = m2.sum();
  EXPECT_NEAR(correct2, res2, 1e-9);
}

TEST(KernelGenerator, nested_rowwise_sum_2d) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res_cl = stan::math::sum_2d(stan::math::rowwise_sum(m_cl));
  MatrixXd raw_res = stan::math::from_matrix_cl(res_cl);
  EXPECT_GE(m.rows(), raw_res.rows());
  EXPECT_GE(m.cols(), raw_res.cols());
  double res = raw_res.sum();
  double correct = m.sum();
  EXPECT_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, sum_2d_test_large) {
  for (int M : {1, 2, 5, 9, 63, 64, 65, 4095, 4096, 4967, 4096 * 4}) {
    for (int N : {1, 2, 5, 9, 63, 64, 65, 4095, 4096, 4967, 4096 * 4}) {
      if (N * M > 1e6) {
        continue;
      }
      MatrixXd m = MatrixXd::Random(N, M);

      matrix_cl<double> m_cl(m);
      matrix_cl<double> res_cl = stan::math::sum_2d(m_cl);
      MatrixXd raw_res = stan::math::from_matrix_cl(res_cl);
      EXPECT_GE(m.rows(), raw_res.rows());
      EXPECT_GE(m.cols(), raw_res.cols());
      double res = raw_res.sum();
      double correct = m.sum();
      EXPECT_NEAR(correct, res, 1e-9);
    }
  }
}

TEST(KernelGenerator, sum_2d_and_id_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl, res2_cl;
  stan::math::results(res1_cl, res2_cl)
      = stan::math::expressions(m_cl, stan::math::sum_2d(m_cl));
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  EXPECT_MATRIX_NEAR(m, res1, 1e-9);
  MatrixXd raw_res2 = stan::math::from_matrix_cl(res2_cl);
  EXPECT_GE(m.rows(), raw_res2.rows());
  EXPECT_GE(m.cols(), raw_res2.cols());
  double res2 = raw_res2.sum();
  double correct2 = m.sum();
  EXPECT_NEAR(correct2, res2, 1e-9);
}

TEST(KernelGenerator, reduction_2d_of_rowwise_broadcast_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;
  VectorXd v(3);
  v << 4, 5, 6;

  matrix_cl<double> m_cl(m);
  matrix_cl<double> v_cl(v);

  matrix_cl<double> res1_cl;
  matrix_cl<double> res2_cl;

  stan::math::results(res1_cl, res2_cl) = stan::math::expressions(
      m_cl, stan::math::sum_2d(stan::math::rowwise_broadcast(v_cl)));
  MatrixXd raw_res = stan::math::from_matrix_cl(res2_cl);
  EXPECT_GE(v.rows(), raw_res.rows());
  EXPECT_GE(m.cols(), raw_res.cols());
  double res = raw_res.sum();
  double correct = v.replicate(1, 2).sum();
  EXPECT_NEAR(correct, res, 1e-9);
}

#endif
