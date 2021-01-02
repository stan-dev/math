#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixOpenCL, cholesky_decompose_cpu_vs_cl_small) {
  stan::math::matrix_d m0(3, 3);
  m0 << 25, 15, -5, 15, 18, 0, -5, 0, 11;

  stan::math::matrix_d m1(4, 4);
  m1 << 18, 22, 54, 42, 22, 70, 86, 62, 54, 86, 174, 134, 42, 62, 134, 106;

  stan::math::matrix_cl<double> m0_cl(m0);
  stan::math::matrix_cl<double> m1_cl(m1);

  stan::math::matrix_d m0_res = stan::math::cholesky_decompose(m0);
  stan::math::matrix_d m1_res = stan::math::cholesky_decompose(m1);

  stan::math::opencl::cholesky_decompose(m0_cl);
  stan::math::opencl::cholesky_decompose(m1_cl);

  m0 = stan::math::from_matrix_cl(m0_cl);
  m1 = stan::math::from_matrix_cl(m1_cl);

  EXPECT_MATRIX_NEAR(m0, m0_res, 1e-8);
  EXPECT_MATRIX_NEAR(m1, m1_res, 1e-8);
}

namespace {
void cholesky_decompose_test(int size) {
  stan::math::matrix_d m1 = stan::math::matrix_d::Random(size, size);
  stan::math::matrix_d m1_pos_def
      = m1 * m1.transpose() + size * Eigen::MatrixXd::Identity(size, size);

  stan::math::matrix_d m1_cpu(size, size);
  stan::math::matrix_d m1_cl(size, size);

  stan::math::check_symmetric("cholesky_decompose", "m", m1_pos_def);
  Eigen::LLT<stan::math::matrix_d> llt(m1_pos_def.rows());
  llt.compute(m1_pos_def);
  stan::math::check_pos_definite("cholesky_decompose", "m", llt);
  m1_cpu = llt.matrixL();

  m1_cl = stan::math::cholesky_decompose(m1_pos_def);

  double max_error = 0;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j <= i; j++) {
      double abs_err = std::fabs(m1_cpu(i, j) - m1_cl(i, j));
      double a = std::max(abs_err / m1_cpu(i, j), abs_err / m1_cl(i, j));
      max_error = std::max(max_error, a);
    }
  }
  EXPECT_LT(max_error, 1e-8);
}
}  // namespace

TEST(MathMatrixOpenCL, cholesky_decompose_small) {
  cholesky_decompose_test(10);
  cholesky_decompose_test(50);
  cholesky_decompose_test(100);
}

TEST(MathMatrixOpenCL, cholesky_decompose_big) {
  cholesky_decompose_test(1251);
  cholesky_decompose_test(1704);
  cholesky_decompose_test(2000);
}

TEST(MathMatrixOpenCL, cholesky_decompose_big_tuning_opts) {
  std::vector<int> size_transfer({256, 512, 1300});
  std::vector<int> cholesky_min_size({64, 256});
  std::vector<int> cholesky_part({2, 4});
  for (auto&& size_t_ : size_transfer) {
    for (auto&& min_size_ : cholesky_min_size) {
      for (auto&& part_ : cholesky_part) {
        stan::math::opencl_context.tuning_opts().cholesky_size_worth_transfer
            = size_t_;
        stan::math::opencl_context.tuning_opts().cholesky_min_L11_size
            = min_size_;
        stan::math::opencl_context.tuning_opts().cholesky_partition = part_;
        cholesky_decompose_test(1300);
      }
    }
    stan::math::opencl_context.tuning_opts().cholesky_size_worth_transfer = 128;
    stan::math::opencl_context.tuning_opts().cholesky_min_L11_size = 64;
    stan::math::opencl_context.tuning_opts().cholesky_partition = 4;
    cholesky_decompose_test(128);
    cholesky_decompose_test(65);
    cholesky_decompose_test(130);
    cholesky_decompose_test(128 * 4);
    cholesky_decompose_test(128 * 4 - 1);
  }
}
#endif
