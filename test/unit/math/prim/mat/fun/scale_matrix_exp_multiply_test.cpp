#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/util.hpp>
#include <stan/math/prim/mat/fun/matrix_exp.hpp>
#include <stan/math/prim/mat/fun/scale_matrix_exp_multiply.hpp>
#include <vector>

template <int N, int M>
inline void test_scale_matrix_exp_multiply() {
  using stan::math::scale_matrix_exp_multiply;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  std::srand(1999);

  Eigen::Matrix<double, N, N> A = Eigen::Matrix<double, N, N>::Random();
  Eigen::Matrix<double, N, M> B = Eigen::Matrix<double, N, M>::Random();
  Eigen::Matrix<double, Dynamic, Dynamic> A0 = A;

  // brute force
  Eigen::Matrix<double, N, N> expA = stan::math::matrix_exp(A0);
  Eigen::Matrix<double, N, M> expAB;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      expAB(i, j) = 0.0;
      for (int k = 0; k < N; ++k) {
        expAB(i, j) += expA(i, k) * B(k, j);
      }
    }
  }

  // matrix_exp_multiply
  const double t = 1.0;
  Eigen::Matrix<double, N, M> res_dv = scale_matrix_exp_multiply(t, A, B);
  for (int l = 0; l < res_dv.size(); ++l) {
    EXPECT_FLOAT_EQ(res_dv(l), expAB(l));
  }
}

TEST(MathMatrix, scale_matrix_exp_multiply) {
  test_scale_matrix_exp_multiply<1, 1>();
  test_scale_matrix_exp_multiply<1, 5>();
  test_scale_matrix_exp_multiply<5, 1>();
  test_scale_matrix_exp_multiply<5, 5>();
  test_scale_matrix_exp_multiply<10, 10>();
}
