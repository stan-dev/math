#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/util.hpp>
#include <stan/math/prim/mat/fun/matrix_exp.hpp>
#include <stan/math/prim/mat/fun/scale_matrix_exp_multiply.hpp>
#include <vector>

template <int N, int M>
inline void test_scale_matrix_exp_multiply() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::scale_matrix_exp_multiply;

  std::srand(1999);

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(N, N);
  Eigen::Matrix<double, -1, M> B = Eigen::MatrixXd::Random(N, M);
  Eigen::Matrix<double, Dynamic, Dynamic> A0 = A;

  // brute force
  Eigen::Matrix<double, N, M> expAB
      = stan::math::multiply(stan::math::matrix_exp(A0), B);

  // matrix_exp_multiply
  const double t = 1.0;
  Eigen::Matrix<double, N, M> res = scale_matrix_exp_multiply(t, A, B);
  EXPECT_EQ(res.size(), expAB.size());
  for (int l = 0; l < res.size(); ++l) {
    EXPECT_FLOAT_EQ(res(l), expAB(l));
  }
}

TEST(MathMatrix, scale_matrix_exp_multiply) {
  test_scale_matrix_exp_multiply<1, 1>();
  test_scale_matrix_exp_multiply<1, 5>();
  test_scale_matrix_exp_multiply<5, 1>();
  test_scale_matrix_exp_multiply<5, 5>();
  test_scale_matrix_exp_multiply<20, 2>();
}

TEST(MathMatrix, scale_matrix_exp_multiply_exception) {
  using stan::math::scale_matrix_exp_multiply;
  const double t = 1.0;
  {  // nonzero size
    Eigen::MatrixXd A(0, 0);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(1, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
    EXPECT_THROW(scale_matrix_exp_multiply(t, B, A), std::invalid_argument);
  }

  {  // multiplicable
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(2, 2);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
  }

  {  // square
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(2, 3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
  }
}
