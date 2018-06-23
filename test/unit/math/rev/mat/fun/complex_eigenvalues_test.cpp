#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

struct make_I {
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& a) const {
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
    const int K = static_cast<int>(sqrt(a.rows()));
    matrix_t A(K, K);
    int pos = 0;
    for (int i = 0; i < K; i++)
      for (int j = 0; j < K; j++)
        A(i, j) = a[pos++];
    Eigen::EigenSolver<matrix_t> es(A);
    matrix_t I = (es.eigenvectors().inverse() * A * es.eigenvectors()
                  * es.eigenvalues().asDiagonal().inverse())
                     .real();
    I.resize(K * K, 1);
    return I;
  }
};

TEST(MathMatrix, eigen_decomposition) {
  using stan::math::matrix_d;
  int K = 9;
  matrix_d A(K, K);
  A.setRandom();
  matrix_d J;
  Eigen::VectorXd f_x;
  A.resize(K * K, 1);
  make_I f;
  stan::math::jacobian(f, A, f_x, J);
  const double TOL = 7e-13;
  int pos = 0;
  for (int j = 0; j < K; j++) {
    for (int i = 0; i < j; i++)
      EXPECT_NEAR(f_x[pos++], 0, TOL);
    EXPECT_NEAR(f_x[pos++], 1, TOL);
    for (int i = j + 1; i < K; i++)
      EXPECT_NEAR(f_x[pos++], 0, TOL);
  }
  for (int i = 0; i < J.rows(); i++)
    for (int j = 0; j < J.cols(); j++)
      EXPECT_NEAR(J(i, j), 0, TOL);
  EXPECT_TRUE(J.array().any());
}
