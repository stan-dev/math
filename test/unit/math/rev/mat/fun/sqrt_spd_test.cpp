#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::matrix_v;
  using stan::math::sqrt_spd;

  matrix_v a(3, 3);
  a << 1.0, 2.0, 3.0, 2.0, 5.0, 7.9, 3.0, 7.9, 1.08;
  test::check_varis_on_stack(sqrt_spd(a));
}

struct make_zero {
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

    matrix_t sqrt_A = stan::math::sqrt_spd(A);
    matrix_t zero = A - sqrt_A * sqrt_A;
    zero.resize(K * K, 1);
    return zero;
  }
};

TEST(AgradRevMatrix, sqrt_spd) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::sqrt_spd;
  using stan::math::vector_v;

  int K = 2;
  matrix_d A(K, K);
  A.setRandom();
  A = A.transpose() * A;
  matrix_d J;
  Eigen::VectorXd f_x;
  A.resize(K * K, 1);
  make_zero f;
  stan::math::jacobian(f, A, f_x, J);
  const double TOL = 1e-14;
  int pos = 0;
  for (int i = 0; i < K; i++)
    for (int j = 0; j < K; j++)
      EXPECT_NEAR(f_x(pos++), 0, TOL);
  for (int i = 0; i < J.rows(); i++)
    for (int j = 0; j < J.cols(); j++)
      if (i != j)
        EXPECT_NEAR(J(i, j), 0, TOL);
  EXPECT_TRUE(J.array().any());
}
