#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::inverse_spd;
  using stan::math::matrix_v;

  matrix_v a(2, 2);
  a << 2.0, 3.0, 3.0, 7.0;

  test::check_varis_on_stack(stan::math::inverse_spd(a));
}

struct make_I {
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1> &a) const {
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
    const int K = static_cast<int>(sqrt(a.rows()));
    matrix_t A(K, K);
    int pos = 0;
    for (int i = 0; i < K; i++)
      for (int j = 0; j < K; j++)
        A(i, j) = a[pos++];

    matrix_t I = A * stan::math::inverse_spd(A);
    I.resize(K * K, 1);
    return I;
  }
};

TEST(AgradRevMatrix, inverse) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::sqrt_spd;
  using stan::math::vector_v;

  int K = 11;
  matrix_d A(K, K);
  A.setRandom();
  A = A.transpose() * A;
  matrix_d J;
  Eigen::VectorXd f_x;
  A.resize(K * K, 1);
  make_I f;
  stan::math::jacobian(f, A, f_x, J);
  const double TOL = 1e-13;
  int pos = 0;
  for (int i = 0; i < K; i++)
    for (int j = 0; j < K; j++)
      EXPECT_NEAR(f_x(pos++), i == j, TOL);
  for (int i = 0; i < J.rows(); i++)
    for (int j = 0; j < J.cols(); j++)
      EXPECT_NEAR(J(i, j), 0, 10 * TOL);
  EXPECT_TRUE(J.array().any());
}

TEST(AgradRevMatrix, inverse_spd_val) {
  using stan::math::inverse_spd;
  using stan::math::matrix_v;

  matrix_v a(2, 2);
  a << 2.0, 3.0, 3.0, 7.0;

  matrix_v a_inv = inverse_spd(a);

  matrix_v I = multiply(a, a_inv);

  EXPECT_NEAR(1.0, I(0, 0).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0).val(), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1).val(), 1.0e-12);

  EXPECT_THROW(inverse_spd(matrix_v(2, 3)), std::invalid_argument);

  a << 2.0, 3.0, 1.0, 7.0;
  EXPECT_THROW(inverse_spd(a), std::domain_error);
  a << 1.0, -1.0, -1.0, -1.0;
  EXPECT_THROW(inverse_spd(a), std::domain_error);
}
/*
TEST(AgradRevMatrix, inverse_spd_grad) {
  using stan::math::inverse_spd;
  using stan::math::matrix_v;

  for (size_t k = 0; k < 2; ++k) {
    for (size_t l = 0; l < 2; ++l) {
      matrix_v ad(2, 2);
      ad << 2.0, 3.0, 3.0, 7.0;

      AVEC x = createAVEC(ad(0, 0), ad(0, 1), ad(1, 0), ad(1, 1));

      matrix_v ad_inv = inverse_spd(ad);
      std::cout << "ad_inv = " << std::endl << ad_inv << std::endl;

      // int k = 0;
      // int l = 1;
      VEC g(4);
      (0.5 * (ad_inv(k, l) + ad_inv(l, k))).grad(x, g);
      std::cout << "passed" << std::endl;
      int idx = 0;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          EXPECT_FLOAT_EQ(-0.5 * ad_inv(k, i).val() * ad_inv(j, l).val()
                              - 0.5 * ad_inv(l, i).val() * ad_inv(j, k).val(),
                          g[idx]);
          ++idx;
        }
      }
    }
  }
}

TEST(AgradRevMatrix, inverse_spd_inverse_spd_sum) {
  using stan::math::inverse_spd;
  using stan::math::matrix_v;
  using stan::math::sum;

  matrix_v a(4, 4);
  a << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
      0.0, 1.0;

  AVEC x;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      x.push_back(a(i, j));

  AVAR a_inv_inv_sum = sum(inverse_spd(inverse_spd(a)));

  VEC g;
  a_inv_inv_sum.grad(x, g);

  size_t k = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      EXPECT_FLOAT_EQ(1.0, g[k]);
      k++;
    }
  }
}
*/
