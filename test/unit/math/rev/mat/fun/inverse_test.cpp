#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>

TEST(AgradRevMatrix, inverse_val) {
  using stan::math::inverse;
  using stan::math::matrix_v;

  matrix_v a(2, 2);
  a << 2.0, 3.0, 5.0, 7.0;

  matrix_v a_inv = inverse(a);

  matrix_v I = multiply(a, a_inv);

  EXPECT_NEAR(1.0, I(0, 0).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1).val(), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0).val(), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1).val(), 1.0e-12);

  EXPECT_THROW(inverse(matrix_v(2, 3)), std::invalid_argument);
}
TEST(AgradRevMatrix, inverse_grad) {
  using stan::math::inverse;
  using stan::math::matrix_v;

  for (size_t k = 0; k < 2; ++k) {
    for (size_t l = 0; l < 2; ++l) {
      matrix_v ad(2, 2);
      ad << 2.0, 3.0, 5.0, 7.0;

      AVEC x = createAVEC(ad(0, 0), ad(0, 1), ad(1, 0), ad(1, 1));

      matrix_v ad_inv = inverse(ad);

      // int k = 0;
      // int l = 1;
      VEC g;
      ad_inv(k, l).grad(x, g);

      int idx = 0;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          EXPECT_FLOAT_EQ(-ad_inv(k, i).val() * ad_inv(j, l).val(), g[idx]);
          ++idx;
        }
      }
    }
  }
}
TEST(AgradRevMatrix, inverse_inverse_sum) {
  using stan::math::inverse;
  using stan::math::matrix_v;
  using stan::math::sum;

  matrix_v a(4, 4);
  a << 2.0, 3.0, 4.0, 5.0, 9.0, -1.0, 2.0, 2.0, 4.0, 3.0, 7.0, -1.0, 0.0, 1.0,
      19.0, 112.0;

  AVEC x;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      x.push_back(a(i, j));

  AVAR a_inv_inv_sum = sum(inverse(inverse(a)));

  VEC g = cgradvec(a_inv_inv_sum, x);

  for (size_t k = 0; k < x.size(); ++k)
    EXPECT_FLOAT_EQ(1.0, g[k]);
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::inverse;
  using stan::math::matrix_v;

  matrix_v a(2, 2);
  a << 2.0, 3.0, 5.0, 7.0;

  test::check_varis_on_stack(stan::math::inverse(a));
}

TEST(AgradRevMatrix, infinite_input_inverse) {
  using stan::math::inverse;
  using stan::math::matrix_v;

  matrix_v a(2, 2);
  a << 2.0, 3.0, 5.0, std::numeric_limits<double>::infinity();

  // inf input is allowed
  EXPECT_NO_THROW(inverse(a));

  // a_inv contains nan
  matrix_v a_inv = inverse(a);
  EXPECT_THROW(stan::math::multiply(a, a_inv), std::domain_error);
}

TEST(AgradRevMatrix, nan_input_inverse) {
  using stan::math::inverse;
  using stan::math::matrix_v;

  matrix_v a(2, 2);
  a << 2.0, 3.0, 5.0, std::numeric_limits<double>::quiet_NaN();

  // nan input is allowed
  EXPECT_NO_THROW(inverse(a));

  // a_inv contains nan
  matrix_v a_inv = inverse(a);
  EXPECT_THROW(stan::math::multiply(a, a_inv), std::domain_error);
}

TEST(AgradRevMatrix, exception_mat_inverse) {
  using stan::math::inverse;
  using stan::math::matrix_v;

  matrix_v a;

  // singular matrix, inverse does not exist.
  // Eigen does not throw errors, but returns inf
  a.resize(2, 2);
  a << 6.0, 4.0, 3.0, 2.0;
  EXPECT_NO_THROW(stan::math::inverse(a));

  // zero size
  a.resize(0, 0);
  EXPECT_THROW(stan::math::inverse(a), std::invalid_argument);

  // not square
  a.resize(2, 3);
  EXPECT_THROW(stan::math::inverse(a), std::invalid_argument);
}

TEST(AgradRevMatrix, inverse_1x1) {
  using stan::math::inverse;
  using stan::math::matrix_v;

  matrix_v a;
  a.resize(1, 1);
  a << 2.0;

  matrix_v a_inv = inverse(a);
  EXPECT_NEAR(0.5, a_inv(0, 0).val(), 1.0E-12);
}

TEST(AgradRevMatrix, inverse_ad) {
  Eigen::MatrixXd x(2, 2);
  x << 1.9, 0.3, 0.3, 1.7;

  auto g = [](const auto& u) { return stan::math::inverse(u); };
  stan::test::expect_ad(g, x);

  x << 1.9, 0.3, 0.3, std::numeric_limits<double>::infinity();
  stan::test::expect_ad(g, x);

  x << 1.9, 0.3, 0.3, std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::test::expect_ad(g, x), std::domain_error);

  x.resize(1, 1);
  x << 2.0;
  stan::test::expect_ad(g, x);
}
