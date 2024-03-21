#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

TEST(ProbDistributionsMultiNormal, MultiNormalVar) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  Matrix<var, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<var, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<var, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  EXPECT_FLOAT_EQ(-11.73908, stan::math::multi_normal_lpdf(y, mu, Sigma).val());
}

TEST(ProbDistributionsMultiNormal, check_varis_on_stack) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::to_var;
  using std::vector;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  test::check_varis_on_stack(stan::math::multi_normal_lpdf<true>(
      to_var(y), to_var(mu), to_var(Sigma)));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<true>(to_var(y), to_var(mu), Sigma));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<true>(to_var(y), mu, to_var(Sigma)));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<true>(to_var(y), mu, Sigma));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<true>(y, to_var(mu), to_var(Sigma)));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<true>(y, to_var(mu), Sigma));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<true>(y, mu, to_var(Sigma)));

  test::check_varis_on_stack(stan::math::multi_normal_lpdf<false>(
      to_var(y), to_var(mu), to_var(Sigma)));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<false>(to_var(y), to_var(mu), Sigma));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<false>(to_var(y), mu, to_var(Sigma)));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<false>(to_var(y), mu, Sigma));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<false>(y, to_var(mu), to_var(Sigma)));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<false>(y, to_var(mu), Sigma));
  test::check_varis_on_stack(
      stan::math::multi_normal_lpdf<false>(y, mu, to_var(Sigma)));
}
