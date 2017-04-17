#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/util.hpp>

using Eigen::Dynamic;
using Eigen::Matrix;
using std::vector;

TEST(ProbDistributionsMultiNormalPrec,MultiNormalVar) {
  using stan::math::var;
  Matrix<var,Dynamic,1> y(3,1);
  y << 2.0, -2.0, 11.0;
  Matrix<var,Dynamic,1> mu(3,1);
  mu << 1.0, -1.0, 3.0;
  Matrix<var,Dynamic,Dynamic> Sigma(3,3);
  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;
  Matrix<var,Dynamic,Dynamic> L = Sigma.inverse();
  EXPECT_FLOAT_EQ(-11.73908, stan::math::multi_normal_prec_log(y,mu,L).val());
}

TEST(AgradRev, check_varis_on_stack) {
  using stan::math::to_var;
  Matrix<double,Dynamic,1> y(3,1);
  y << 2.0, -2.0, 11.0;
  Matrix<double,Dynamic,1> mu(3,1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double,Dynamic,Dynamic> Sigma(3,3);
  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;
  Matrix<double,Dynamic,Dynamic> L = Sigma.inverse();
  
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<true>(to_var(y), to_var(mu), to_var(L)));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<true>(to_var(y), to_var(mu), L));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<true>(to_var(y), mu, to_var(L)));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<true>(to_var(y), mu, L));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<true>(y, to_var(mu), to_var(L)));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<true>(y, to_var(mu), L));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<true>(y, mu, to_var(L)));

  test::check_varis_on_stack(stan::math::multi_normal_prec_log<false>(to_var(y), to_var(mu), to_var(L)));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<false>(to_var(y), to_var(mu), L));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<false>(to_var(y), mu, to_var(L)));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<false>(to_var(y), mu, L));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<false>(y, to_var(mu), to_var(L)));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<false>(y, to_var(mu), L));
  test::check_varis_on_stack(stan::math::multi_normal_prec_log<false>(y, mu, to_var(L)));
}
