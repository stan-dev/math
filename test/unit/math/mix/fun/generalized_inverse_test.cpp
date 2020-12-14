#include <test/unit/math/test_ad.hpp>
#include <vector>
#include <gtest/gtest.h>

TEST(mathMixMatFun, ad_tests) {
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto f = [](const auto& G) { return stan::math::generalized_inverse(G); };
  auto fjit
      = [](const auto& G) { return stan::math::generalized_inverse(G, 0.0); };
  auto fjit2
      = [](const auto& G) { return stan::math::generalized_inverse(G, 1e-4); };

  /*
    Eigen::MatrixXd t(0, 0);
    expect_ad(f, t);
    expect_ad(fjit, t);
    expect_ad_matvar(f, t);

  Eigen::MatrixXd u(1, 1);
  u << 2;
  expect_ad(f, u);
  expect_ad(fjit, u);
  expect_ad_matvar(f, u);
*/
  Eigen::MatrixXd v(2, 3);
  v << 1, 3, 5, 2, 4, 6;
  expect_ad(f, v);
  /*
    v << 1.9, 1.3, 2.5, 0.4, 1.7, 0.1;
    expect_ad(f, v);
    expect_ad(fjit, v);
    expect_ad_matvar(f, v);

    // issues around zero require looser tolerances for hessians
    stan::test::ad_tolerances tols;
    tols.hessian_hessian_ = 2.0;
    tols.hessian_fvar_hessian_ = 2.0;

    Eigen::MatrixXd w(3, 4);
    w << 2, 3, 5, 7, 11, 13, 17, 19, 23, 25, 27, 29;
    expect_ad(tols, f, w);
    expect_ad(tols, fjit, w);
    expect_ad_matvar(f, w);

    Eigen::MatrixXd z(2, 2);
    z << 1, 2, 5, std::numeric_limits<double>::quiet_NaN();
    EXPECT_NO_THROW(stan::math::generalized_inverse(z));
    EXPECT_NO_THROW(stan::math::generalized_inverse(z, 0.0));

    // autodiff throws, so following fails (throw behavior must match to pass)

    Eigen::MatrixXd a(2, 2);
    a << 1.9, 0.3, 0.3, std::numeric_limits<double>::infinity();
    expect_ad(f, a);
    expect_ad(fjit, a);
    expect_ad_matvar(f, a);

    // singular matrix, should use the
    // alias to input small amount of jitter on the diagonal
    Eigen::MatrixXd m(3, 2);
    m << 1, 2, 2, 4, 1, 2;
    EXPECT_THROW(stan::math::generalized_inverse(m), std::domain_error);

    // should work with jittered version
    stan::test::ad_tolerances tols3;
    tols3.hessian_hessian_ = 0.01;
    tols3.hessian_fvar_hessian_ = 0.01;
    expect_ad(tols3, fjit2, m);
    */
}
