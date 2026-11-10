#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <stan/math/mix.hpp>
#include <vector>

TEST(mathMixScalFun, gammaP) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::gamma_p(x1, x2);
  };
  // common nonzero double arguments excluding positive infinity
  // due to issues with finite diff gradient
  const std::vector<double> args{-1.3,
                                 0.49,
                                 0.99,
                                 1.01,
                                 stan::math::negative_infinity(),
                                 stan::math::not_a_number()};
  auto int_args = stan::test::internal::common_nonzero_int_args();
  for (double x1 : args)
    for (double x2 : args) {
      stan::test::expect_ad(f, x1, x2);
    }
  for (double x1 : args)
    for (int x2 : int_args) {
      stan::test::expect_ad(f, x1, x2);
    }
  for (int x1 : int_args)
    for (double x2 : args) {
      stan::test::expect_ad(f, x1, x2);
    }
  for (int x1 : int_args)
    for (int x2 : int_args) {
      stan::test::expect_ad(f, x1, x2);
    }
  for (double x2 : args) {
    stan::test::expect_ad(f, stan::math::positive_infinity(), x2);
  }
  for (int x2 : int_args) {
    stan::test::expect_ad(f, stan::math::positive_infinity(), x2);
  }
}

// separate tests when a is positive_infinity
TEST(mathMixScalFun, gammaP_pos_inf) {
  auto g = [](const auto& x) { return stan::math::gamma_p(x(0), x(1)); };
  stan::math::vector_d x(2);
  x << 0.5001, stan::math::positive_infinity();
  Eigen::Matrix<double, Eigen::Dynamic, 1> grad_ad;
  double fx_ad = 1;

  stan::math::gradient(g, x, fx_ad, grad_ad);
  stan::test::expect_near_rel("gradient() val", 1, fx_ad, 1e-8);
  stan::test::expect_near_rel("gradient() grad z", stan::math::not_a_number(),
                              grad_ad(0), 1e-3);
  stan::test::expect_near_rel("gradient() grad a", stan::math::not_a_number(),
                              grad_ad(1), 1e-3);

  Eigen::Matrix<double, Eigen::Dynamic, 1> grad_ad_fvar;
  stan::math::gradient<double>(g, x, fx_ad, grad_ad_fvar);
  stan::test::expect_near_rel("gradient_fvar() val", 1, fx_ad, 1e-8);
  stan::test::expect_near_rel("gradient_fvar() grad z",
                              stan::math::not_a_number(), grad_ad_fvar(0),
                              1e-3);
  stan::test::expect_near_rel("gradient_fvar() grad a",
                              stan::math::not_a_number(), grad_ad_fvar(1),
                              1e-3);

  Eigen::Matrix<double, Eigen::Dynamic, 1> hessian_grad_ad;
  Eigen::MatrixXd H_ad;
  stan::math::hessian(g, x, fx_ad, hessian_grad_ad, H_ad);
  stan::test::expect_near_rel("hessian() val", 1, fx_ad, 1e-8);
  stan::test::expect_near_rel("hessian() grad z", stan::math::not_a_number(),
                              hessian_grad_ad(0), 1e-3);
  stan::test::expect_near_rel("hessian() grad a", stan::math::not_a_number(),
                              hessian_grad_ad(1), 1e-3);
  stan::test::expect_near_rel("hessian() Hessian z", stan::math::not_a_number(),
                              H_ad(0), 1e-3);
  stan::test::expect_near_rel("hessian() Hessian a", stan::math::not_a_number(),
                              H_ad(1), 1e-3);

  std::vector<Eigen::MatrixXd> grad_H_ad;
  stan::math::grad_hessian(g, x, fx_ad, H_ad, grad_H_ad);
  stan::test::expect_near_rel("grad_hessian() val", 1, fx_ad, 1e-8);
  for (int i = 0; i < H_ad.size(); i++) {
    stan::test::expect_near_rel("grad_hessian() Hessian", 0, H_ad(i), 1e-3);
  }
  for (int i = 0; i < grad_H_ad.size(); i++) {
    for (int j = 0; j < grad_H_ad[i].size(); j++) {
      stan::test::expect_near_rel("gradient() grad Hessian",
                                  stan::math::not_a_number(), grad_H_ad[i](i),
                                  1e-3);
    }
  }
}

TEST(mathMixScalFun, gammaP_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::gamma_p;
    return gamma_p(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}

// Reverse-mode derivatives far in the upper tail, z / a > 10, where both
// derivatives are about 1e-13. References from mpmath at 50 digits.
TEST(mathMixScalFun, gammaP_rev_far_upper_tail) {
  using stan::math::gamma_p;
  using stan::math::var;
  const double d_a = -4.2728631691092249555e-13;
  const double d_z = 1.5030507939978329636e-13;
  {
    var a = 3.0;
    var z = 36.0;
    var p = gamma_p(a, z);
    p.grad();
    stan::test::expect_near_rel("vv d/da", d_a, a.adj(), 1e-10);
    stan::test::expect_near_rel("vv d/dz", d_z, z.adj(), 1e-10);
    stan::math::recover_memory();
  }
  {
    var a = 3.0;
    var p = gamma_p(a, 36.0);
    p.grad();
    stan::test::expect_near_rel("vd d/da", d_a, a.adj(), 1e-10);
    stan::math::recover_memory();
  }
  {
    var z = 36.0;
    var p = gamma_p(3.0, z);
    p.grad();
    stan::test::expect_near_rel("dv d/dz", d_z, z.adj(), 1e-10);
    stan::math::recover_memory();
  }
  // log1m(gamma_p(alpha, beta * y)) = log Q at y = 12, alpha = 3, beta = 3.
  // 1 - P is about 1.6e-13 and keeps only about three digits in double.
  {
    var y = 12.0;
    var alpha = 3.0;
    var beta = 3.0;
    var lp = stan::math::log1m(gamma_p(alpha, beta * y));
    lp.grad();
    stan::test::expect_near_rel("log1m d/dy", -2.8379562043795620438, y.adj(),
                                1e-2);
    stan::test::expect_near_rel("log1m d/dalpha", 2.6892412395869051037,
                                alpha.adj(), 1e-2);
    stan::test::expect_near_rel("log1m d/dbeta", -11.351824817518248175,
                                beta.adj(), 1e-2);
    stan::math::recover_memory();
  }
}
