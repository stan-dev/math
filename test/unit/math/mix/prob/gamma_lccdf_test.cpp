#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/rev/prob/gamma_lccdf_tail_reference.hpp>
#include <cmath>
#include <string>

TEST_F(AgradRev, mathMixScalFun_gamma_lccdf) {
  auto f = [](const auto& y, const auto& alpha, const auto& beta) {
    return stan::math::gamma_lccdf(y, alpha, beta);
  };

  stan::test::expect_ad(f, 0.8, 1.1, 2.3);
  stan::test::expect_ad(f, 2.5, 3.0, 1.5);
  stan::test::expect_ad(f, 6.5, 5.0, 1.0);
}

TEST_F(AgradRev, mathMixScalFun_gamma_lccdf_upper_tail_reference) {
  using stan::math::fvar;
  using stan::math::gamma_lccdf;
  using stan::math::var;
  namespace ref = stan::test::gamma_lccdf_tail;

  for (const auto& r : ref::rows()) {
    const std::string at = ref::label(r);

    fvar<double> a1(r.alpha, 1.0);
    ref::expect_rel(r.d_alpha, gamma_lccdf(r.y, a1, r.beta).d_, 1e-11,
                    "fvar<double> d/dalpha " + at);

    var a2_val = r.alpha;
    fvar<var> a2(a2_val, 1.0);
    fvar<var> lp2 = gamma_lccdf(r.y, a2, r.beta);
    ref::expect_rel(r.d_alpha, lp2.d_.val(), 1e-11, "fvar<var> d/dalpha " + at);
    lp2.d_.grad();
    ref::expect_rel(r.d2_alpha, a2_val.adj(), ref::d2_tol(r),
                    "fvar<var> d2/dalpha2 " + at);
    stan::math::recover_memory();

    fvar<fvar<double>> a3(fvar<double>(r.alpha, 1.0), fvar<double>(1.0, 0.0));
    fvar<fvar<double>> lp3 = gamma_lccdf(r.y, a3, r.beta);
    ref::expect_rel(r.value, lp3.val_.val_, 1e-10,
                    "fvar<fvar<double>> value " + at);
    ref::expect_rel(r.d_alpha, lp3.val_.d_, 1e-11,
                    "fvar<fvar<double>> inner d/dalpha " + at);
    ref::expect_rel(r.d_alpha, lp3.d_.val_, 1e-11,
                    "fvar<fvar<double>> outer d/dalpha " + at);
    ref::expect_rel(r.d2_alpha, lp3.d_.d_, ref::d2_tol(r),
                    "fvar<fvar<double>> d2/dalpha2 " + at);

    var a4_val = r.alpha;
    fvar<fvar<var>> a4(fvar<var>(a4_val, 1.0), fvar<var>(1.0, 0.0));
    fvar<fvar<var>> lp4 = gamma_lccdf(r.y, a4, r.beta);
    ref::expect_rel(r.d_alpha, lp4.d_.val_.val(), 1e-11,
                    "fvar<fvar<var>> d/dalpha " + at);
    ref::expect_rel(r.d2_alpha, lp4.d_.d_.val(), ref::d2_tol(r),
                    "fvar<fvar<var>> d2/dalpha2 " + at);
    lp4.d_.d_.grad();
    EXPECT_TRUE(std::isfinite(a4_val.adj()))
        << "fvar<fvar<var>> d3/dalpha3 " << at;
    stan::math::recover_memory();
  }
}
