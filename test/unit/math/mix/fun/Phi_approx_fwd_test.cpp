#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/inv_logit_tail_refs.hpp>
#include <vector>

// Forward-mode derivatives of Phi_approx(x) = inv_logit(0.07056 x^3 + 1.5976 x)
// come from the fvar overload of inv_logit; they were 0 wherever its own
// derivative was.
namespace Phi_approx_fwd_test {
inline std::vector<double> args() {
  return {-40, -10, -8, -7, -1, 0, 1, 7, 8, 10, 40};
}
// first derivatives by mpmath at 60 digits; at x = 40 and -40 the true
// value, about 4e-1987, underflows to 0
inline std::vector<double> derivs() {
  return {0.0,
          5.9589784505245725e-37,
          8.7097694481792379e-21,
          5.1340418352583332e-15,
          0.24152728992165639,
          0.3994,
          0.24152728992165639,
          5.1340418352583332e-15,
          8.7097694481792379e-21,
          5.9589784505245725e-37,
          0.0};
}
// second derivatives, computed the same way
inline std::vector<double> second_derivs() {
  return {0.0,
          1.3455155964317485e-35,
          1.2996275105706113e-19,
          6.018298087621098e-14,
          0.24180038682517271,
          0.0,
          -0.24180038682517271,
          -6.018298087621098e-14,
          -1.2996275105706113e-19,
          -1.3455155964317485e-35,
          0.0};
}
}  // namespace Phi_approx_fwd_test

TEST(mathMixMatFun, PhiApproxFwdDerivativeTails) {
  using inv_logit_tail_refs::expect_rel;
  using stan::math::fvar;
  using stan::math::Phi_approx;
  using stan::math::var;
  std::vector<double> xs = Phi_approx_fwd_test::args();
  std::vector<double> ds = Phi_approx_fwd_test::derivs();
  std::vector<double> d2s = Phi_approx_fwd_test::second_derivs();
  for (size_t i = 0; i < xs.size(); ++i) {
    expect_rel(ds[i], Phi_approx(fvar<double>(xs[i], 1)).d_, xs[i]);

    fvar<fvar<double>> ff(fvar<double>(xs[i], 1), fvar<double>(1, 0));
    fvar<fvar<double>> y_ff = Phi_approx(ff);
    expect_rel(ds[i], y_ff.d_.val_, xs[i]);
    expect_rel(d2s[i], y_ff.d_.d_, xs[i]);

    fvar<var> fv(xs[i], 1);
    fvar<var> y_fv = Phi_approx(fv);
    expect_rel(ds[i], y_fv.d_.val(), xs[i]);
    y_fv.d_.grad();
    expect_rel(d2s[i], fv.val_.adj(), xs[i]);
    stan::math::recover_memory();
  }
}
