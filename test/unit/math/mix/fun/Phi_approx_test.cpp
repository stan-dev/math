#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/Phi_approx_tail_refs.hpp>
#include <vector>

TEST(mathMixMatFun, PhiApprox) {
  auto f = [](const auto& x1) { return stan::math::Phi_approx(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -3.0, 1, 1.3, 3);
}

TEST(mathMixMatFun, PhiApprox_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::Phi_approx;
    return Phi_approx(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-3.0, 1, 1.3, 3};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}

TEST(mathMixMatFun, PhiApproxScalarDerivativeTails) {
  using stan::math::var;
  std::vector<double> xs = Phi_approx_tail_refs::args();
  std::vector<double> ds = Phi_approx_tail_refs::derivs();
  for (size_t i = 0; i < xs.size(); ++i) {
    var a = xs[i];
    stan::math::Phi_approx(a).grad();
    Phi_approx_tail_refs::expect_rel(ds[i], a.adj(), xs[i]);
    stan::math::recover_memory();
  }
}

TEST(mathMixMatFun, PhiApproxDerivativeTails) {
  using Phi_approx_tail_refs::expect_rel;
  using stan::math::Phi_approx;
  using stan::math::sum;
  using stan::math::var;
  std::vector<double> xs = Phi_approx_tail_refs::args();
  std::vector<double> ds = Phi_approx_tail_refs::derivs();
  Eigen::VectorXd x = stan::math::to_vector(xs);

  stan::math::var_value<Eigen::VectorXd> xv(x);
  sum(Phi_approx(xv)).grad();
  for (int i = 0; i < x.size(); ++i) {
    expect_rel(ds[i], xv.adj()(i), xs[i]);
  }
  stan::math::recover_memory();

  Eigen::Matrix<var, Eigen::Dynamic, 1> xm = x;
  sum(Phi_approx(xm)).grad();
  for (int i = 0; i < x.size(); ++i) {
    expect_rel(ds[i], xm(i).adj(), xs[i]);
  }
  stan::math::recover_memory();

  std::vector<var> xs_v(xs.begin(), xs.end());
  sum(Phi_approx(xs_v)).grad();
  for (size_t i = 0; i < xs.size(); ++i) {
    expect_rel(ds[i], xs_v[i].adj(), xs[i]);
  }
  stan::math::recover_memory();
}

// Forward mode differentiates through the fvar overload of inv_logit, whose
// derivative was 0 wherever the reverse-mode one was.
TEST(mathMixMatFun, PhiApproxFwdDerivativeTails) {
  using Phi_approx_tail_refs::expect_rel;
  using stan::math::fvar;
  using stan::math::Phi_approx;
  using stan::math::var;
  std::vector<double> xs = Phi_approx_tail_refs::args();
  std::vector<double> ds = Phi_approx_tail_refs::derivs();
  std::vector<double> d2s = Phi_approx_tail_refs::second_derivs();
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
