#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/inv_logit_tail_refs.hpp>

TEST(mathMixMatFun, invLogit) {
  auto f = [](const auto& x1) { return stan::math::inv_logit(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1.2, -0.2, 0.5, 1, 1.3, 1.5,
                                      3);

  std::vector<double> com_args = stan::test::internal::common_nonzero_args();
  std::vector<double> args{-2.6, -0.5, 0.5, 1.5};

  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(com_args));
  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(args));
}

TEST(mathMixMatFun, invLogitDerivativesLargeArgs) {
  using inv_logit_tail_refs::ref;
  for (const auto& r : inv_logit_tail_refs::refs()) {
    stan::math::var a = r.x;
    stan::math::inv_logit(a).grad();
    inv_logit_tail_refs::expect_rel(r.d1, a.adj(), r.x);
    stan::math::recover_memory();
  }
  inv_logit_tail_refs::expect_fwd_derivatives(
      [](const auto& x) { return stan::math::inv_logit(x); },
      [](const ref& r) { return r.d1; }, [](const ref& r) { return r.d2; });
}

TEST(mathMixMatFun, invLogitDerivativesLargeArgsContainers) {
  using inv_logit_tail_refs::expect_rel;
  using stan::math::fvar;
  using stan::math::inv_logit;
  using stan::math::sum;
  using stan::math::var;
  std::vector<inv_logit_tail_refs::ref> ps = inv_logit_tail_refs::refs();
  int n = ps.size();
  Eigen::VectorXd x(n);
  Eigen::VectorXd d1(n);
  Eigen::VectorXd d2(n);
  for (int i = 0; i < n; ++i) {
    x(i) = ps[i].x;
    d1(i) = ps[i].d1;
    d2(i) = ps[i].d2;
  }

  stan::math::var_value<Eigen::VectorXd> xv(x);
  sum(inv_logit(xv)).grad();
  for (int i = 0; i < n; ++i) {
    expect_rel(d1(i), xv.adj()(i), x(i));
  }
  stan::math::recover_memory();

  stan::math::var_value<Eigen::RowVectorXd> xr(x.transpose());
  sum(inv_logit(xr)).grad();
  for (int i = 0; i < n; ++i) {
    expect_rel(d1(i), xr.adj()(i), x(i));
  }
  stan::math::recover_memory();

  // inv_logit' is even
  Eigen::MatrixXd x2(n, 2);
  x2 << x, -x;
  stan::math::var_value<Eigen::MatrixXd> xmv(x2);
  sum(inv_logit(xmv)).grad();
  for (int i = 0; i < n; ++i) {
    expect_rel(d1(i), xmv.adj()(i, 0), x(i));
    expect_rel(d1(i), xmv.adj()(i, 1), -x(i));
  }
  stan::math::recover_memory();

  Eigen::Matrix<var, Eigen::Dynamic, 1> xm = x;
  sum(inv_logit(xm)).grad();
  for (int i = 0; i < n; ++i) {
    expect_rel(d1(i), xm(i).adj(), x(i));
  }
  stan::math::recover_memory();

  std::vector<var> xs(x.data(), x.data() + n);
  sum(inv_logit(xs)).grad();
  for (int i = 0; i < n; ++i) {
    expect_rel(d1(i), xs[i].adj(), x(i));
  }
  stan::math::recover_memory();

  Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> xf(n);
  for (int i = 0; i < n; ++i) {
    xf(i) = fvar<var>(x(i), 1);
  }
  Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> yf = inv_logit(xf);
  for (int i = 0; i < n; ++i) {
    expect_rel(d1(i), yf(i).d_.val(), x(i));
  }
  sum(yf).d_.grad();
  for (int i = 0; i < n; ++i) {
    expect_rel(d2(i), xf(i).val_.adj(), x(i));
  }
  stan::math::recover_memory();
}
