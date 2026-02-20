#include <stan/math/mix.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <vector>

template <typename F>
inline void expect_match_autodiff(const F& f, Eigen::VectorXd x) {
  double fx_fd;
  Eigen::VectorXd grad_fx_fd;
  stan::math::finite_diff_gradient_auto(f, x, fx_fd, grad_fx_fd);

  double fx;
  Eigen::VectorXd grad_fx;
  stan::math::gradient(f, x, fx, grad_fx);

  EXPECT_FLOAT_EQ(fx, fx_fd);
  EXPECT_EQ(grad_fx.size(), grad_fx_fd.size());
  stan::test::expect_near_rel("expect_match_autodiff", grad_fx, grad_fx_fd,
                              stan::test::relative_tolerance(1e-7, 1e-9));
}

TEST(MathMixMatFunctor, FiniteDiffGradientAuto) {
  auto norm_fun
      = [](const auto& x) { return stan::math::normal_lpdf(x(0), x(1), x(2)); };

  for (double frac = 1e-40; frac <= 1e+40; frac *= 1e10) {
    Eigen::VectorXd y(3);
    y << 1 * frac, 2 * frac, 3;
    expect_match_autodiff(norm_fun, y);
  }

  for (double frac = 1; frac <= 1e+10; frac *= 10) {
    Eigen::VectorXd y(3);
    y << 0, 0, frac;
    expect_match_autodiff(norm_fun, y);
  }

  auto log_fun
      = [](const auto& x) { return stan::math::sum(stan::math::log(x)); };
  Eigen::VectorXd y(0);
  expect_match_autodiff(log_fun, y);

  Eigen::VectorXd z(1);
  z << 2;
  expect_match_autodiff(log_fun, z);

  Eigen::VectorXd w(5);
  w << 1, 2, 3, 4, 5;
  expect_match_autodiff(log_fun, w);
}

TEST(MathMixMatFunctor, FiniteDiffGradientAutoTupleNested) {
  using stan::math::fvar;
  using fvar_stdvec = std::vector<fvar<double>>;
  fvar_stdvec a(2);
  a[0] = fvar<double>(1.3, 0.0);
  a[1] = fvar<double>(-0.7, 0.0);
  using fvar_mat = Eigen::Matrix<fvar<double>, -1, 1>;
  fvar_mat b{{fvar<double>(2.0, 0.0), fvar<double>(-3.0, 0.0)}};
  const std::vector<int> x_i{5, 6};
  Eigen::VectorXd c{{4.0}};
  auto args = std::make_tuple(a,
    std::make_tuple(b, x_i),
    c);
  fvar_mat b_grad = fvar_mat::Zero(b.size());
  auto grads = std::make_tuple(
      fvar_stdvec(a.size(), 0.0),
      std::make_tuple(b_grad, std::vector<int>(x_i.size(), 0)),
      Eigen::VectorXd{{0.0}});

  auto f = [](const auto& arg_a, const auto& arg_b_tuple, const auto& arg_c) {
    const auto& arg_b = std::get<0>(arg_b_tuple);
    const auto& arg_i = std::get<1>(arg_b_tuple);
    return arg_a[0] * arg_a[1] + arg_b(0) * arg_b(1) + arg_c(0) * arg_i[0];
  };

  fvar<double> fx = 0;
  stan::math::finite_diff_gradient_auto(f, fx, args, grads);

  EXPECT_NEAR(13.09, fx.val(), 1e-8);
  EXPECT_NEAR(a[1].val_, std::get<0>(grads)[0].val_, 1e-8);
  EXPECT_NEAR(a[0].val_, std::get<0>(grads)[1].val_, 1e-8);
  EXPECT_NEAR(b(1).val_, std::get<0>(std::get<1>(grads))(0).val_, 1e-8);
  EXPECT_NEAR(b(0).val_, std::get<0>(std::get<1>(grads))(1).val_, 1e-8);
  EXPECT_EQ(0, std::get<1>(std::get<1>(grads))[0]);
  EXPECT_EQ(0, std::get<1>(std::get<1>(grads))[1]);
}

TEST(MathMixMatFunctor, FiniteDiffGradientAutoTupleRestoresInputs) {
  using stan::math::fvar;
  using fvar_stdvec = std::vector<fvar<double>>;
  fvar_stdvec a(2);
  a[0] = fvar<double>(1.3, 0.0);
  a[1] = fvar<double>(-0.7, 0.0);
  using fvar_mat = Eigen::Matrix<fvar<double>, -1, 1>;
  fvar_mat b{{fvar<double>(2.0, 0.0), fvar<double>(-3.0, 0.0)}};
  std::vector<int> x_i{5, 6};
  Eigen::VectorXd c{{4.0}};

  auto args = std::make_tuple(a, std::make_tuple(b, x_i), c);
  auto args_before = args;
  fvar_mat b_grad = fvar_mat::Zero(b.size());
  auto grads = std::make_tuple(
      fvar_stdvec(a.size(), 0.0), std::make_tuple(b_grad, x_i),
      Eigen::VectorXd::Zero(c.size()));

  auto f = [](const auto& arg_a, const auto& arg_b_tuple, const auto& arg_c) {
    const auto& arg_b = std::get<0>(arg_b_tuple);
    const auto& arg_i = std::get<1>(arg_b_tuple);
    return arg_a[0] * arg_a[1] + arg_b(0) * arg_b(1) + arg_c(0) * arg_i[0];
  };

  fvar<double> fx = 0;
  stan::math::finite_diff_gradient_auto(f, fx, args, grads);

  const auto& a_after = std::get<0>(args);
  const auto& a_orig = std::get<0>(args_before);
  for (std::size_t i = 0; i < a_after.size(); ++i) {
    EXPECT_FLOAT_EQ(a_orig[i].val_, a_after[i].val_);
    EXPECT_FLOAT_EQ(a_orig[i].d_, a_after[i].d_);
  }

  const auto& b_after = std::get<0>(std::get<1>(args));
  const auto& b_orig = std::get<0>(std::get<1>(args_before));
  for (Eigen::Index i = 0; i < b_after.size(); ++i) {
    EXPECT_FLOAT_EQ(b_orig(i).val_, b_after(i).val_);
    EXPECT_FLOAT_EQ(b_orig(i).d_, b_after(i).d_);
  }

  const auto& x_i_after = std::get<1>(std::get<1>(args));
  const auto& x_i_orig = std::get<1>(std::get<1>(args_before));
  EXPECT_EQ(x_i_orig, x_i_after);

  const auto& c_after = std::get<2>(args);
  const auto& c_orig = std::get<2>(args_before);
  ASSERT_EQ(c_orig.size(), c_after.size());
  for (Eigen::Index i = 0; i < c_after.size(); ++i) {
    EXPECT_FLOAT_EQ(c_orig(i), c_after(i));
  }
}
