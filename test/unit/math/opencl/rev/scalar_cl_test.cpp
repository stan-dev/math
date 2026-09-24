#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <type_traits>

using stan::math::var;
using stan::math::opencl::ScalarCl;
using stan::math::opencl::to_host;

TEST(ScalarClVar, traits) {
  EXPECT_TRUE((std::is_same<stan::scalar_type_t<ScalarCl<var>>, var>::value));
  EXPECT_TRUE(
      (std::is_same<stan::value_type_t<const ScalarCl<var>&>, var>::value));
  EXPECT_TRUE((stan::is_rev_kernel_expression<ScalarCl<var>>::value));
  EXPECT_TRUE((stan::is_prim_or_rev_kernel_expression<ScalarCl<var>>::value));
  EXPECT_FALSE((stan::is_kernel_expression<ScalarCl<var>>::value));
  EXPECT_FALSE(
      (stan::is_nonscalar_prim_or_rev_kernel_expression<ScalarCl<var>>::value));
  EXPECT_FALSE((stan::is_var<ScalarCl<var>>::value));
  EXPECT_FALSE((stan::is_stan_scalar<ScalarCl<var>>::value));
  EXPECT_TRUE((stan::is_autodiff_v<ScalarCl<var>>));
  EXPECT_FALSE((std::is_convertible<ScalarCl<var>, var>::value));
  EXPECT_FALSE((std::is_convertible<var, ScalarCl<var>>::value));
  EXPECT_FALSE((std::is_convertible<double, ScalarCl<var>>::value));
}

TEST(ScalarClVar, views_are_prim_device_scalars) {
  EXPECT_TRUE((stan::is_prim_scalar_cl<ScalarCl<const double&>>::value));
  EXPECT_TRUE((stan::is_prim_scalar_cl<ScalarCl<double&>>::value));
  EXPECT_TRUE(
      (std::is_same<stan::scalar_type_t<ScalarCl<double&>>, double>::value));
  EXPECT_TRUE((stan::is_kernel_expression<ScalarCl<const double&>>::value));
}

TEST_F(AgradRev, ScalarClVar_constructors) {
  ScalarCl<var> zero;
  EXPECT_EQ(to_host(zero.val()), 0.0);
  EXPECT_EQ(to_host(zero.adj()), 0.0);
  ScalarCl<var> x(2.5);
  EXPECT_EQ(to_host(x.val()), 2.5);
  ScalarCl<double> d(-3.0);
  ScalarCl<var> y(d);
  EXPECT_EQ(to_host(y.val()), -3.0);
  EXPECT_NE(y.val().buffer()(), d.buffer()());
}

TEST_F(AgradRev, ScalarClVar_cpu_var_round_trip) {
  var a = 1.5;
  ScalarCl<var> x(a);
  var out = to_host(x);
  EXPECT_EQ(out.val(), 1.5);
  out.grad();
  EXPECT_EQ(to_host(x.adj()), 1.0);
  EXPECT_EQ(a.adj(), 1.0);
}

TEST_F(AgradRev, ScalarClVar_value_and_adjoint_views) {
  ScalarCl<var> x(4.0);
  EXPECT_TRUE(
      (std::is_same<decltype(value_of(x)), ScalarCl<const double&>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(adjoint_of(x)), ScalarCl<double&>>::value));
  EXPECT_EQ(value_of(x).buffer()(), x.val().buffer()());
  adjoint_of(x) += 2.0;
  EXPECT_EQ(to_host(x.adj()), 2.0);
  EXPECT_EQ(stan::math::value_of(ScalarCl<double>(1.0)).matrix().rows(), 1);
}

namespace {
template <typename F>
void expect_binary_grad(F f, double a_val, double b_val) {
  // CPU reference
  var a_ref = a_val;
  var b_ref = b_val;
  var res_ref = f(a_ref, b_ref);
  res_ref.grad();
  double res_ref_val = res_ref.val();
  double a_ref_adj = a_ref.adj();
  double b_ref_adj = b_ref.adj();
  stan::math::recover_memory();

  // device var, device var
  {
    var a = a_val;
    var b = b_val;
    var res = to_host(f(ScalarCl<var>(a), ScalarCl<var>(b)));
    res.grad();
    EXPECT_NEAR(res.val(), res_ref_val, 1e-12);
    EXPECT_NEAR(a.adj(), a_ref_adj, 1e-12);
    EXPECT_NEAR(b.adj(), b_ref_adj, 1e-12);
    stan::math::recover_memory();
  }
  // device var, host double
  {
    var a = a_val;
    var res = to_host(f(ScalarCl<var>(a), b_val));
    res.grad();
    EXPECT_NEAR(res.val(), res_ref_val, 1e-12);
    EXPECT_NEAR(a.adj(), a_ref_adj, 1e-12);
    stan::math::recover_memory();
  }
  // device double, device var
  {
    var b = b_val;
    var res = to_host(f(ScalarCl<double>(a_val), ScalarCl<var>(b)));
    res.grad();
    EXPECT_NEAR(res.val(), res_ref_val, 1e-12);
    EXPECT_NEAR(b.adj(), b_ref_adj, 1e-12);
    stan::math::recover_memory();
  }
  // CPU var, device double
  {
    var a = a_val;
    var res = to_host(f(a, ScalarCl<double>(b_val)));
    res.grad();
    EXPECT_NEAR(res.val(), res_ref_val, 1e-12);
    EXPECT_NEAR(a.adj(), a_ref_adj, 1e-12);
    stan::math::recover_memory();
  }
  // device var, CPU var
  {
    var a = a_val;
    var b = b_val;
    var res = to_host(f(ScalarCl<var>(a), b));
    res.grad();
    EXPECT_NEAR(res.val(), res_ref_val, 1e-12);
    EXPECT_NEAR(a.adj(), a_ref_adj, 1e-12);
    EXPECT_NEAR(b.adj(), b_ref_adj, 1e-12);
    stan::math::recover_memory();
  }
}

template <typename F>
void expect_unary_grad(F f, double a_val) {
  var a_ref = a_val;
  var res_ref = f(a_ref);
  res_ref.grad();
  double res_ref_val = res_ref.val();
  double a_ref_adj = a_ref.adj();
  stan::math::recover_memory();

  var a = a_val;
  var res = to_host(f(ScalarCl<var>(a)));
  res.grad();
  EXPECT_NEAR(res.val(), res_ref_val, 1e-12);
  EXPECT_NEAR(a.adj(), a_ref_adj, 1e-12);
  stan::math::recover_memory();
}
}  // namespace

TEST(ScalarClVar, binary_operators) {
  expect_binary_grad([](const auto& a, const auto& b) { return a + b; }, 1.5,
                     -0.5);
  expect_binary_grad([](const auto& a, const auto& b) { return a - b; }, 1.5,
                     -0.5);
  expect_binary_grad([](const auto& a, const auto& b) { return a * b; }, 1.5,
                     -0.5);
  expect_binary_grad([](const auto& a, const auto& b) { return a / b; }, 1.5,
                     -0.5);
  expect_binary_grad(
      [](const auto& a, const auto& b) { return a * b + a / b - b; }, 2.0, 3.0);
}

TEST(ScalarClVar, unary_functions) {
  using stan::math::exp;
  using stan::math::log;
  using stan::math::sqrt;
  using stan::math::square;
  expect_unary_grad([](const auto& a) { return -a; }, 1.25);
  expect_unary_grad([](const auto& a) { return exp(a); }, 0.75);
  expect_unary_grad([](const auto& a) { return log(a); }, 0.75);
  expect_unary_grad([](const auto& a) { return sqrt(a); }, 2.0);
  expect_unary_grad([](const auto& a) { return square(a); }, -1.5);
  expect_unary_grad(
      [](const auto& a) { return exp(square(a)) * log(a) + 2.0 * a; }, 1.1);
}

TEST_F(AgradRev, ScalarClVar_copies_share_identity) {
  var a = 3.0;
  ScalarCl<var> x(a);
  ScalarCl<var> y = x;
  EXPECT_EQ(x.node().vi_, y.node().vi_);
  var res = to_host(x * y);
  res.grad();
  EXPECT_NEAR(a.adj(), 6.0, 1e-12);
}

TEST_F(AgradRev, ScalarClVar_rebinding_keeps_earlier_uses) {
  var a = 2.0;
  ScalarCl<var> x(a);
  ScalarCl<var> z = x * 3.0;  // uses the node x refers to now
  auto* old_node = x.node().vi_;
  x += 1.0;  // rebinds x to a new node
  EXPECT_NE(x.node().vi_, old_node);
  x *= x;  // x = (a + 1)^2
  var res = to_host(z + x);
  EXPECT_NEAR(res.val(), 6.0 + 9.0, 1e-12);
  res.grad();
  EXPECT_NEAR(a.adj(), 3.0 + 2.0 * 3.0, 1e-12);
}

TEST_F(AgradRev, ScalarClVar_repeated_conversions) {
  var a = 1.5;
  ScalarCl<var> x(a);
  var r1 = to_host(x);
  var r2 = to_host(x * 2.0);
  var res = r1 * r2;  // 2 a^2
  res.grad();
  EXPECT_NEAR(a.adj(), 4.0 * 1.5, 1e-12);
}

TEST_F(AgradRev, ScalarClVar_repeated_gradients) {
  var a = 0.5;
  ScalarCl<var> x(a);
  var res = to_host(exp(x) * x);
  double expected = std::exp(0.5) * 1.5;
  res.grad();
  EXPECT_NEAR(a.adj(), expected, 1e-12);
  stan::math::set_zero_all_adjoints();
  EXPECT_EQ(to_host(x.adj()), 0.0);
  EXPECT_EQ(a.adj(), 0.0);
  res.grad();
  EXPECT_NEAR(a.adj(), expected, 1e-12);
}

TEST_F(AgradRev, ScalarClVar_nested) {
  var a = 2.0;
  ScalarCl<var> x(a);
  stan::math::start_nested();
  {
    var b = 3.0;
    var inner = to_host(ScalarCl<var>(b) * x);
    inner.grad();
    EXPECT_NEAR(b.adj(), 2.0, 1e-12);
  }
  stan::math::recover_memory_nested();
  stan::math::set_zero_all_adjoints();
  EXPECT_EQ(to_host(x.adj()), 0.0);
  var res = to_host(square(x));
  res.grad();
  EXPECT_NEAR(a.adj(), 4.0, 1e-12);
}

TEST_F(AgradRev, ScalarClVar_captured_device_double_is_a_copy) {
  var a = 2.0;
  ScalarCl<double> d(5.0);
  ScalarCl<var> x(a);
  ScalarCl<var> y = x * d;
  d = ScalarCl<double>(100.0);  // must not change the derivative of y
  var res = to_host(y);
  res.grad();
  EXPECT_NEAR(a.adj(), 5.0, 1e-12);
}
#endif
