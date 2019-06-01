#ifndef TEST_UNIT_MATH_TEST_AD_HPP
#define TEST_UNIT_MATH_TEST_AD_HPP

#include <stan/math.hpp>
#include <test/unit/math/util.hpp>
#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>
#include <string>
#include <vector>

namespace stan {
namespace test {

template <typename G>
void expect_ad_derivatives(const G& g, const Eigen::VectorXd& x) {
  double gx = g(x);
  stan::math::test::test_gradient(g, x, gx);
  stan::math::test::test_gradient_fvar(g, x, gx);
  stan::math::test::test_hessian(g, x, gx);
  stan::math::test::test_hessian_fvar(g, x, gx);
  stan::math::test::test_grad_hessian(g, x, gx);
}

template <typename F, typename H, typename... Ts>
void expect_ad_helper(const F& f, const H& h, const Eigen::VectorXd& x,
                      Ts... xs) {
  size_t result_size;
  try {
    auto y = f(xs...);
    result_size = serialize<double>(y).size();
  } catch (...) {
    stan::math::test::expect_all_throw(h(0), x);
    return;
  }
  for (size_t i = 0; i < result_size; ++i)
    expect_ad_derivatives(h(i), x);
}

template <typename F, typename T1, typename T2>
void expect_ad_vv(const F& f, const T1& x1, const T2& x2) {
  auto h = [&](const int& i) {
    return [&](const auto& v) {
      auto ds = to_deserializer(v);
      auto x1ds = ds.read(x1);
      auto x2ds = ds.read(x2);
      return serialize_return(f(x1ds, x2ds))[i];
    };
  };
  expect_ad_helper(f, h, serialize_args(x1, x2), x1, x2);
}

template <typename F, typename T1, typename T2>
void expect_ad_vd(const F& f, const T1& x1, const T2& x2) {
  auto h = [&](const int& i) {
    return [&](const auto& v) {
      auto ds = to_deserializer(v);
      auto x1ds = ds.read(x1);
      return serialize_return(f(x1ds, x2))[i];
    };
  };
  Eigen::VectorXd x = serialize_args(x1);
  expect_ad_helper(f, h, serialize_args(x1), x1, x2);
}

template <typename F, typename T1, typename T2>
void expect_ad_dv(const F& f, const T1& x1, const T2& x2) {
  auto h = [&](const int& i) {
    return [&](const auto& v) {
      auto ds = to_deserializer(v);
      auto x2ds = ds.read(x2);
      return serialize_return(f(x1, x2ds))[i];
    };
  };
  expect_ad_helper(f, h, serialize_args(x2), x1, x2);
}

template <typename F, typename T>
void expect_ad_v(const F& f, const T& x1) {
  auto h = [&](const int& i) {
    return [&](const auto& v) {
      auto ds = to_deserializer(v);
      auto x1ds = ds.read(x1);
      return serialize_return(f(x1ds))[i];
    };
  };
  expect_ad_helper(f, h, serialize_args(x1), x1);
}

// CLIENT AUTODIFF TEST FUNCTIONS
//
// These are the public test functions that expect a functor
// f encapsulating a polymorphic call to a Stan function and
// a sequence of arguments with double-based scalars.
//
// These functions evaluate that all levels of functional autodiff
// provide the same answers as finite differences on the double-based
// implementations, including exception behavior.  To completely
// test a new differentiable function, a developer need only
//
// (a) independently test the double-based implementation, and
//
// (b) use these functions to provide arguments to test at all
//     levels of autodiff.
//
// The underlying heavy lifting is provided in the file
//     unit/math/mix/mat/util/autodiff_tester.hpp
//
// Example use cases are provided in this directory in file
//     test_ad_test.cpp.

// Unary function autodiff tester
template <typename F, typename T>
void expect_ad(const F& f, const T& x) {
  expect_ad_v(f, x);
}

// Binary function autodiff tester
template <typename F, typename T1, typename T2>
void expect_ad(const F& f, const T1& x1, const T2& x2) {
  expect_ad_vv(f, x1, x2);
  expect_ad_vd(f, x1, x2);
  expect_ad_dv(f, x1, x2);
}

}  // namespace test
}  // namespace stan

#endif
