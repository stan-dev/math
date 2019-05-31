#ifndef TEST_UNIT_MATH_TEST_AD_HPP
#define TEST_UNIT_MATH_TEST_AD_HPP

#include <stan/math.hpp>
#include <test/unit/math/util.hpp>
#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>
#include <string>
#include <vector>

namespace stan {
namespace test {

template <typename T>
Eigen::Matrix<T, -1, 1> to_eigen_vector(const std::vector<T>& x) {
  return Eigen::Map<const Eigen::Matrix<T, -1, 1>>(x.data(), x.size());
}

template <typename T, int R, int C>
std::vector<T> to_std_vector(const Eigen::Matrix<T, R, C>& x) {
  std::vector<T> y;
  y.reserve(x.size());
  for (int i = 0; i < x.size(); ++i)
    y.push_back(x(i));
  return y;
}

// g : Eigen::Matrix<T, -1, 1> -> T
template <typename G>
void expect_ad_helper(const G& g, const std::vector<double>& xs) {
  static const bool TEST_DERIVS = true;
  Eigen::VectorXd x = to_eigen_vector(xs);
  double gx;
  try {
    gx = g(x);
    stan::math::test::test_gradient(g, x, gx, TEST_DERIVS);
    stan::math::test::test_gradient_fvar(g, x, gx, TEST_DERIVS);
    stan::math::test::test_hessian(g, x, gx, TEST_DERIVS);
    stan::math::test::test_hessian_fvar(g, x, gx, TEST_DERIVS);
    stan::math::test::test_grad_hessian(g, x, gx, TEST_DERIVS);
  } catch (...) {
    stan::math::test::expect_all_throw(g, x);
  }
}

template <typename F, typename T1, typename T2>
void expect_ad_vv(const F& f, const T1& x1, const T2& x2) {
  std::vector<double> xs = serialize<double>(x1, x2);  // *****
  try {
    auto y = f(x1, x2);  // *****
    std::vector<double> ys = serialize<double>(y);
    for (size_t i = 0; i < ys.size(); ++i) {
      auto g = [&](const auto& v) {
        typedef typename scalar_type<decltype(v)>::type scalar_t;
        deserializer<scalar_t> ds(to_std_vector(v));
        auto x1ds = ds.read(x1);                       // *****
        auto x2ds = ds.read(x2);                       // *****
        return serialize<scalar_t>(f(x1ds, x2ds))[i];  // *****
      };
      expect_ad_helper(g, xs);
    }
  } catch (...) {
    std::cout << "Errorful inputs not yet testable.";
  }
}

template <typename F, typename T1, typename T2>
void expect_ad_dv(const F& f, const T1& x1, const T2& x2) {
  std::vector<double> xs = serialize<double>(x2);
  try {
    auto y = f(x1, x2);
    std::vector<double> ys = serialize<double>(y);
    for (size_t i = 0; i < ys.size(); ++i) {
      auto g = [&](const auto& v) {
        typedef typename scalar_type<decltype(v)>::type scalar_t;
        deserializer<scalar_t> ds(to_std_vector(v));
        auto x2ds = ds.read(x2);
        return serialize<scalar_t>(f(x1, x2ds))[i];
      };
      expect_ad_helper(g, xs);
    }
  } catch (...) {
    std::cout << "Errorful inputs not yet testable.";
  }
}

template <typename F, typename T1, typename T2>
void expect_ad_vd(const F& f, const T1& x1, const T2& x2) {
  std::vector<double> xs = serialize<double>(x1);
  try {
    auto y = f(x1, x2);
    std::vector<double> ys = serialize<double>(y);
    for (size_t i = 0; i < ys.size(); ++i) {
      auto g = [&](const auto& v) {
        typedef typename scalar_type<decltype(v)>::type scalar_t;
        deserializer<scalar_t> ds(to_std_vector(v));
        auto x1ds = ds.read(x1);
        return serialize<scalar_t>(f(x1ds, x2))[i];
      };
      expect_ad_helper(g, xs);
    }
  } catch (...) {
    std::cout << "Errorful inputs not yet testable.";
  }
}

// f : T -> U  and  x : T  for any Stan-legal types T, U
// x will be double-based in the call, but f must be polymorphic
template <typename F, typename T>
void expect_ad_v(const F& f, const T& x) {
  std::vector<double> xs = serialize<double>(x);
  try {
    auto y = f(x);
    std::vector<double> ys = serialize<double>(y);
    for (size_t i = 0; i < ys.size(); ++i) {
      // v : Matrix<T, -1, 1>
      // g : Matrix<T, -1, 1> -> T
      // T == decltype(v)::Scalar (will be instantiated to autodiff var)
      auto g = [&](const auto& v) {
        typedef typename scalar_type<decltype(v)>::type scalar_t;
        deserializer<scalar_t> ds(to_std_vector(v));
        return serialize<scalar_t>(f(ds.read(x)))[i];
      };
      expect_ad_helper(g, xs);
    }
  } catch (...) {
    std::cout << "Errorful inputs not yet testable.";
  }
}

// UNARY TEST FUNCTION
template <typename F, typename T>
void expect_ad(const F& f, const T& x) {
  expect_ad_v(f, x);
}

// BINARY TEST FUNCTION
template <typename F, typename T1, typename T2>
void expect_ad(const F& f, const T1& x1, const T2& x2) {
  expect_ad_vv(f, x1, x2);
  // expect_ad_vd(f, x1, x2);
  // expect_ad_dv(f, x1, x2);
}

}  // namespace test
}  // namespace stan

#endif
