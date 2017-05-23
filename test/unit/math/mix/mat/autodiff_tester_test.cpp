#include <stan/math.hpp>
#include <test/unit/math/mix/mat/autodiff_tester.hpp>
#include <gtest/gtest.h>
#include <iostream>

struct foo_fun {
  template <typename T>
  T operator()(const Eigen::Matrix<T, -1, 1>& x) const {
    return 1.75 * x[0] * x[0];
  }
};

TEST(foo, bar) {
  foo_fun f;
  Eigen::VectorXd x(1);
  x << 2.0;
  stan::math::test::test_functor(f, x, 4 * 1.75);
}

struct bar_fun {
  template <typename T>
  T operator()(const Eigen::Matrix<T, -1, 1>& x) const {
    return 1.2 * x[0] * exp(x[1]);
  }
};

TEST(foo, bar2) {
  bar_fun f;
  Eigen::VectorXd x(2);
  x << 2.0, 3.0;
  stan::math::test::test_functor(f, x, 1.2 * 2.0 * std::exp(3.0));
}

TEST(foo, oldFramework) {
  bar_fun f;
  stan::math::test::ad_tester<bar_fun> t(f);
  Eigen::VectorXd x(2);
  x << 2.0, 3.0;
  t.good(x, 1.2 * 2.0 * std::exp(3.0));

  Eigen::VectorXd b(2);
  b << 2.0, 3.0;
  t.good(b, f(b));

  t.test();
}

struct baz_fun {
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type
  apply(const T1& x1, const T2& x2) {
    using std::exp;
    return -2.5 * x1 * x2 * exp(x2);
  }
};

TEST(foo, framework) {
  stan::math::test::ad_tester_binary<baz_fun, double, double> t;
  t.good(2.1, -1.3, -2.5 * 2.1 * (-1.3) * std::exp(-1.3));
  t.run_tests();
}
