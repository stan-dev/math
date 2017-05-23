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

TEST(foo, framework) {
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
