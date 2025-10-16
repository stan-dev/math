#ifndef TEST_UNIT_MATH_REV_CORE_GRADABLE_HPP
#define TEST_UNIT_MATH_REV_CORE_GRADABLE_HPP

#include <stan/math/rev.hpp>
#include <stan/math/rev/fun/quad_form.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <vector>

struct gradable {
  std::vector<stan::math::var> x_;
  stan::math::var f_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> g_expected_;
  gradable(const std::vector<stan::math::var>& x, const stan::math::var& f,
           const Eigen::Matrix<double, Eigen::Dynamic, 1>& g_expected)
      : x_(x), f_(f), g_expected_(g_expected) {}
  void test() {
    std::vector<double> g;
    f_.grad(x_, g);
    EXPECT_EQ(g_expected_.size(), static_cast<int>(g.size()));
    for (int i = 0; i < g_expected_.size(); ++i)
      EXPECT_FLOAT_EQ(g_expected_(i), g[i]);
  }

  double adj() { return f_.adj(); }
};

gradable setup_quad_form() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::quad_form;
  using stan::math::var;
  using std::vector;

  Matrix<var, Dynamic, 1> u(3);
  u << 2, 3, 5;

  Matrix<var, Dynamic, Dynamic> S(3, 3);
  S << 7, 11, 13, 17, 19, 23, 29, 31, 37;

  vector<var> x;
  for (int i = 0; i < u.size(); ++i)
    x.push_back(u(i));
  for (int i = 0; i < S.size(); ++i)
    x.push_back(S(i));

  var f = quad_form(S, u);

  Matrix<double, 1, Dynamic> g_expected(12);
  g_expected << 322, 440, 616, 4, 6, 10, 6, 9, 15, 10, 15, 25;

  return gradable(x, f, g_expected);
}

gradable setup_simple() {
  stan::math::var a = 3;
  stan::math::var b = 7;
  std::vector<stan::math::var> x;
  x.push_back(a);
  x.push_back(b);
  stan::math::var f = 2 * a * b;
  Eigen::Matrix<double, Eigen::Dynamic, 1> g_expected(2);
  g_expected << 14, 6;
  return gradable(x, f, g_expected);
}

#endif
