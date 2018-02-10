#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>

struct f1 {
  template <typename T1, typename T2, typename T3>
  inline
  typename stan::return_type<T1, T2, T3, int>::type
  operator()(const T1& x, const T2& y, const T3& x_r, int x_i,
             std::ostream& msgs) const {
    return exp(x) + y;
  }
};

struct f2 {
  template <typename T1, typename T2, typename T3>
  inline
  typename stan::return_type<T1, T2, T3, int>::type
  operator()(const T1& x, const T2& y, const T3& x_r, int x_i,
             std::ostream& msgs) const {
    return exp(y * cos(2 * 3.141593 * x)) + y;
  }
};

struct f3 {
  template <typename T1, typename T2, typename T3>
  inline
  typename stan::return_type<T1, T2, T3, int>::type
  operator()(const T1& x, const std::vector<T2>& y, const T3& x_r,
             int x_i, std::ostream& msgs) const {
    return exp(x) + pow(y[0], 2.5) + 2 * pow(y[1], 3) + 2 * y[2];
  }
};

struct g3 {
  template <typename T1, typename T2, typename T3>
  inline
  typename stan::return_type<T1, T2, T3, int>::type
  operator()(const T1& x, const std::vector<T2>& y,
             const T3& x_r, int x_i,
             const int ii, std::ostream& msgs) const {
    if (ii == 1)
      return 2.5 * pow(y[0], 1.5);
    else if (ii == 2)
      return 6 * pow(y[1], 2.);
    else
      return 2.;
  }
};


struct f4 {
  template <typename T1, typename T2, typename T3>
  inline
  typename stan::return_type<T1, T2, T3, int>::type
  operator()(const T1& x, const std::vector<T2>& y,
             const std::vector<T3>& x_r, const std::vector<int>& x_i,
             std::ostream& msgs) const {
    return exp(x) + pow(y[0], 2.5) + 2 * pow(y[1], 3) + 2 * y[2]
           + x_r[0] + 2*x_r[1] + 3*x_i[0] + 4*x_i[1];
  }
};


std::ostringstream msgs;

TEST(StanMath_integrate_1d_tsc, test1) {
  using stan::math::integrate_1d_tsc;

  f1 if1;
  double x_r;
  int x_i;

  EXPECT_FLOAT_EQ(
      integrate_1d_tsc(if1, .2, .7, stan::math::var(.5), x_r, x_i, msgs).val(),
      0.7923499 + .25);
}

TEST(StanMath_integrate_1d_tsc, finite_diff) {
  using stan::math::integrate_1d_tsc;

  double x_r;
  int x_i;

  {
    f1 if1;

    AVAR a = .6;
    AVAR f = integrate_1d_tsc(if1, .2, .7, a, x_r, x_i, msgs);
    EXPECT_FLOAT_EQ(integrate_1d_tsc(if1, .2, .7, .6, x_r, x_i, msgs), f.val());

    AVEC x = createAVEC(a);
    VEC g;
    f.grad(x, g);

    EXPECT_FLOAT_EQ((integrate_1d_tsc(if1, .2, .7, .6 + 1e-6, x_r, x_i, msgs)
                     - integrate_1d_tsc(if1, .2, .7, .6 - 1e-6, x_r, x_i, msgs))
                        / 2e-6,
                    g[0]);
  }
  {
    f2 if2;

    AVAR a = 0.68;
    AVAR f = integrate_1d_tsc(if2, 0., 1.1, a, x_r, x_i, msgs);
    EXPECT_FLOAT_EQ(integrate_1d_tsc(if2, 0., 1.1, .68, x_r, x_i, msgs), f.val());

    AVEC x = createAVEC(a);
    VEC g;
    f.grad(x, g);

    EXPECT_FLOAT_EQ((integrate_1d_tsc(if2, 0., 1.1, .68 + 1e-6, x_r, x_i, msgs)
                     - integrate_1d_tsc(if2, 0., 1.1, .68 - 1e-6, x_r, x_i, msgs))
                        / 2e-6,
                    g[0]);
  }
  {
    f3 if3;

    AVAR a = 0.68;
    AVAR b = 0.38;
    AVAR c = 0.78;
    AVEC vec = createAVEC(a, b, c);
    AVAR f = integrate_1d_tsc(if3, 0., 1.1, vec, x_r, x_i, msgs);

    VEC g;
    double p1;
    double p2;
    f.grad(vec, g);

    std::vector<double> vecd = value_of(vec);
    EXPECT_FLOAT_EQ(integrate_1d_tsc(if3, 0., 1.1, vecd, x_r, x_i, msgs), f.val());

    vecd[0] += 1e-6;
    p1 = integrate_1d_tsc(if3, 0., 1.1, vecd, x_r, x_i, msgs);
    vecd[0] -= 2e-6;
    p2 = integrate_1d_tsc(if3, 0., 1.1, vecd, x_r, x_i, msgs);

    EXPECT_FLOAT_EQ((p1 - p2) / 2e-6, g[0]);

    g3 ig3;
    stan::math::set_zero_all_adjoints();
    integrate_1d_tscg(if3, ig3, 0., 1.1, vec, x_r, x_i, msgs).grad();
    EXPECT_FLOAT_EQ((p1 - p2) / 2e-6, a.adj());
  }
  {
    f4 if4;

    std::vector<double> x_r(2);
    std::vector<int> x_i(2);
    x_r[0] = 3.2;
    x_r[1] = 6.33;
    x_i[0] = 13;
    x_i[1] = -12;

    AVAR a = 0.68;
    AVAR b = 0.38;
    AVAR c = 0.78;
    AVEC vec = createAVEC(a, b, c);

    std::vector<double> vecd = value_of(vec);

    AVAR f = integrate_1d_tsc(if4, 0., 1.1, vec, x_r, x_i, msgs);
    EXPECT_FLOAT_EQ(integrate_1d_tsc(if4, 0., 1.1, vecd, x_r, x_i, msgs), f.val());
  }
}
