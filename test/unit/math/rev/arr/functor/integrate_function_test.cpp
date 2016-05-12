#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>

struct f1 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const T2& y, std::ostream* msgs) const {
    return exp(x) + y;
  }
};


struct f2 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const T2& y, std::ostream* msgs) const {
    return exp(y*cos(2*3.141593*x)) + y;
  }
};

struct f3 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const std::vector<T2>& y,
             std::ostream* msgs) const {
    return exp(x) + pow(y[0], 2.5) + 2*pow(y[1], 3) + 2*y[2];
  }
};

struct g3 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const std::vector<T2>& y, const int ii,
             std::ostream* msgs) const {
    if (ii == 1)
      return 2.5 * pow(y[0], 1.5);
    else if (ii == 2)
      return 6 * pow(y[1], 2.);
    else
      return 2.;
  }
};


TEST(StanMath_integrate_function, test1) {
  using stan::math::integrate_function;

  f1 if1;

  EXPECT_FLOAT_EQ(integrate_function(if1, .2, .7, stan::math::var(.5), 0).val(), 0.7923499+.25);

}


TEST(StanMath_integrate_function, finite_diff) {
  using stan::math::integrate_function;

  {
  f1 if1;

  AVAR a = .6;
  AVAR f = integrate_function(if1, .2, .7, a, 0);
  EXPECT_FLOAT_EQ(integrate_function(if1, .2, .7, .6, 0), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x,g);

  EXPECT_FLOAT_EQ((integrate_function(if1, .2, .7, .6+1e-6, 0) -
    integrate_function(if1, .2, .7, .6-1e-6, 0))/2e-6, g[0]);
  }
  {
  f2 if2;

  AVAR a = 0.68;
  AVAR f = integrate_function(if2, 0., 1.1, a, 0);
  EXPECT_FLOAT_EQ(integrate_function(if2, 0., 1.1, .68, 0), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x,g);

  EXPECT_FLOAT_EQ((integrate_function(if2, 0., 1.1, .68+1e-6, 0) -
    integrate_function(if2, 0., 1.1, .68-1e-6, 0))/2e-6, g[0]);
  }
  {
  f3 if3;

  AVAR a = 0.68;
  AVAR b = 0.38;
  AVAR c = 0.78;
  AVEC vec = createAVEC(a, b, c);
  AVAR f = integrate_function(if3, 0., 1.1, vec, 0);

  VEC g;
  double p1;
  double p2;
  f.grad(vec, g);

  std::vector<double> vecd = value_of(vec);
  EXPECT_FLOAT_EQ(integrate_function(if3, 0., 1.1, vecd, 0), f.val());

  vecd[0] += 1e-6;
  p1 = integrate_function(if3, 0., 1.1, vecd, 0);
  vecd[0] -= 2e-6;
  p2 = integrate_function(if3, 0., 1.1, vecd, 0);

  EXPECT_FLOAT_EQ((p1 - p2)/2e-6, g[0]);


  g3 ig3;
  stan::math::set_zero_all_adjoints();
  integrate_function_grad(if3, ig3, 0., 1.1, vec, 0).grad();
  EXPECT_FLOAT_EQ((p1 - p2)/2e-6, a.adj());
  }
}
