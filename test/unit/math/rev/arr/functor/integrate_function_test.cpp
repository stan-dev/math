#include <stan/math/prim/mat/meta/get.hpp>
#include <stan/math/prim/arr/meta/get.hpp>
#include <stan/math/prim/mat/meta/length.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <stan/math/prim/arr/meta/is_vector.hpp>
#include <stan/math/prim/mat/meta/is_vector_like.hpp>
#include <stan/math/rev/scal/fun/pow.hpp>
#include <stan/math/rev/scal/fun/cos.hpp>
#include <stan/math/rev/scal/fun/exp.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <gtest/gtest.h>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <stan/math/prim/arr/functor/integrate_function.hpp>
#include <test/unit/util.hpp>


struct f1 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const T2& y) const {
    return exp(x) + y;
  }
};


struct f2 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const T2& y) const {
    return exp(y*cos(2*3.141593*x)) + y;
  }
};

struct f3 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const std::vector<T2>& y) const {
    return exp(x) + pow(y[0], 2.5) + 2*pow(y[1], 3) + 2*y[2];
  }
};


TEST(StanMath_integrate_function, test1) {
  using stan::math::integrate_function;

  f1 if1;

  EXPECT_FLOAT_EQ(integrate_function(if1, .2, .7, stan::math::var(.5)).val(), 0.7923499+.25);

}


TEST(StanMath_integrate_function, finite_diff) {
  using stan::math::integrate_function;

  {
  f1 if1;

  AVAR a = .6;
  AVAR f = integrate_function(if1, .2, .7, a);
  EXPECT_FLOAT_EQ(integrate_function(if1, .2, .7, .6), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x,g);

  EXPECT_FLOAT_EQ((integrate_function(if1, .2, .7, .6+1e-6) -
    integrate_function(if1, .2, .7, .6-1e-6))/2e-6, g[0]);
  }
  {
  f2 if2;

  AVAR a = 0.68;
  AVAR f = integrate_function(if2, 0., 1.1, a);
  EXPECT_FLOAT_EQ(integrate_function(if2, 0., 1.1, .68), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x,g);

  EXPECT_FLOAT_EQ((integrate_function(if2, 0., 1.1, .68+1e-6) -
    integrate_function(if2, 0., 1.1, .68-1e-6))/2e-6, g[0]);
  }
  {
  f3 if3;

  AVAR a = 0.68;
  AVAR b = 0.38;
  AVAR c = 0.78;
  AVEC vec = createAVEC(a, b, c);
  AVAR f = integrate_function(if3, 0., 1.1, vec);

  VEC g;
  double p1;
  double p2;
  f.grad(vec, g);

  std::vector<double> vecd = value_of(vec);
  EXPECT_FLOAT_EQ(integrate_function(if3, 0., 1.1, vecd), f.val());

  vecd[0] += 1e-6;
  p1 = integrate_function(if3, 0., 1.1, vecd);
  vecd[0] -= 2e-6;
  p2 = integrate_function(if3, 0., 1.1, vecd);

  EXPECT_FLOAT_EQ((p1 - p2)/2e-6, g[0]);
  }
}
