#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>

TEST(AgradMixOperatorEqual, FvarVar) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> x(0.5,1.3);
  fvar<var> y(1.5,1.0);
  fvar<var> z(0.5,1.3);

  EXPECT_TRUE(x == z);
  EXPECT_FALSE(x == y);
  EXPECT_FALSE(z == y);
}

#include <stan/math.hpp>
#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>
#include <vector>

struct op_equal_f {
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type
  apply(const T1& x1, const T2& x2) {
    return (x1 == x2);
  }
};

TEST(mathMixCore, operatorEqual) {
  using stan::math::test::test_ad;

  test_ad<op_equal_f>(1.2, 1.0, 1.2 == 1.0);

  std::vector<double> xs;
  xs.push_back(0.5);
  xs.push_back(0);
  xs.push_back(-1.3);
  xs.push_back(stan::math::positive_infinity());
  xs.push_back(stan::math::negative_infinity());
  xs.push_back(stan::math::not_a_number());

  for (size_t i = 0; i < xs.size(); ++i) {
    for (size_t j = 0; j < xs.size(); ++j) {
      bool test_derivs = xs[i] != xs[j];
      test_ad<op_equal_f>(xs[i], xs[j], xs[i] == xs[j], test_derivs);
    }
  }
}

