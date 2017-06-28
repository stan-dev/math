#include <stan/math.hpp>
#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>
#include <vector>
#include <iostream>

struct op_divide_f {
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type
  apply(const T1& x1, const T2& x2) {
    return x1 / x2;
  }
};

TEST(mathMixCore, operatorDivision) {
  using stan::math::test::test_ad;

  std::vector<double> xs;
  xs.push_back(0.5);
  xs.push_back(0);
  xs.push_back(-1.3);
  xs.push_back(stan::math::positive_infinity());
  xs.push_back(stan::math::negative_infinity());
  xs.push_back(stan::math::not_a_number());

  for (size_t i = 0; i < xs.size(); ++i)
    for (size_t j = 0; j < xs.size(); ++j)
      test_ad<op_divide_f>(xs[i], xs[j], xs[i] / xs[j]);
}

