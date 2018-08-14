#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>

#include <vector>

struct op_gamma_p_f {
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type apply(
      const T1& x1, const T2& x2) {
    return stan::math::gamma_p(x1, x2);
  }
};

TEST(mathMixScalFunGammaP, GoodGammaP) {
  std::vector<double> xs;
  xs.push_back(2.5);
  xs.push_back(3.1);
  xs.push_back(1.3);
  xs.push_back(stan::math::not_a_number());
  stan::math::test::test_args<op_gamma_p_f, false>(xs, xs);
}

TEST(mathMixScalFunGammaP, PosInfGammaP) {
  std::vector<double> xs;
  xs.push_back(stan::math::positive_infinity());
  stan::math::test::test_args<op_gamma_p_f, false>(xs, xs);
}

TEST(mathMixScalFunGammaP, NaNGammaP) {
  std::vector<double> xs;
  xs.push_back(stan::math::not_a_number());
  stan::math::test::test_args<op_gamma_p_f, false>(xs, xs);
}

TEST(mathMixScalFunGammaP, NegInfGammaP) {
  std::vector<double> xs;
  xs.push_back(stan::math::negative_infinity());
  stan::math::test::test_args<op_gamma_p_f, false>(xs, xs);
}
