#ifndef STAN_TEST_UNIT_MATH_MIX_FUN_LUB_CONSTRAIN_HELPERS_HPP
#define STAN_TEST_UNIT_MATH_MIX_FUN_LUB_CONSTRAIN_HELPERS_HPP

namespace lub_constrain_tests {
template <typename T1, typename T2, typename T3>
void expect(const T1& x, const T2& lb, const T3& ub) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<false>(x, lb, ub, lp);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<true>(x, lb, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain<true>(x, lb, ub, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    auto xx = stan::math::lub_constrain<true>(x, lb, ub, lp);
    return stan::math::add(lp, stan::math::sum(xx));
  };

  stan::test::expect_ad(f1, x, lb, ub);
  stan::test::expect_ad(f2, x, lb, ub);
  stan::test::expect_ad(f3, x, lb, ub);
  stan::test::expect_ad(f4, x, lb, ub);
}
template <typename T1, typename T2, typename T3>
void expect_vec(const T1& x, const T2& lb, const T3& ub) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<false>(x, lb, ub, lp);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<true>(x, lb, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain<true>(x, lb, ub, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    auto xx = stan::math::lub_constrain<true>(x, lb, ub, lp);
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> xx_acc = 0;
    for (size_t i = 0; i < xx.size(); ++i) {
      xx_acc += stan::math::sum(xx[i]);
    }
    return stan::math::add(lp, xx_acc);
  };

  stan::test::expect_ad(f1, x, lb, ub);
  stan::test::expect_ad(f2, x, lb, ub);
  stan::test::expect_ad(f3, x, lb, ub);
  stan::test::expect_ad(f4, x, lb, ub);
}
}  // namespace lub_constrain_tests

#endif
