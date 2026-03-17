#ifndef STAN_TEST_UNIT_MATH_MIX_CONSTRAINT_LUB_CONSTRAIN_HELPERS_HPP
#define STAN_TEST_UNIT_MATH_MIX_CONSTRAINT_LUB_CONSTRAIN_HELPERS_HPP

namespace lub_constrain_tests {
template <typename T1, typename T2, typename T3>
inline void expect(const T1& x, const T2& lb, const T3& ub) {
  auto f1 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<false>(std::forward<decltype(x)>(x), std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub), lp);
  };
  auto f2 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<true>(std::forward<decltype(x)>(x), std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub), lp);
  };
  auto f3 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain<true>(std::forward<decltype(x)>(x), std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub), lp);
    return lp;
  };
  auto f4 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    auto xx = stan::math::lub_constrain<true>(std::forward<decltype(x)>(x), std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub), lp);
    return stan::math::add(lp, stan::math::sum(std::move(xx)));
  };
  auto f5 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<false>(std::forward<decltype(x)>(x), std::forward_as_tuple(std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub)), lp);
  };
  auto f6 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<true>(std::forward<decltype(x)>(x), std::forward_as_tuple(std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub)), lp);
  };

  stan::test::expect_ad(f1, x, lb, ub);
  stan::test::expect_ad(f2, x, lb, ub);
  stan::test::expect_ad(f3, x, lb, ub);
  stan::test::expect_ad(f4, x, lb, ub);
  stan::test::expect_ad(f5, x, lb, ub);
  stan::test::expect_ad(f6, x, lb, ub);
}
template <typename T1, typename T2, typename T3>
inline void expect_vec(const T1& x, const T2& lb, const T3& ub) {
  auto f1 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<false>(std::forward<decltype(x)>(x), std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub), lp);
  };
  auto f2 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<true>(std::forward<decltype(x)>(x), std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub), lp);
  };
  auto f3 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain<true>(std::forward<decltype(x)>(x), std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub), lp);
    return lp;
  };
  auto f4 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    auto xx = stan::math::eval(stan::math::lub_constrain<true>(std::forward<decltype(x)>(x), std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub), lp));
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> xx_acc = 0;
    for (size_t i = 0; i < xx.size(); ++i) {
      xx_acc += stan::math::sum(xx[i]);
    }
    return stan::math::add(lp, xx_acc);
  };
  auto f5 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<false>(std::forward<decltype(x)>(x), std::make_tuple(std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub)), lp);
  };
  auto f6 = [](auto&& x, auto&& lb, auto&& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<true>(std::forward<decltype(x)>(x), std::make_tuple(std::forward<decltype(lb)>(lb), std::forward<decltype(ub)>(ub)), lp);
  };

  stan::test::expect_ad(f1, x, lb, ub);
  stan::test::expect_ad(f2, x, lb, ub);
  stan::test::expect_ad(f3, x, lb, ub);
  stan::test::expect_ad(f4, x, lb, ub);
  stan::test::expect_ad(f5, x, lb, ub);
  stan::test::expect_ad(f6, x, lb, ub);
}
}  // namespace lub_constrain_tests

#endif
