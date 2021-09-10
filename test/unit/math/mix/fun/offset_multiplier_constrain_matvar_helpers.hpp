#ifndef STAN_TEST_UNIT_MATH_MIX_FUN_OFFSET_MULTIPLIER_MAT_CONSTRAIN_HELPERS_HPP
#define STAN_TEST_UNIT_MATH_MIX_FUN_OFFSET_MULTIPLIER_MAT_CONSTRAIN_HELPERS_HPP

namespace offset_multiplier_constrain_tests {
template <typename T1, typename T2, typename T3>
void expect_matvar(const T1& x, const T2& mu, const T3& sigma) {
  auto f1 = [](const auto& x, const auto& mu, const auto& sigma) {
    stan::return_type_t<decltype(x), decltype(mu), decltype(sigma)> lp = 0;
    return stan::math::offset_multiplier_constrain<false>(x, mu, sigma, lp);
  };
  auto f2 = [](const auto& x, const auto& mu, const auto& sigma) {
    stan::return_type_t<decltype(x), decltype(mu), decltype(sigma)> lp = 0;
    return stan::math::offset_multiplier_constrain<true>(x, mu, sigma, lp);
  };
  auto f3 = [](const auto& x, const auto& mu, const auto& sigma) {
    stan::return_type_t<decltype(x), decltype(mu), decltype(sigma)> lp = 0;
    stan::math::offset_multiplier_constrain<true>(x, mu, sigma, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& mu, const auto& sigma) {
    using lub_t
        = stan::return_type_t<decltype(x), decltype(mu), decltype(sigma)>;
    lub_t lp = 0;
    auto xx = stan::math::offset_multiplier_constrain<true>(x, mu, sigma, lp);
    return stan::math::add(lp, lub_t(stan::math::sum(xx)));
  };

  stan::test::expect_ad_matvar(f1, x, mu, sigma);
  stan::test::expect_ad_matvar(f2, x, mu, sigma);
  stan::test::expect_ad_matvar(f3, x, mu, sigma);
  stan::test::expect_ad_matvar(f4, x, mu, sigma);
}
template <typename T1, typename T2, typename T3>
void expect_vec_matvar(const T1& x, const T2& mu, const T3& sigma) {
  auto f1 = [](const auto& x, const auto& mu, const auto& sigma) {
    stan::return_type_t<decltype(x), decltype(mu), decltype(sigma)> lp = 0;
    return stan::math::offset_multiplier_constrain<false>(x, mu, sigma, lp);
  };
  auto f2 = [](const auto& x, const auto& mu, const auto& sigma) {
    stan::return_type_t<decltype(x), decltype(mu), decltype(sigma)> lp = 0;
    return stan::math::offset_multiplier_constrain<true>(x, mu, sigma, lp);
  };
  auto f3 = [](const auto& x, const auto& mu, const auto& sigma) {
    stan::return_type_t<decltype(x), decltype(mu), decltype(sigma)> lp = 0;
    stan::math::offset_multiplier_constrain<true>(x, mu, sigma, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& mu, const auto& sigma) {
    using lub_t
        = stan::return_type_t<decltype(x), decltype(mu), decltype(sigma)>;
    lub_t lp = 0;
    auto xx = stan::math::offset_multiplier_constrain<true>(x, mu, sigma, lp);
    lub_t xx_acc = 0;
    for (size_t i = 0; i < xx.size(); ++i) {
      xx_acc += stan::math::sum(xx[i]);
    }
    return stan::math::add(lp, xx_acc);
  };

  stan::test::expect_ad_matvar(f1, x, mu, sigma);
  stan::test::expect_ad_matvar(f2, x, mu, sigma);
  stan::test::expect_ad_matvar(f3, x, mu, sigma);
  stan::test::expect_ad_matvar(f4, x, mu, sigma);
}
}  // namespace offset_multiplier_constrain_tests

#endif
