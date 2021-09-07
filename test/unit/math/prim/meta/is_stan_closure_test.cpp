#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/closure_adapter.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, IsStanClosure) {
  auto lambda = [](auto msg) { return 0.0; };
  auto cl = stan::math::from_lambda(lambda);
  EXPECT_FALSE((stan::is_stan_closure<decltype(lambda)>::value));
  EXPECT_TRUE((stan::is_stan_closure<decltype(cl)>::value));
}

TEST(MathMetaPrim, ClosureReturnType) {
  EXPECT_SAME_TYPE(const std::vector<double>&,
                   stan::closure_return_type<std::vector<double>, true>::type);
  EXPECT_SAME_TYPE(std::vector<double>,
                   stan::closure_return_type<std::vector<double>, false>::type);
}
