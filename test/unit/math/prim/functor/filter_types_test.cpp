#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/functor/filter_types.hpp>
#include <stan/math/prim/meta/contains_autodiff.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

TEST(MathPrimFunctorFilterTypes, autodiff_predicate_preserves_nested_tuples) {
  using stan::math::fvar;
  std::vector<fvar<double>> a{
      fvar<double>(1.0, 0.0),
      fvar<double>(2.0, 0.0),
  };
  std::vector<int> b{3, 4};
  fvar<double> c(5.0, 0.0);
  int x = 1;
  auto nested = std::forward_as_tuple(b, c);
  auto args = std::forward_as_tuple(x, a, nested);

  auto filtered = stan::math::filter_types<stan::contains_autodiff>(args);
  using expected_t = std::tuple<std::vector<fvar<double>>&, std::tuple<fvar<double>&>>;
  static_assert(std::is_same_v<decltype(filtered), expected_t>,
                "filter_types should preserve nested tuple structure.");

  std::get<0>(filtered)[0].val_ = 7.0;
  EXPECT_FLOAT_EQ(7.0, a[0].val_);
  std::get<0>(std::get<1>(filtered)).val_ = -9.0;
  EXPECT_FLOAT_EQ(-9.0, c.val_);
}

TEST(MathPrimFunctorFilterTypes, no_match_returns_empty_tuple) {
  int x = 1;
  std::tuple<double, double> y{2.0, 3.0};
  std::vector<int> z{4, 5};
  auto args = std::forward_as_tuple(x, y, z);
  auto filtered = stan::math::filter_types<stan::contains_autodiff>(args);
  static_assert(std::tuple_size<std::decay_t<decltype(filtered)>>::value == 0,
                "Expected empty tuple when predicate matches no entries.");
}
