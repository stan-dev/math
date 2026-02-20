#include <stan/math/prim/fun/zeroed_filtered_tuple.hpp>
#include <stan/math/prim/meta/contains_autodiff.hpp>
#include <test/unit/math/prim/util.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

namespace {

template <typename T>
struct is_any_fp_scalar_impl {
  static constexpr bool value
      = std::is_floating_point_v<stan::scalar_type_t<std::decay_t<T>>>;
};

template <typename... Types>
struct is_any_fp_scalar_impl<std::tuple<Types...>> {
  static constexpr bool value
      = (is_any_fp_scalar_impl<std::decay_t<Types>>::value || ...);
};

template <typename T, typename... VecArgs>
struct is_any_fp_scalar_impl<std::vector<T, VecArgs...>> {
  static constexpr bool value = is_any_fp_scalar_impl<std::decay_t<T>>::value;
};

template <typename T>
struct is_any_fp_scalar {
  static constexpr bool value = is_any_fp_scalar_impl<std::decay_t<T>>::value;
};

}  // namespace

TEST(MathPrimFunZeroedContainer, floating_point_predicate_filters_and_zeros) {
  std::vector<double> a{
      1.3,
      -0.7,
  };
  Eigen::Matrix<double, -1, 1> b(2);
  b << 2.0, -3.0;
  std::vector<int> x_i{5, 6};
  Eigen::VectorXd c{{4.0}};

  auto args = std::make_tuple(a, std::make_tuple(b, x_i), c);
  auto grads = stan::math::zeroed_filtered_tuple<is_any_fp_scalar>(args);

  using grads_t = std::decay_t<decltype(grads)>;
  EXPECT_EQ(3, std::tuple_size<grads_t>::value);

  const auto& a_grads = std::get<0>(grads);
  EXPECT_EQ(a.size(), a_grads.size());
  for (std::size_t i = 0; i < a_grads.size(); ++i) {
    EXPECT_FLOAT_EQ(0.0, a_grads[i]);
  }

  const auto& b_grads = std::get<0>(std::get<1>(grads));
  EXPECT_EQ(b.size(), b_grads.size());
  for (Eigen::Index i = 0; i < b_grads.size(); ++i) {
    EXPECT_FLOAT_EQ(0.0, b_grads(i));
  }

  const auto& c_grads = std::get<2>(grads);
  EXPECT_EQ(c.size(), c_grads.size());
  for (Eigen::Index i = 0; i < c_grads.size(); ++i) {
    EXPECT_FLOAT_EQ(0.0, c_grads(i));
  }
}

TEST(MathPrimFunZeroedContainer, integral_predicate_filters_and_zeros) {
  auto args = std::make_tuple(std::vector<int>{1, 2, 3}, 1.5,
                              std::make_tuple(2.5, 7));
  auto zeros = stan::math::zeroed_filtered_tuple<std::is_integral>(args);

  using zeros_t = std::decay_t<decltype(zeros)>;
  EXPECT_EQ(2, std::tuple_size<zeros_t>::value);

  const auto& vec_zero = std::get<0>(zeros);
  EXPECT_EQ(3, vec_zero.size());
  EXPECT_EQ(0, vec_zero[0]);
  EXPECT_EQ(0, vec_zero[1]);
  EXPECT_EQ(0, vec_zero[2]);

  const auto& nested_zero = std::get<0>(std::get<1>(zeros));
  EXPECT_EQ(0, nested_zero);
}

TEST(MathPrimFunZeroedContainer, nested_tuple_values_are_zeroed) {
  auto args = std::make_tuple(
      1.25, std::make_tuple(2.5, -3.0, 1, std::make_tuple(4.75)));
  auto zeros = stan::math::zeroed_filtered_tuple<std::is_floating_point>(args);

  using expected_t
      = std::tuple<double, std::tuple<double, double, std::tuple<double>>>;
  static_assert(std::is_same_v<std::decay_t<decltype(zeros)>, expected_t>,
                "zeroed_filtered_tuple should preserve nested tuple structure.");

  EXPECT_FLOAT_EQ(0.0, std::get<0>(zeros));
  EXPECT_FLOAT_EQ(0.0, std::get<0>(std::get<1>(zeros)));
  EXPECT_FLOAT_EQ(0.0, std::get<1>(std::get<1>(zeros)));
  EXPECT_FLOAT_EQ(0.0, std::get<0>(std::get<2>(std::get<1>(zeros))));
}

TEST(MathPrimFunZeroedContainer, tuple_containing_tuple_is_zeroed) {
  auto args = std::make_tuple(std::make_tuple(3.14), 3.14, std::make_tuple(1));
  auto zeros = stan::math::zeroed_filtered_tuple<std::is_floating_point>(args);

  using expected_t = std::tuple<std::tuple<double>, double>;
  static_assert(
      std::is_same_v<std::decay_t<decltype(zeros)>, expected_t>,
      "zeroed_filtered_tuple should support tuple<tuple<...>> inputs.");

  EXPECT_FLOAT_EQ(0.0, std::get<0>(std::get<0>(zeros)));
}

TEST(MathPrimFunZeroedContainer, empty_tuple_returns_empty_tuple) {
  auto zeros = stan::math::zeroed_filtered_tuple<std::is_integral>(std::make_tuple());
  static_assert(std::tuple_size<std::decay_t<decltype(zeros)>>::value == 0,
                "Expected empty tuple output for empty tuple input.");
}

TEST(MathPrimFunZeroedContainer, no_match_returns_empty_tuple) {
  auto args = std::make_tuple(1.0, std::make_tuple(2.0, 3.0));
  auto zeros = stan::math::zeroed_filtered_tuple<std::is_integral>(args);
  static_assert(std::tuple_size<std::decay_t<decltype(zeros)>>::value == 0,
                "Expected empty tuple when predicate matches no entries.");
}

TEST(MathPrimFunZeroedContainer, vector_of_tuples_preserves_shape) {
  using vec_tuple_t = std::vector<std::tuple<double, int>>;
  vec_tuple_t vec{{1.2, 3}, {-0.4, 5}};
  auto args = std::make_tuple(vec);
  auto zeros = stan::math::zeroed_filtered_tuple<is_any_fp_scalar>(args);

  using expected_t = std::tuple<std::vector<std::tuple<double>>>;
  static_assert(std::is_same_v<std::decay_t<decltype(zeros)>, expected_t>,
                "Expected compact tuple shape for vector-of-tuples.");

  const auto& vec_zero = std::get<0>(zeros);
  ASSERT_EQ(2, vec_zero.size());
  EXPECT_FLOAT_EQ(0.0, std::get<0>(vec_zero[0]));
  EXPECT_FLOAT_EQ(0.0, std::get<0>(vec_zero[1]));
}
