#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <vector>

TEST(MathMetaPrim, is_tuple) {
  EXPECT_TRUE((stan::math::is_tuple<std::tuple<double>>::value));
  EXPECT_TRUE((stan::math::is_tuple<std::tuple<double>&>::value));
  EXPECT_TRUE((stan::math::is_tuple<std::tuple<double>&&>::value));
  EXPECT_FALSE((stan::math::is_tuple<double>::value));
  EXPECT_FALSE((stan::math::is_tuple<std::vector<double>>::value));
}

TEST(MathMetaPrim, all_tuple_elements) {
  using stan::is_std_vector;
  using stan::math::internal::all_tuple_elements;
  EXPECT_TRUE((all_tuple_elements<
               is_std_vector,
               std::tuple<std::vector<int>, std::vector<double>>>::value));
  EXPECT_FALSE(
      (all_tuple_elements<is_std_vector,
                          std::tuple<std::vector<int>, double>>::value));
  EXPECT_TRUE((all_tuple_elements<
               is_std_vector,
               std::tuple<std::vector<double>&, std::vector<double>&>>::value));
  EXPECT_FALSE((all_tuple_elements<is_std_vector, std::vector<double>>::value));
  EXPECT_FALSE((all_tuple_elements<is_std_vector, double>::value));
}
