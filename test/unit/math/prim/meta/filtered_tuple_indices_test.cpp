#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/meta/contains_autodiff.hpp>
#include <stan/math/prim/meta/filtered_tuple_indices.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

TEST(MathMetaFilteredTupleIndices, autodiff_indices) {
  using stan::math::fvar;
  using tuple_t
      = std::tuple<int, fvar<double>, std::vector<double>,
                   std::tuple<fvar<double>, double>>;
  using expected_t = std::index_sequence<1, 3>;
  using idxs_t = stan::math::filtered_tuple_indices_t<stan::contains_autodiff,
                                                      tuple_t>;
  EXPECT_SAME_TYPE(expected_t, idxs_t);
}

TEST(MathMetaFilteredTupleIndices, integral_indices_and_size) {
  using tuple_t = std::tuple<double, int, float, long>;
  using expected_t = std::index_sequence<1, 3>;
  using idxs_t
      = stan::math::filtered_tuple_indices_t<std::is_integral, tuple_t>;
  EXPECT_SAME_TYPE(expected_t, idxs_t);
  EXPECT_EQ(2, (stan::math::filtered_tuple_size_v<std::is_integral, tuple_t>));
}

TEST(MathMetaFilteredTupleIndices, tuple_decay) {
  using tuple_ref_t = const std::tuple<const int, double, const long>&;
  using expected_t = std::index_sequence<0, 2>;
  using idxs_t
      = stan::math::filtered_tuple_indices_t<std::is_integral, tuple_ref_t>;
  EXPECT_SAME_TYPE(expected_t, idxs_t);
}
