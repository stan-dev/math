#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/meta/contains_autodiff.hpp>
#include <stan/math/prim/meta/filter_index_sequence.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <utility>
#include <vector>

TEST(MathMetaFilterIndexSequence, integral_filtering) {
  using in_seq_t = std::index_sequence<0, 1, 2, 3>;
  using expected_t = std::index_sequence<0, 2>;
  using filtered_t = stan::math::filter_index_sequence_t<
      std::is_integral, in_seq_t, int, double, long, float>;
  EXPECT_SAME_TYPE(expected_t, filtered_t);
}

TEST(MathMetaFilterIndexSequence, empty_filtering) {
  using in_seq_t = std::index_sequence<0, 1, 2>;
  using expected_t = std::index_sequence<>;
  using filtered_t = stan::math::filter_index_sequence_t<
      std::is_integral, in_seq_t, double, float, long double>;
  EXPECT_SAME_TYPE(expected_t, filtered_t);
}

TEST(MathMetaFilterIndexSequence, autodiff_filtering) {
  using stan::math::fvar;
  using in_seq_t = std::index_sequence<0, 1, 2>;
  using expected_t = std::index_sequence<1>;
  using filtered_t = stan::math::filter_index_sequence_t<
      stan::contains_autodiff, in_seq_t, double, fvar<double>,
      std::vector<int>>;
  EXPECT_SAME_TYPE(expected_t, filtered_t);
}
