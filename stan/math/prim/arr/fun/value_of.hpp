#ifndef STAN_MATH_PRIM_ARR_FUN_VALUE_OF_HPP
#define STAN_MATH_PRIM_ARR_FUN_VALUE_OF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <algorithm>
#include <vector>
#include <utility>
#include <cstddef>

namespace stan {
namespace math {

/**
 * Convert a std::vector of type T to a std::vector of doubles.
 *
 * T must implement value_of_rec. See
 * test/math/fwd/mat/fun/value_of_rec.cpp for fvar and var usage.
 *
 * @tparam T Scalar type in std::vector
 * @param[in] x std::vector to be converted
 * @return std::vector of values
 **/
template <typename T, enable_if_std_vector<T>...,
          enable_if_double_or_int<scalar_type_decay_t<T>>...>
inline auto&& value_of(T&& x) {
  return std::forward<T>(x);
}

template <typename T, enable_if_std_vector<T>...,
          enable_if_arithmetic<scalar_type_decay_t<T>>...,
          enable_if_not_double_or_int<scalar_type_decay_t<T>>...>
inline auto value_of(T&& x) {
  std::vector<double> x_dbl(x.size());
  std::copy(x.begin(), x.end(), x_dbl.begin());
  return x_dbl;
}

}  // namespace math
}  // namespace stan

#endif
