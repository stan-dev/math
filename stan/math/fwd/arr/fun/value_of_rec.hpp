#ifndef STAN_MATH_FWD_ARR_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_FWD_ARR_FUN_VALUE_OF_REC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/value_of_rec.hpp>
#include <algorithm>
#include <utility>
#include <vector>
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
template <typename T, enable_if_std_vector<T>* = nullptr,
          enable_if_fvar<scalar_type_decay_t<T>>* = nullptr>
inline auto value_of_rec(T&& x) {
  std::vector<double> result(x.size());
  std::transform(std::forward<T>(x).begin(), std::forward<T>(x).end(),
                 result.begin(), [](auto&& x) { return value_of_rec(x.val_); });
  return result;
}

}  // namespace math
}  // namespace stan

#endif
