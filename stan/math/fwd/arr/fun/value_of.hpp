#ifndef STAN_MATH_FWD_ARR_FUN_VALUE_OF_HPP
#define STAN_MATH_FWD_ARR_FUN_VALUE_OF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/scal/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <algorithm>
#include <utility>
#include <vector>
#include <cstddef>

namespace stan {
namespace math {

/**
 * Convert a std::vector of type T to a std::vector of the partial type.
 *
 * T must implement value_of_rec. See
 * test/math/fwd/mat/fun/value_of_rec.cpp for fvar and var usage.
 *
 * @tparam T Scalar type in std::vector
 * @param[in] x std::vector to be converted
 * @return std::vector of partial type.
 **/
template <typename T, require_std_vector_fvar<T>...>
inline auto value_of(T&& x) {
  std::vector<partials_type_t<scalar_type_decay_t<T>>> result(x.size());
  std::transform(std::forward<T>(x).begin(), std::forward<T>(x).end(),
                 result.begin(), [](auto&& x) { return x.val_; });
  return result;
}

}  // namespace math
}  // namespace stan

#endif
