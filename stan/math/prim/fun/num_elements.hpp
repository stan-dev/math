#ifndef STAN_MATH_PRIM_FUN_NUM_ELEMENTS_HPP
#define STAN_MATH_PRIM_FUN_NUM_ELEMENTS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Returns the number of elements in the specified, non-nested input
 *
 * @tparam T type of input
 * @param v argument
 * @return number of contained arguments
 */
template <typename T, require_not_std_vector_vt<is_container, T>* = nullptr,
          require_not_tuple_t<T>* = nullptr>
inline size_t num_elements(const T& v) {
  return math::size(v);
}

/**
 * Returns the number of elements in the specified vector
 *
 * @tparam T type of elements in the vector
 * @param v argument vector
 * @return number of contained arguments
 */
template <typename T, require_std_vector_vt<is_container, T>* = nullptr>
inline size_t num_elements(const T& v) {
  size_t size = 0;
  for (auto&& v_val : v) {
    size += num_elements(v_val);
  }
  return size;
}

template <typename T, require_tuple_t<T>* = nullptr>
inline size_t num_elements(const T& v) {
  size_t size = 0;
  math::apply(
      [&size](const auto&... args) {
        static_cast<void>(
            std::initializer_list<int>{(size += num_elements(args), 0)...});
      },
      v);
  return size;
}

}  // namespace math
}  // namespace stan

#endif
