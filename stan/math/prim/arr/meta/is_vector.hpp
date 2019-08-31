#ifndef STAN_MATH_PRIM_ARR_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_ARR_META_IS_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <type_traits>
#include <vector>

namespace stan {

template <typename T>
struct is_vector<std::vector<T>> : std::true_type {
  typedef T type;
};

template <typename T>
struct is_std_vector<std::vector<T>> : std::true_type {
  typedef T type;
};

}  // namespace stan
#endif
