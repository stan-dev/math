#ifndef STAN_MATH_PRIM_MAT_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_MAT_META_IS_VECTOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <type_traits>

namespace stan {

template <typename T>
struct is_vector<T, std::enable_if_t<is_eigen_vector<T>::value>>
    : std::true_type {
  typedef T type;
};

}  // namespace stan
#endif
