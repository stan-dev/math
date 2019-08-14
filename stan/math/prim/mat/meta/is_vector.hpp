#ifndef STAN_MATH_PRIM_MAT_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_MAT_META_IS_VECTOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/enable_if_eigen_vector.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>

namespace stan {

/**
  * Specialization of is_vector for Eigen vectors
  */
template <typename T>
struct is_vector<T, enable_if_eigen_vector<T>> : std::true_type {
  typedef std::decay_t<T> type;
};

}  // namespace stan
#endif
