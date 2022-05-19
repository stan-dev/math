#ifndef STAN_MATH_PRIM_FUN_AS_ARRAY_OR_SCALAR_HPP
#define STAN_MATH_PRIM_FUN_AS_ARRAY_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns specified input value.
 *
 * @tparam Scalar Type of element.
 * @param v Specified value.
 * @return Same value.
 */
template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
inline Scalar as_array_or_scalar(Scalar&& v) {
  return std::forward<Scalar>(v);
}

/**
 * Returns a reference to rvalue specified input value.
 *
 * @tparam Scalar Type of element.
 * @param v Specified value.
 * @return Same value.
 */
template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
inline Scalar& as_array_or_scalar(Scalar& v) {
  return v;
}

/**
 * Returns specified input value.
 *
 * @tparam EigArray Type of element.
 * @param v Specified value.
 * @return Same value.
 */
template <typename EigArray, require_eigen_array_t<EigArray>* = nullptr>
inline EigArray as_array_or_scalar(EigArray&& v) {
  return std::forward<EigArray>(v);
}

/**
 * Converts a matrix type to an array.
 *
 * @tparam Eig Type of \c Eigen \c Matrix or expression
 * @param v Specified \c Eigen \c Matrix or expression.
 * @return Matrix converted to an array.
 */
template <typename Eig, typename = require_eigen_t<Eig>,
          require_not_eigen_array_t<Eig>* = nullptr>
inline auto as_array_or_scalar(Eig&& v) {
  return make_holder([](auto& x) { return x.array(); }, std::forward<Eig>(v));
}

/**
 * Converts a std::vector type to an array.
 *
 * @tparam StdVec Type of scalar element.
 * @param v Specified vector.
 * @return Matrix converted to an array.
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr,
          require_not_std_vector_t<value_type_t<StdVec>>* = nullptr>
inline auto as_array_or_scalar(StdVec&& v) {
  using T_map
      = Eigen::Map<const Eigen::Array<value_type_t<StdVec>, Eigen::Dynamic, 1>>;
  return make_holder([](auto& x) { return T_map(x.data(), x.size()); },
                     std::forward<StdVec>(v));
}

/**
 * Converts an std::vector<std::vector> to an Eigen Array.
 * @tparam Container A standard vector with inner container of a standard vector
 *  with an inner stan scalar.
 * @param v specified vector of vectorised
 * @return An Eigen Array with dynamic rows and columns.
 */
template <typename Container, require_std_vector_vt<is_std_vector, Container>* = nullptr,
          require_std_vector_vt<is_stan_scalar, value_type_t<Container>>* = nullptr>
inline auto as_array_or_scalar(Container&& v) {
  Eigen::Array<scalar_type_t<Container>, -1, -1> ret(v.size(), v[0].size());
  for (size_t i = 0; i < v.size(); ++i) {
    ret.row(i) = Eigen::Map<const Eigen::Array<scalar_type_t<Container>, 1, -1>>(
        v[i].data(), v[i].size());
  }
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
