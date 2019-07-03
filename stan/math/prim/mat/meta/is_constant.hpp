#ifndef STAN_MATH_PRIM_MAT_META_IS_CONSTANT_HPP
#define STAN_MATH_PRIM_MAT_META_IS_CONSTANT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/is_constant.hpp>

namespace stan {
/**
 * Defines a public enum named value and sets it to true
 * if the type of the elements in the provided Eigen Matrix
 * is constant, false otherwise. This is used in
 * the is_constant_all metaprogram.
 *
 * @tparam T type of the elements in the Eigen Matrix
 * @tparam R number of rows in the Eigen Matrix
 * @tparam C number of cols in the eigen Matrix
 */
template <typename T, int R, int C>
struct is_constant<Eigen::Matrix<T, R, C> > {
  enum { value = is_constant<T>::value };
};

/**
 * Defines a public enum named value and sets it to true
 * if the type of the elements in the provided Eigen Block
 * is constant, false otherwise. This is used in
 * the is_constant_all metaprogram.
 * @tparam type of the elements in the Eigen Block
 */
template <typename T>
struct is_constant<Eigen::Block<T> > {
  enum { value = is_constant<T>::value };
};

}  // namespace stan
#endif
