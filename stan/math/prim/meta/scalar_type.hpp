#ifndef STAN_MATH_PRIM_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_META_SCALAR_TYPE_HPP

#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <vector>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
/**
 * Metaprogram structure to determine the base scalar type
 * of a template argument.
 *
 * <p>This base class should be specialized for structured types.
 *
 * @tparam T Type of object.
 */
template <typename T>
struct scalar_type {
  typedef T type;
};

template <typename T>
struct scalar_type<T*> {
  typedef typename scalar_type<T>::type type;
};

}  // namespace stan

namespace stan {
template <typename T>
struct scalar_type<std::vector<T> > {
  typedef typename scalar_type<T>::type type;
};

template <typename T>
struct scalar_type<const std::vector<T> > {
  typedef typename scalar_type<T>::type type;
};

template <typename T>
struct scalar_type<std::vector<T>&> {
  typedef typename scalar_type<T>::type type;
};

template <typename T>
struct scalar_type<const std::vector<T>&> {
  typedef typename scalar_type<T>::type type;
};
}  // namespace stan

namespace stan {

template <typename T, int R, int C>
struct scalar_type<Eigen::Matrix<T, R, C> > {
  typedef typename scalar_type<T>::type type;
};

template <typename T, int R, int C>
struct scalar_type<const Eigen::Matrix<T, R, C> > {
  typedef typename scalar_type<T>::type type;
};

template <typename T, int R, int C>
struct scalar_type<Eigen::Matrix<T, R, C>&> {
  typedef typename scalar_type<T>::type type;
};

template <typename T, int R, int C>
struct scalar_type<const Eigen::Matrix<T, R, C>&> {
  typedef typename scalar_type<T>::type type;
};

template <typename T>
struct scalar_type<Eigen::Block<T> > {
  typedef typename scalar_type<T>::type type;
};
}  // namespace stan
#endif
