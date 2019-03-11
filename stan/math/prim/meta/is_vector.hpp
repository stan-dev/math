#ifndef STAN_MATH_PRIM_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_META_IS_VECTOR_HPP

#include <stan/math/prim/meta/is_vector.hpp>
#include <vector>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {

// FIXME: use boost::type_traits::remove_all_extents to
//        extend to array/ptr types

template <typename T>
struct is_vector {
  enum { value = 0 };
  typedef T type;
};
}  // namespace stan

namespace stan {

// FIXME: use boost::type_traits::remove_all_extents to
//   extend to array/ptr types

template <typename T>
struct is_vector<const T> {
  enum { value = is_vector<T>::value };
  typedef T type;
};
template <typename T>
struct is_vector<std::vector<T> > {
  enum { value = 1 };
  typedef T type;
};
}  // namespace stan

namespace stan {

// FIXME: use boost::type_traits::remove_all_extents to
//        extend to array/ptr types

template <typename T>
struct is_vector<Eigen::Matrix<T, Eigen::Dynamic, 1> > {
  enum { value = 1 };
  typedef T type;
};
template <typename T>
struct is_vector<Eigen::Matrix<T, 1, Eigen::Dynamic> > {
  enum { value = 1 };
  typedef T type;
};
template <typename T>
struct is_vector<Eigen::Block<T> > {
  enum { value = 1 };
  typedef T type;
};
}  // namespace stan
#endif
