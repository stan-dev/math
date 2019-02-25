#ifndef STAN_MATH_PRIM_META_IS_CONSTANT_STRUCT_HPP
#define STAN_MATH_PRIM_META_IS_CONSTANT_STRUCT_HPP

#include <stanh/prim/meta/is_constant_struct.hpp>
#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/meta/is_constant.hpp>

namespace stan {

template <typename T, int R, int C>
struct is_constant_struct<Eigen::Matrix<T, R, C> > {
  enum { value = is_constant_struct<T>::value };
};

template <typename T>
struct is_constant_struct<Eigen::Block<T> > {
  enum { value = is_constant_struct<T>::value };
};

}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_META_IS_CONSTANT_STRUCT_HPP
#define STAN_MATH_PRIM_META_IS_CONSTANT_STRUCT_HPP

#include <stanh/prim/meta/is_constant.hpp>
#include <stanh/prim/meta/is_constant_struct.hpp>
#include <vector>

namespace stan {

template <typename T>
struct is_constant_struct<std::vector<T> > {
  enum { value = is_constant_struct<T>::value };
};

}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_META_IS_CONSTANT_STRUCT_HPP
#define STAN_MATH_PRIM_META_IS_CONSTANT_STRUCT_HPP

#include <stanh/prim/meta/is_constant_struct.hpp>
#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/meta/is_constant.hpp>

namespace stan {

template <typename T, int R, int C>
struct is_constant_struct<Eigen::Matrix<T, R, C> > {
  enum { value = is_constant_struct<T>::value };
};

template <typename T>
struct is_constant_struct<Eigen::Block<T> > {
  enum { value = is_constant_struct<T>::value };
};

}  // namespace stan
#endif
