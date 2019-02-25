#ifndef STAN_MATH_PRIM_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_META_SCALAR_TYPE_HPP

#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/metaar_type.hpp>

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
#ifndef STAN_MATH_PRIM_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_META_SCALAR_TYPE_HPP

#include <stanh/prim/metaar_type.hpp>
#include <vector>

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
#endif
#ifndef STAN_MATH_PRIM_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_META_SCALAR_TYPE_HPP

#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/metaar_type.hpp>

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
