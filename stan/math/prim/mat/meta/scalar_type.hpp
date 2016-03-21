#ifndef STAN_MATH_PRIM_MAT_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_MAT_META_SCALAR_TYPE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <stan/math/prim/mat/meta/value_type.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>

namespace stan {

  template <typename T>
  struct scalar_type<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > {
    typedef typename scalar_type<T>::type type;
  };

}
#endif
