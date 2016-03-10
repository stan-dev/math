#ifndef STAN_MATH_MAT_SCAL_META_VECTORVIEW_HPP
#define STAN_MATH_MAT_SCAL_META_VECTORVIEW_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/VectorView.hpp>
#include <boost/type_traits.hpp>

namespace stan {

  template <typename T, int R, int C>
  class VectorView<Eigen::Matrix<T, R, C>, true, false> {
  public:
    typedef typename scalar_type<T>::type scalar_t;

    template <typename X>
    explicit VectorView(X& x) : x_(x.data()) { }

    scalar_t& operator[](int i) {
      return x_[i];
    }
  private:
    scalar_t* x_;
  };

  template <typename T, int R, int C>
  class VectorView<const Eigen::Matrix<T, R, C>, true, false> {
  public:
    typedef typename boost::add_const<typename scalar_type<T>::type>::type
    scalar_t;

    template <typename X>
    explicit VectorView(X& x) : x_(x.data()) { }

    scalar_t& operator[](int i) const {
      return x_[i];
    }
  private:
    scalar_t* x_;
  };

}
#endif
