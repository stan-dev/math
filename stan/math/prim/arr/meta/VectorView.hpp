#ifndef STAN_MATH_ARR_SCAL_META_VECTORVIEW_HPP
#define STAN_MATH_ARR_SCAL_META_VECTORVIEW_HPP

#include <stan/math/prim/scal/meta/VectorView.hpp>
#include <vector>

namespace stan {

  template <typename T>
  class VectorView<std::vector<T>, true, false> {
  public:
    typedef typename scalar_type<T>::type scalar_t;

    template <typename X>
    explicit VectorView(X& x) : x_(&x[0]) { }

    scalar_t& operator[](int i) {
      return x_[i];
    }

  private:
    scalar_t* x_;
  };

  template <typename T>
  class VectorView<const std::vector<T>, true, false> {
  public:
    typedef typename boost::add_const<typename scalar_type<T>::type>::type
    scalar_t;

    template <typename X>
    explicit VectorView(X& x) : x_(&x[0]) { }

    scalar_t& operator[](int i) const {
      return x_[i];
    }
  private:
    scalar_t* x_;
  };

}
#endif
