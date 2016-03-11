#ifndef STAN_MATH_PRIM_SCAL_META_VECTORVIEW_HPP
#define STAN_MATH_PRIM_SCAL_META_VECTORVIEW_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <boost/type_traits.hpp>
#include <stdexcept>

namespace stan {

  /**
   *  VectorView is a template metaprogram that takes its argument and
   *  allows it to be used like a vector. There are three template parameters
   *  - T: Type of the thing to be wrapped. For example,
   *       double, var, vector<double>, etc.
   *  - is_array: Boolean variable indicating whether the underlying type is
   *      an array.
   *  - throw_if_accessed: Boolean variable indicating whether this instance
   *      should not be used and should throw if operator[] is used.
   *
   *  For a scalar value, it broadcasts the single value when using
   *  operator[].
   *
   *  For a vector, operator[] looks into the value passed in.
   *  Note: this is not safe. It is possible to read past the size of
   *  an array.
   *
   *  Uses:
   *    Read arguments to prob functions as vectors, even if scalars, so
   *    they can be read by common code (and scalars automatically
   *    broadcast up to behave like vectors) : VectorView of immutable
   *    const array of double* (no allocation)
   *
   *    Build up derivatives into common storage : VectorView of
   *    mutable shared array (no allocation because allocated on
   *    auto-diff arena memory)
   */
  template <typename T,
            bool is_array = stan::is_vector_like<T>::value,
            bool throw_if_accessed = false>
  class VectorView {
  public:
    typedef typename
    boost::conditional<boost::is_const<T>::value,
                       typename boost::add_const<
                         typename scalar_type<T>::type>::type,
                       typename scalar_type<T>::type>::type scalar_t;

    template <typename X>
    explicit VectorView(X x) {
      throw std::logic_error("VectorView: the default template "
                             "specialization not implemented");
    }

    scalar_t& operator[](int i) {
      throw std::logic_error("VectorView: the default template "
                             "specialization not implemented");
    }

    scalar_t& operator[](int i) const {
      throw std::logic_error("VectorView: the default template "
                             "specialization not implemented");
    }
  };


  template <typename T, bool is_array>
  class VectorView<T, is_array, true> {
  public:
    typedef typename
    boost::conditional<boost::is_const<T>::value,
                       typename boost::add_const<
                         typename scalar_type<T>::type>::type,
                       typename scalar_type<T>::type>::type scalar_t;
    VectorView() { }

    template <typename X>
    explicit VectorView(X x) { }

    scalar_t& operator[](int i) {
      throw std::logic_error("VectorView: this cannot be accessed");
    }

    scalar_t& operator[](int i) const {
      throw std::logic_error("VectorView: this cannot be accessed");
    }
  };

  // this covers non-vectors: double
  template <typename T>
  class VectorView<T, false, false> {
  public:
    typedef typename
    boost::conditional<boost::is_const<T>::value,
                       typename boost::add_const<
                         typename scalar_type<T>::type>::type,
                       typename scalar_type<T>::type>::type scalar_t;

    explicit VectorView(scalar_t& x) : x_(&x) { }

    explicit VectorView(scalar_t* x) : x_(x) { }

    scalar_t& operator[](int i) {
      return *x_;
    }

    scalar_t& operator[](int i) const {
      return *x_;
    }
  private:
    scalar_t* x_;
  };


  // this covers raw memory: double*
  template <typename T>
  class VectorView<T, true, false> {
  public:
    typedef typename
    boost::conditional<boost::is_const<T>::value,
                       typename boost::add_const<
                         typename scalar_type<T>::type>::type,
                       typename scalar_type<T>::type>::type scalar_t;

    explicit VectorView(scalar_t* x) : x_(x) { }

    scalar_t& operator[](int i) {
      return x_[i];
    }

    scalar_t& operator[](int i) const {
      return x_[i];
    }

  private:
    scalar_t* x_;
  };
}
#endif
