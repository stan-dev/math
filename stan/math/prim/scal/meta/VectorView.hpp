#ifndef STAN_MATH_PRIM_SCAL_META_VECTORVIEW_HPP
#define STAN_MATH_PRIM_SCAL_META_VECTORVIEW_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <boost/type_traits.hpp>
#include <stdexcept>

namespace stan {

  /**
   * VectorView is a template expression that is constructed with a
   * container or scalar, which it then allows to be used as an array
   * using <code>operator[]</code>.
   *
   * For a scalar value, any index returns the reference or pointer
   * used to construct the view.
   *
   * For a container, the index returns a reference to the position in
   * the underlying container used to construct the view.  WARNING: 
   * There is no bounds checking for container indices and they will
   * segfault if accessed beyond their boundaries.
   *
   * The first use is to read arguments to prob functions as vectors,
   * even if scalars, so they can be read by common code (and scalars
   * automatically broadcast up to behave like vectors) : VectorView
   * of immutable const array of double* (no allocation).
   *
   * The second use is to build up derivatives into common storage :
   * VectorView of mutable shared array (no allocation because
   * allocated on auto-diff arena memory).
   *
   * Because it deals with references to its inputs, it is up to the
   * client of VectorView to ensure that the container being wrapped
   * is not modified while the VectorView is in use in such a way as
   * to disrupt the indexing.  Similarly, because it deals with
   * references, it cannot be constructed with a literal or expression.
   *
   * @tparam T  Type of scalar or container being wrapped.
   * @tparam is_array True if underlying type T can be indexed with
   * operator[].
   * @tparam throw_if_accessed True if the behavior is to throw an
   * exception whenever <code>operator[]</code> is called.
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
