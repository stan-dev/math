#ifndef STAN_MATH_FWD_CORE_FVAR_HPP
#define STAN_MATH_FWD_CORE_FVAR_HPP

#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/fwd/scal/meta/ad_promotable.hpp>
#include <stan/math/cplx/complex.hpp>
#include <ostream>
#include <type_traits>

namespace stan {
namespace math {

/**
 * This template class represents scalars used in forward-mode
 * automatic differentiation, which consist of values and
 * directional derivatives of the specified template type.  When
 * performing operations on instances of this class, all operands
 * should be either primitive integer or double values or dual
 * numbers representing derivatives in the same direction.  The
 * typical use case is to have a unit length directional
 * derivative in the direction of a single independent variable.
 *
 * By using reverse-mode automatic derivative variables,
 * second-order derivatives may
 * be calculated.  By using fvar&lt;<var&gt; instances,
 * third-order derivatives may be calculated.  These are called
 * mixed-mode automatic differentiation variable in Stan.
 *
 * Specialized functionals that perform differentiation on
 * functors may be found in the matrix subdirectories of the
 * reverse, forward, and mixed-mode directories.
 *
 * The <a
 * href="https://en.wikipedia.org/wiki/Automatic_differentiation">Wikipedia
 * page on automatic differentiation</a> describes how
 * forward-mode automatic differentiation works mathematically in
 * terms of dual numbers.
 *
 * @tparam T type of value and tangent
 */
template <typename T>
struct fvar {
  /**
   * The value of this variable.
   */
  T val_;

  /**
   * The tangent (derivative) of this variable.
   */
  T d_;

  /**
   * Return the value of this variable.
   *
   * @return value of this variable
   */
  T val() const { return val_; }

  /**
   * Return the tangent (derivative) of this variable.
   *
   * @return tangent of this variable
   */
  T tangent() const { return d_; }

  /**
   * Construct a forward variable with zero value and tangent.
   */
  fvar() : val_(0.0), d_(0.0) {}

  /**
   * Construct a forward variable with value and tangent set to
   * the value and tangent of the specified variable.
   *
   * @param[in] x variable to be copied
   */
  fvar(const fvar<T>& x) : val_(x.val_), d_(x.d_) {}

  /**
   * Construct a forward variable with the specified value and
   * zero tangent.
   *
   * @tparam V type of value (must be assignable to the value and
   *   tangent type T)
   * @param[in] v value
   */
  fvar(const T& v) : val_(v), d_(0.0) {  // NOLINT(runtime/explicit)
    if (unlikely(is_nan(v)))
      d_ = v;
  }

  /**
   * Construct a forward variable with the specified value and
   * zero tangent.
   *
   * @tparam V type of value (must be assignable to the value and
   *   tangent type T)
   * @param[in] v value
   */
  template <typename V, typename std::enable_if_t<
                            ad_promotable<V, T>::value>* = nullptr>
  fvar(const V& v) : val_(v), d_(0.0) {  // NOLINT(runtime/explicit)
    if (unlikely(is_nan(v)))
      d_ = v;
  }

  /**
   * Construct a forward variable with the specified value and
   * tangent.
   *
   * @tparam V type of value (must be assignable to the value and
   *   tangent type T)
   * @tparam D type of tangent (must be assignable to the value and
   *   tangent type T)
   * @param[in] v value
   * @param[in] d tangent
   */
  template <typename V, typename D>
  fvar(const V& v, const D& d) : val_(v), d_(d) {
    if (unlikely(is_nan(v)))
      d_ = v;
  }

  /**
   * Add the specified variable to this variable and return a
   * reference to this variable.
   *
   * @param[in] x2 variable to add
   * @return reference to this variable after addition
   */
  inline fvar<T>& operator+=(const fvar<T>& x2) {
    val_ += x2.val_;
    d_ += x2.d_;
    return *this;
  }

  /**
   * Add the specified value to this variable and return a
   * reference to this variable.
   *
   * @param[in] x2 value to add
   * @return reference to this variable after addition
   */
  inline fvar<T>& operator+=(double x2) {
    val_ += x2;
    return *this;
  }

  /**
   * Subtract the specified variable from this variable and return a
   * reference to this variable.
   *
   * @param[in] x2 variable to subtract
   * @return reference to this variable after subtraction
   */
  inline fvar<T>& operator-=(const fvar<T>& x2) {
    val_ -= x2.val_;
    d_ -= x2.d_;
    return *this;
  }

  /**
   * Subtract the specified value from this variable and return a
   * reference to this variable.
   *
   * @param[in] x2 value to add
   * @return reference to this variable after subtraction
   */
  inline fvar<T>& operator-=(double x2) {
    val_ -= x2;
    return *this;
  }

  /**
   * Multiply this variable by the the specified variable and
   * return a reference to this variable.
   *
   * @param[in] x2 variable to multiply
   * @return reference to this variable after multiplication
   */
  inline fvar<T>& operator*=(const fvar<T>& x2) {
    d_ = d_ * x2.val_ + val_ * x2.d_;
    val_ *= x2.val_;
    return *this;
  }

  /**
   * Multiply this variable by the the specified value and
   * return a reference to this variable.
   *
   * @param[in] x2 value to multiply
   * @return reference to this variable after multiplication
   */
  inline fvar<T>& operator*=(double x2) {
    val_ *= x2;
    d_ *= x2;
    return *this;
  }

  /**
   * Divide this variable by the the specified variable and
   * return a reference to this variable.
   *
   * @param[in] x2 variable to divide this variable by
   * @return reference to this variable after division
   */
  inline fvar<T>& operator/=(const fvar<T>& x2) {
    d_ = (d_ * x2.val_ - val_ * x2.d_) / (x2.val_ * x2.val_);
    val_ /= x2.val_;
    return *this;
  }

  /**
   * Divide this value by the the specified variable and
   * return a reference to this variable.
   *
   * @param[in] x2 value to divide this variable by
   * @return reference to this variable after division
   */
  inline fvar<T>& operator/=(double x2) {
    val_ /= x2;
    d_ /= x2;
    return *this;
  }

  /**
   * Increment this variable by one and return a reference to this
   * variable after the increment.
   *
   * @return reference to this variable after increment
   */
  inline fvar<T>& operator++() {
    ++val_;
    return *this;
  }

  /**
   * Increment this variable by one and return a reference to a
   * copy of this variable before it was incremented.
   *
   * @return reference to copy of this variable before increment
   */
  inline fvar<T> operator++(int /*dummy*/) {
    fvar<T> result(val_, d_);
    ++val_;
    return result;
  }

  /**
   * Decrement this variable by one and return a reference to this
   * variable after the decrement.
   *
   * @return reference to this variable after decrement
   */
  inline fvar<T>& operator--() {
    --val_;
    return *this;
  }

  /**
   * Decrement this variable by one and return a reference to a
   * copy of this variable before it was decremented.
   *
   * @return reference to copy of this variable before decrement
   */
  inline fvar<T> operator--(int /*dummy*/) {
    fvar<T> result(val_, d_);
    --val_;
    return result;
  }

  /**
   * Write the value of the specified variable to the specified
   * output stream, returning a reference to the output stream.
   *
   * @param[in,out] os stream for writing value
   * @param[in] v variable whose value is written
   * @return reference to the specified output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const fvar<T>& v) {
    return os << v.val_;
  }
};

// this should probably be moved to something _parallel_ to:
// stan/math/rev/scal/fun/boost_isfinite
// needed for fullPivLu() method on non-Hermitian Eigen matrix types
// used in stan's complex test
// called from Eigen::cplx::isfinite_impl via ADL
template <class T>
bool isfinite(fvar<T> const& v) {
  using boost::math::isfinite;
  using std::isfinite;
  return isfinite(v.val());
}

namespace cplx {

/**
 * Fvar for std::complex<fvar>.
 *
 * This class exists to unify fvar with the same extended-but-hidden interface
 * of stan::math::complex's class template that exists behind std::complex<var>.
 */
template <class T = double>
struct z_fvar : fvar<T> {
  template<class Z = double,
   std::enable_if_t<is_arith<Z>::value>* = nullptr>
  z_fvar(Z const& z = 0) : fvar<T>(z) {}  // NOLINT
};

/// helper type traits to avoid forward declarations in other headers
template <class T>
struct is_fr_var_helper<fvar<T>> : std::true_type {};
template <class T>
struct is_fr_var_helper<z_fvar<T>> : std::true_type {};
template <class T>
struct to_arith_helper<z_fvar<T>> {
  typedef fvar<to_arith_t<T>> type;
};

}  // namespace cplx

/// helper functions to avoid forward declarations in other headers
// declared here because it depends on ADL from an uncontrolled type
template <class T>
inline T rval_help(fvar<T> const& f) {
  return f.val();
}

}  // namespace math
}  // namespace stan

namespace std {

/** Template specialization to std::complex<fvar> to inherit the shared
 * extended-but-hidden interface of stan::math::complex<z_fvar>. In this way,
 * z_fvar is hidden from end user code.*/
// template <>
template <class T>
struct complex<stan::math::fvar<T>>
    : stan::math::cplx::complex<stan::math::cplx::z_fvar<T>> {
  /// inherit all ctors
  using stan::math::cplx::complex<stan::math::cplx::z_fvar<T>>::complex;
};

//override clang's division, which uses scalbn and logb
template <class T>
inline std::complex<stan::math::fvar<T>>
operator/ (
    std::complex<stan::math::cplx::z_fvar<T>> const& t,
    std::complex<stan::math::cplx::z_fvar<T>> const& u) {
  return stan::math::cplx::division(t, u);
}

}  // namespace std

namespace stan{
namespace math{

//these forwards allow movement of several (target) function
//template std::complex overloads out of the std namespace
//and into the stan::math::cplx namespace, even though the
//target function templates are already tightly constrained
template <class T, class U,
 std::enable_if_t<std::is_arithmetic<U>::value>* = nullptr>
inline auto
operator+(fvar<T> const& t, std::complex<U> const& u) {
  return stan::math::cplx::operator+(t,u);
}

template <class T, class U,
 std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
inline auto
operator+(std::complex<T> const& t, fvar<U> const& u) {
 return stan::math::cplx::operator+(t,u);
}

template <class T, class U,
 std::enable_if_t<std::is_arithmetic<U>::value>* = nullptr>
inline auto
operator-(fvar<T> const& t, std::complex<U> const& u) {
  return stan::math::cplx::operator-(t,u);
}

template <class T, class U,
 std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
inline auto
operator-(std::complex<T> const& t, fvar<U> const& u) {
 return stan::math::cplx::operator-(t,u);
}

template <class T, class U,
 std::enable_if_t<std::is_arithmetic<U>::value>* = nullptr>
inline auto
operator*(fvar<T> const& t, std::complex<U> const& u) {
  return stan::math::cplx::operator*(t,u);
}

template <class T, class U,
 std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
inline auto
operator*(std::complex<T> const& t, fvar<U> const& u) {
 return stan::math::cplx::operator*(t,u);
}

template <class T, class U,
 std::enable_if_t<std::is_arithmetic<U>::value>* = nullptr>
inline auto
operator/(fvar<T> const& t, std::complex<U> const& u) {
  return stan::math::cplx::operator/(t,u);
}

template <class T, class U,
 std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
inline auto
operator/(std::complex<T> const& t, fvar<U> const& u) {
 return stan::math::cplx::operator/(t,u);
}

}  // namespace math
}  // namespace stan
#endif
