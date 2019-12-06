#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal/fun/value_of_rec.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/scal/fun/atan2.hpp>
#include <stan/math/rev/scal/fun/cos.hpp>
#include <stan/math/rev/scal/fun/exp.hpp>
#include <stan/math/rev/scal/fun/hypot.hpp>
#include <stan/math/rev/scal/fun/is_inf.hpp>
#include <stan/math/rev/scal/fun/is_nan.hpp>
#include <stan/math/rev/scal/fun/log.hpp>
#include <stan/math/rev/scal/fun/sin.hpp>
#include <stan/math/rev/scal/fun/square.hpp>
#include <stan/math/rev/scal/fun/sqrt.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iterator>
#include <limits>

#include <iostream>

// BEGIN PR 1
// ===================================================================
// specialization of std:: and overloads
// -------------------------------------
// std::iterator_traits<var>         stan/math/rev/core/std_iterator_traits.hpp
// std::iterator_traits<fvar<T>>     stan/math/fwd/core/std_iterator_traits.hpp
// stan::math::signbit               stan/math/prim/scal/fun/signbit.hpp
// stan::math::isinf                 stan/math/prim/scal/fun/isinf.hpp
// stan::math::isfinite              stan/math/prim/scal/fun/isfinite.hpp
// stan::math::isnan                 stan/math/prim/scal/fun/isnan.hpp
// stan::math::isnormal              stan/math/prim/scal/fun/isnormal.hpp
// stan::math::copysign              stan/math/prim/scal/fun/copysign.hpp

// std SPECIALIZATIONS (no complex)
namespace std {
/**
 * Specialization of iterator traits for Stan math.  These all take
 * the form of typedefs.
 */
template <>
struct iterator_traits<stan::math::var> {
  /**
   * Iterator category for traits.
   */
  typedef random_access_iterator_tag iterator_category;

  /**
   * Type for difference between pointers.
   */
  typedef ptrdiff_t difference_type;

  /**
   * Type for value of pointer to values.
   */
  typedef stan::math::var value_type;

  /**
   * Type of pointer to variables.
   */
  typedef stan::math::var* pointer;

  /**
   * Type of reference to variables.
   */
  typedef stan::math::var& reference;
};
}  // namespace std

namespace std {
/**
 * Specialization of iterator traits for Stan math.  These all take
 * the form of typedefs.
 */
template <typename T>
struct iterator_traits<stan::math::fvar<T>> {
  /**
   * Iterator category for traits.
   */
  typedef random_access_iterator_tag iterator_category;

  /**
   * Type for difference between pointers.
   */
  typedef ptrdiff_t difference_type;

  /**
   * Type for value of pointer to values.
   */
  typedef stan::math::fvar<T> value_type;

  /**
   * Type of pointer to variables.
   */
  typedef stan::math::fvar<T>* pointer;

  /**
   * Type of reference to variables.
   */
  typedef stan::math::fvar<T>& reference;
};
}  // namespace std

namespace stan {
namespace math {
/**
 * Return `true` if the specified argument is negative and `false`
 * otherwise.
 *
 * Overloads `std::signbit` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return `true` if the argument is negative
 */
template <typename T, require_autodiff_t<T>...>
inline bool signbit(T&& v) {
  using std::signbit;
  return signbit(v.val());
}

/**
 * Return true if specified argument is infinite (positive or
 * negative).
 *
 * Overloads `std::isinf` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return true if argument is infinite
 */
template <typename T, require_autodiff_t<T>...>
inline bool isinf(T&& v) {
  using std::isinf;
  return isinf(v.val());
}

/**
 * Return true if specified argument is finite (not infinite and not
 * not-a-number).
 *
 * Overloads `std::isfinite` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return true if argument is finite
 */
template <typename T, require_autodiff_t<T>...>
inline bool isfinite(T&& v) {
  using std::isfinite;
  return isfinite(v.val());
}

/**
 * Return true if specified argument is not-a-number.
 *
 * Overloads `std::isnan` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return true if argument is not-a-number
 */
template <typename T, require_autodiff_t<T>...>
inline bool isnan(T&& v) {
  using std::isnan;
  return isnan(v.val());
}

/**
 * Return true if specified argument is normal.  A number is normal if
 * it is finite, non-zero and not subnormal.
 *
 * Overloads `std::isnormal` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam T type of argument
 * @param[in] v argument
 * @return true if argument is normal
 */
template <typename T, require_autodiff_t<T>...>
inline bool isnormal(T&& v) {
  using std::isnormal;
  return isnormal(v.val());
}

/**
 * Return the negation of the first argument if the first and second
 * argument have different signs, otherwise return a copy of the first
 * argument.  For the sake of this function, zero is considered
 * positive.  This function uses negation rather than literally copying
 * signs to preserve derivatives.
 *
 * Overload of `std::copysign` from `cmath` for argument-dependent
 * lookup.
 *
 * @tparam T type of first argument
 * @tparam U type of second argument
 * @param[in] x first complex argument
 * @param[in] y second complex argument
 * @return copy of second argument, negated if necessary to match sign
 * of first argument
 */
template <typename T, typename U, require_all_not_complex_t<T, U>...>
inline auto copysign(T&& x, U&& y) {
  // 0 is considered positive
  return (x < 0 && y >= 0) || (x > 0 && y < 0) ? -x : x;
}
}  // namespace math
}  // namespace stan

// PR 2
// =============================================================
// complex base class and specializations for autodiff
// ---------------------------------------------------
// stan::math::complex_base<V>  CRTP base class
// stan/math/prim/cplx/complex_base.hpp std::complex<var>
// specialization    stan/math/rev/cplx/std_complex.hpp std::coplex<fvar<T>>
// specialization    stan/math/fwd/cplx/std_complex.hpp ALT:  cplx -> core

namespace stan {
namespace math {
/**
 * Base class for complex numbers.  The template parameter `V` is the
 * value type, and the derived type for the curiously recursive
 * template pattern (CRTP) is `complex<V>`.  For this pattern to work,
 * the extending class must be of the from `complex<V>`.
 *
 * @tparam V value type for extending complex class
 */
template <typename V>
class complex_base {
 protected:
  /**
   * Real part.
   */
  V re_;
  /**
   * Imaginary part.
   */
  V im_;

 public:
  /**
   * Type of real and imaginary parts.
   */
  using value_type = V;

  /**
   * Derived complex type used for function return types.
   */
  using complex_type = std::complex<value_type>;

  /**
   * Construct a complex base with the specified real and complex
   * parts.
   *
   * @tparam T real type (must be assignable to `V`)
   * @tparam U imaginary type (must be assignable to `V`)
   * @param[in] x real part
   * @param[in] y imaginary part (default 0)
   */
  template <typename T, typename U>
  complex_base(const T& x, const U& y = U(0))  // NOLINT(runtime/explicit)
      : complex_base(value_type(x), value_type(y)) {}

  /**
   * Constructs complex number from real and imaginary parts.
   *
   * @param[in] re the real part
   * @param[in] im the imaginary part
   */
  complex_base(const value_type& re = value_type(0),
               const value_type& im = value_type(0))
      : re_(re), im_(im) {}

  /**
   * Constructs the complex number from the specified complex number
   * of a different type.
   *
   * @param[in] other another complex to use as source
   */
  template <typename T>
  complex_base(const std::complex<T>& other)  // NOLINT(runtime/explicit)
      : complex_base(other.real(), other.imag()) {}

  /**
   * Destroy this complex number.
   */
  ~complex_base() {}

  /**
   * Return a reference to thi class cast to the derived complex
   * type.
   *
   * @return reference to this class cast to the complex return type
   */
  complex_type& derived_complex_ref() {
    return static_cast<complex_type&>(*this);
  }

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @param[in] x value to assign
   * @return this complex number
   */
  template <typename T>
  complex_type& operator=(const T& x) {
    re_ = x;
    im_ = 0;
    return derived_complex_ref();
  }

  /**
   * Assign the real and imaginary parts of the specified complex
   * number to the real and imaginary part of this complex number.
   *
   * @tparam T value type of argument (must be assignable to
   * `value_type`)
   * @param[in] x complex value to assign
   * @return this complex number
   */
  template <typename T>
  complex_type& operator=(const std::complex<T>& x) {
    re_ = x.real();
    im_ = x.imag();
    return derived_complex_ref();
  }

  /**
   * Return the real part.
   *
   * @return the real part
   */
  value_type real() const { return re_; }

  /**
   * Set the real part to the specified value.
   *
   * @param[in] x the value to set the real part to
   */
  void real(const value_type& x) { re_ = x; }

  /**
   * Return the imaginary part.
   *
   * @return the imaginary part
   */
  value_type imag() const { return im_; }

  /**
   * Set the imaginary part to the specified value.
   *
   * @param[in] x the value to set the imaginary part to
   */
  void imag(const value_type& x) { im_ = x; }

  /**
   * Adds other to this.
   *
   * @tparam X type of scalar argument (assignable to value type)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator+=(const X& other) {
    re_ += other;
    return derived_complex_ref();
  }

  /**
   * Adds other to this.
   *
   * @tparam value type of complex argument (assignable to value type)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator+=(const std::complex<X>& other) {
    re_ += other.real();
    im_ += other.imag();
    return derived_complex_ref();
  }

  /**
   * Subtracts other from this.
   *
   * @tparam X type of scalar argument (assignable to vlaue type)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator-=(const X& other) {
    re_ -= other;
    return derived_complex_ref();
  }

  /**
   * Subtracts other from this.
   *
   * @tparam X value type of complex argument (assignable to value type)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator-=(const std::complex<X>& other) {
    re_ -= other.real();
    im_ -= other.imag();
    return derived_complex_ref();
  }

  /**
   * Multiplies this by other.
   *
   * @tparam X type of scalar argument (assignable to value type)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator*=(const X& other) {
    re_ *= other;
    im_ *= other;
    return derived_complex_ref();
  }

  /**
   * Multiplies this by other.
   *
   * @tparam X value type of complex argument (assignable to value type)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator*=(const std::complex<X>& other) {
    value_type re_temp = re_ * other.real() - im_ * other.imag();
    im_ = re_ * other.imag() + other.real() * im_;
    re_ = re_temp;
    return derived_complex_ref();
  }

  /**
   * Divides this by other.
   *
   * @tparam X type of scalar argument (assignable to value type)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator/=(const X& other) {
    re_ /= other;
    im_ /= other;
    return derived_complex_ref();
  }

  /**
   * Divides this by other.
   *
   * @tparam X value type of complex argument (assignable to value type)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename X>
  complex_type& operator/=(const std::complex<X>& other) {
    using stan::math::square;
    value_type sum_sq_im = square(other.real()) + square(other.imag());
    value_type re_temp = (re_ * other.real() + im_ * other.imag()) / sum_sq_im;
    im_ = (im_ * other.real() - re_ * other.imag()) / sum_sq_im;
    re_ = re_temp;
    return derived_complex_ref();
  }
};

}  // namespace math
}  // namespace stan

// SPECIALIZATION std::complex<var>
namespace std {

/**
 * Specialization of the standard libary complex number type for
 * reverse-mode autodiff type `stan::math::var`.
 */
template <>
class complex<stan::math::var>
    : public stan::math::complex_base<stan::math::var> {
 public:
  using base_t = complex_base<stan::math::var>;

  /**
   * Construct complex number from real and imaginary parts.
   *
   * @tparam V1 type of real part
   * @tparam V2 type of imaginary part
   * @param[in] re real part
   * @param[in] im imaginary part
   */
  template <typename V1, typename V2>
  complex(const V1& re, const V2& im) : complex_base(re, im) {}

  /**
   * Constructs complex number from real part or with default zero
   * value, setting imaginary part to zero.
   *
   * @param[in] re the real part
   */
  complex(const value_type& re = value_type(0))  // NOLINT(runtime/explicit)
      : complex_base(re) {}

  /**
   * Constructs the complex number from the specified complex number.
   *
   * @tparam V value type of complex argument
   * @param[in] other another complex to use as source
   */
  template <typename V>
  complex(const complex<V>& other) : complex_base(other) {}

  /**
   * Destroy this complex number.
   */
  ~complex() {}

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @tparam V type of value
   * @param[in] x value to assign
   * @return this complex number
   */
  template <typename V>
  complex_type& operator=(const V& x) {
    return base_t::operator=(x);
  }

  /**
   * Assign the real and imaginary parts of the specified complex
   * number to the real and imaginary part of this complex number.
   *
   * @tparam V value type of argument (assignable to `stan::math::var`)
   * @param[in] x complex value to assign
   * @return this complex number
   */
  template <typename V>
  complex_type& operator=(const complex<V>& x) {
    return base_t::operator=(x);
  }

  /**
   * Adds other to this.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator+=(const V& other) {
    return base_t::operator+=(other);
  }

  /**
   * Adds other to this.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator+=(const complex<V>& other) {
    return base_t::operator+=(other);
  }

  /**
   * Subtracts other from this.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator-=(const V& other) {
    return base_t::operator-=(other);
  }

  /**
   * Subtracts other from this.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator-=(const complex<V>& other) {
    return base_t::operator-=(other);
  }

  /**
   * Multiplies this by other.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator*=(const V& other) {
    return base_t::operator*=(other);
  }

  /**
   * Multiplies this by other.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator*=(const complex<V>& other) {
    return base_t::operator*=(other);
  }

  /**
   * Divides this by other.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::var`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator/=(const V& other) {
    return base_t::operator/=(other);
  }

  /**
   * Divides this by other.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::var`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator/=(const complex<V>& other) {
    return base_t::operator/=(other);
  }
};

}  // namespace std

// SPECIALIZATION std::complex<fvar<T>>
namespace std {

/**
 * Specialization of the standard libary complex number type for
 * reverse-mode autodiff type `stan::math::fvar<T>`.
 */
template <typename T>
class complex<stan::math::fvar<T>>
    : public stan::math::complex_base<stan::math::fvar<T>> {
 public:
  using base_t = stan::math::complex_base<stan::math::fvar<T>>;
  using value_type = stan::math::fvar<T>;
  using complex_type = complex<value_type>;

  /**
   * Construct complex number from real and imaginary parts.
   *
   * @tparam V1 type of real part
   * @tparam V2 type of imaginary part
   * @param[in] re real part
   * @param[in] im imaginary part
   */
  template <typename V1, typename V2>
  complex(const V1& re, const V2& im)
      : stan::math::complex_base<stan::math::fvar<T>>(re, im) {}

  /**
   * Constructs complex number from real part or with default zero
   * value, setting imaginary part to zero.
   *
   * @param[in] re the real part
   */
  complex(const value_type& re = value_type(0))  // NOLINT(runtime/explicit)
      : stan::math::complex_base<stan::math::fvar<T>>(re) {}

  /**
   * Constructs the complex number from the specified complex number.
   *
   * @tparam V value type of complex argument
   * @param[in] other another complex to use as source
   */
  template <typename V>
  complex(const complex<V>& other)
      : stan::math::complex_base<stan::math::fvar<T>>(other) {}

  /**
   * Destroy this complex number.
   */
  ~complex() {}

  /**
   * Assign the specified value to the real part of this complex number
   * and set imaginary part to zero.
   *
   * @tparam V type of value
   * @param[in] x value to assign
   * @return this complex number
   */
  template <typename V>
  complex_type& operator=(const V& x) {
    return base_t::operator=(x);
  }

  /**
   * Assign the real and imaginary parts of the specified complex
   * number to the real and imaginary part of this complex number.
   *
   * @tparam V value type of argument (assignable to `stan::math::fvar<T>`)
   * @param[in] x complex value to assign
   * @return this complex number
   */
  template <typename V>
  complex_type& operator=(const complex<V>& x) {
    return base_t::operator=(x);
  }

  /**
   * Adds other to this.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::fvar<T>`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator+=(const V& other) {
    return base_t::operator+=(other);
  }

  /**
   * Adds other to this.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::fvar<T>`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator+=(const complex<V>& other) {
    return base_t::operator+=(other);
  }

  /**
   * Subtracts other from this.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::fvar<T>`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator-=(const V& other) {
    return base_t::operator-=(other);
  }

  /**
   * Subtracts other from this.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::fvar<T>`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator-=(const complex<V>& other) {
    return base_t::operator-=(other);
  }

  /**
   * Multiplies this by other.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::fvar<T>`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator*=(const V& other) {
    return base_t::operator*=(other);
  }

  /**
   * Multiplies this by other.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::fvar<T>`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator*=(const complex<V>& other) {
    return base_t::operator*=(other);
  }

  /**
   * Divides this by other.
   *
   * @tparam V type of scalar argument (assignable to `stan::math::fvar<T>`)
   * @param[in] other a scalar value of matching type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator/=(const V& other) {
    return base_t::operator/=(other);
  }

  /**
   * Divides this by other.
   *
   * @tparam V value type of complex argument (assignable to
   * `stan::math::fvar<T>`)
   * @param[in] other a complex value of compatible type
   * @return this complex number
   */
  template <typename V>
  complex_type& operator/=(const complex<V>& other) {
    return base_t::operator/=(other);
  }
};

}  // namespace std

// BEGIN PR 3
// =======================================================
// generic std::complex code w/o complex<var> or complex<fvar<T>>

// stan::math::complex_return   stan/math/prim/cplx/meta/complex_return
// stan::math::X                stan/math/prim/cplx/fun/X.hpp
// stan::math::complex_Y        stan/math/prim/cplx/fun/complex_Y.hpp
//   for X = value_of_rec, copysign, i_times, neg_i_times
//   for Y = negate, add, ..., polar

// ALTERNATIVELY: cplx -> scal

namespace stan {
namespace math {
/**
 * Base template class for metaprogram determining type of return for
 * arithmetic operation involving two template types through defined
 * type `complex_t`.  Works for binary functions only.
 *
 * @tparam U type of first argument
 * @tparam V type of second argument
 */
template <typename U, typename V>
struct complex_return {};
template <typename U, typename V>
struct complex_return<std::complex<U>, std::complex<V>> {
  using complex_t = std::complex<return_type_t<U, V>>;
};
// doc inherited from base function template
template <typename U, typename V>
struct complex_return<std::complex<U>, V> {
  using complex_t = std::complex<return_type_t<U, V>>;
};
// doc inherited from base function template
template <typename U, typename V>
struct complex_return<U, std::complex<V>> {
  using complex_t = std::complex<return_type_t<U, V>>;
};

/**
 * The return type of complex function given a sequence of template
 * parameters for the argument types.
 *
 * @tparam Args argument types
 */
template <typename... Args>
using complex_return_t = typename complex_return<std::decay_t<Args>...>::complex_t;

/**
 * Return the complex number with real and complex parts given by the
 * `double` value of the real and complex parts of the argument.
 *
 * @tparam T value type of complex argument
 * @param[in] complex argument
 * @return complex number with same `double` real and imaginary parts as the
 * argument
 */
template <typename T>
inline std::complex<double> value_of_rec(const std::complex<T>& z) {
  using stan::math::value_of_rec;
  return {value_of_rec(z.real()), value_of_rec(z.imag())};
}

/**
 * Return the complex number composed of the real and complex parts
 * with signs copied from the real and complex parts of the first
 * arguments to the real and complex parts of the second.
 *
 * This is an overload of the standard libary `copysign` for complex
 * numbers that will be used with argument-dependent lookup.
 *
 * @tparam T value type of first argument
 * @tparam U value type of second argument
 * @param[in] x first complex argument
 * @param[in] y second complex argument
 * @return copy of second argument, with components negated if
 * necessary to match sign of first argument
 */
template <typename T, typename U>
inline std::complex<T> copysign(const std::complex<T>& y,
                                const std::complex<U>& x) {
  return {copysign(y.real(), x.real()), copysign(y.imag(), x.imag())};
}

/**
 * Return the specified complex number multiplied by `i`.
 *
 * This compound function is more efficient than mulitplying by a
 * constant `i` because it involves only a single arithmetic negation.
 *
 * @tparam value type of complex argument
 * @param[in] z complex argument
 * @return argument multipled by `i`
 */
template <typename T>
inline std::complex<T> i_times(const std::complex<T>& z) {
  return {-z.imag(), z.real()};
}

/**
 * Return the specified complex number multiplied by `-i`.
 *
 * This compound function is more efficient than mulitplying by the
 * constant `-i` because it involves only a single arithmetic
 * negation.
 *
 * @tparam value type of complex argument
 * @param[in] z complex argument
 * @return argument multipled by `-i`
 */
template <typename T>
inline std::complex<T> neg_i_times(const std::complex<T>& z) {
  return {z.imag(), -z.real()};
}

/**
 * Return the complex negation of the specified complex argument.
 *
 * @tparam V value type of complex argument
 * @param[in] z argument
 * @return negation of argument
 */
template <typename V>
inline std::complex<V> complex_negate(const std::complex<V>& z) {
  return {-z.real(), -z.imag()};
}

/**
 * Return the sum of the specified arguments.  At least one of the
 * arguments must be a complex number.
 *
 * @tparam U type of first argument
 * @tparam V type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return sum of the arguments
 */
template <typename U, typename V, require_any_complex_t<U, V>...>
inline complex_return_t<U, V> add(U&& lhs, V&& rhs) {
  complex_return_t<U, V> y(std::forward<U>(lhs));
  y += rhs;
  return y;
}

/**
 * Return the difference of the specified arguments.  At least one of
 * the arguments must be a complex number.
 *
 * @tparam U type of first argument
 * @tparam V type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return difference of the arguments
 */
template <typename U, typename V, require_any_complex_t<U, V>...>
inline complex_return_t<U, V> subtract(U&& lhs, V&& rhs) {
  complex_return_t<U, V> y(std::forward<U>(lhs));
  y -= rhs;
  return y;
}

/**
 * Return the product of the specified arguments.  At least one of the
 * arguments must be a complex number.
 *
 * @tparam U type of first argument
 * @tparam V type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return sum of the arguments
 */
template <typename U, typename V, require_any_complex_t<U, V>...>
inline complex_return_t<U, V> multiply(U&& lhs, V&& rhs) {
  complex_return_t<U, V> y(std::forward<U>(lhs));
  y *= rhs;
  return y;
}

/**
 * Return the quotient of the specified arguments.  At least one of
 * the arguments must be a complex number.
 *
 * @tparam U type of first argument
 * @tparam V type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return sum of the arguments
 */
template <typename U, typename V, require_any_complex_t<U, V>...>
inline complex_return_t<U, V> divide(U&& lhs, V&& rhs) {
  complex_return_t<U, V> y(std::forward<U>(lhs));
  y /= rhs;
  return y;
}

/**
 * Return `true` if the complex arguments have equal real and
 * imaginary components.
 *
 * @tparam U value type of first argument
 * @tparam V value type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return `true` if the arguments are equal
 */
template <typename U, typename V, require_all_complex_t<U, V>...>
inline bool equal_equal(U&& lhs, V&& rhs) {
  return lhs.real() == rhs.real() && lhs.imag() == rhs.imag();
}

/**
 * Return `true` if the complex argument has a real component equal to
 * the real argument.
 *
 * @tparam U value type of first argument
 * @tparam V type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return `true` if the complex argument has real component equal to
 * the real argument
 */
template <typename U, typename V, require_not_complex_t<U>..., require_complex_t<V>...>
inline bool equal_equal(U&& lhs, V&& rhs) {
  return lhs == rhs.real() && rhs.imag() == 0;
}

/**
 * Return `true` if the complex argument has a real component equal to
 * the real argument.
 *
 * @tparam U value type of first argument
 * @tparam V type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return `true` if the complex argument has real component equal to
 * the real argument
 */
template <typename U, typename V, require_not_complex_t<V>..., require_complex_t<U>...>
inline bool equal_equal(U&& lhs, V&& rhs) {
  return lhs.real() == rhs && lhs.imag() == 0;
}

/**
 * Return `true` if the arguments are not equal.  At least one of the
 * arguments must be a complex number.
 *
 * @tparam U type of first argument
 * @tparam V type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return `true` if the arguments are not equal
 */
template <typename U, typename V, require_any_complex_t<U, V>...>
inline bool not_equal(U&& lhs, V&& rhs) {
  return !equal_equal(std::forward<U>(lhs), std::forward<V>(rhs));
}

/**
 * Return the real part of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return real part of argument
 */
template <typename V, require_complex_t<V>...>
inline auto real(V&& z) {
  return z.real();
}

/**
 * Return the imaginary part of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return imaginary part of argument
 */
template <typename V, require_complex_t<V>...>
inline auto imag(V&& z) {
  return z.imag();
}

/**
 * Return the absolute value of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return absolute value of the argument
 */
template <typename V, require_complex_t<V>...>
inline auto complex_abs(V&& z) {
  return hypot(real(z), imag(z));
}

/**
 * Return the phase angle of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return phase angle of the argument
 */
template <typename V, require_complex_t<V>...>
inline auto complex_arg(V&& z) {
  return atan2(imag(z), real(z));
}

/**
 * Return the squared magnitude of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return squared magnitude of the argument
 */
template <typename V, require_complex_t<V>...>
inline auto complex_norm(V&& z) {
  return square(real(z)) + square(imag(z));
}

/**
 * Return the complex conjugate the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return complex conjugate of the argument
 */
template <typename V, require_complex_t<V>...>
inline auto complex_conj(V&& z) {
  return std::complex<value_type_t<V>>(z.real(), -z.imag());
}

/**
 * Return the projection of the complex argument onto the Riemann
 * sphere.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return projection of the argument onto the Riemann sphere
 */
template <typename V>
inline std::complex<V> complex_proj(const std::complex<V>& z) {
  if (is_inf(z.real()) || is_inf(z.imag())) {
    return {std::numeric_limits<V>::infinity(), z.imag() < 0 ? -0.0 : 0.0};
  }
  return z;
}

/**
 * Return the natural exponent of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return exponential of the argument
 */
template <typename V>
inline std::complex<V> complex_exp(const std::complex<V>& z) {
  if (is_inf(z.real()) && z.real() > 0) {
    if (is_nan(z.imag()) || z.imag() == 0) {
      // (+inf, nan), (+inf, 0)
      return z;
    } else if (is_inf(z.imag()) && z.imag() > 0) {
      // (+inf, +inf)
      return {z.real(), std::numeric_limits<double>::quiet_NaN()};
    } else if (is_inf(z.imag()) && z.imag() < 0) {
      // (+inf, -inf)
      return {std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN()};
    }
  }
  if (is_inf(z.real()) && z.real() < 0
      && (is_nan(z.imag()) || is_inf(z.imag()))) {
    // (-inf, nan), (-inf, -inf), (-inf, inf)
    return {0, 0};
  }
  if (is_nan(z.real()) && z.imag() == -0.0) {
    // (nan, -0)
    return z;
  }
  V exp_re = exp(z.real());
  return {exp_re * cos(z.imag()), exp_re * sin(z.imag())};
}

/**
 * Return the natural logarithm of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return natural logarithm of the argument
 */
template <typename V>
inline std::complex<V> complex_log(const std::complex<V>& z) {
  static const double inf = std::numeric_limits<double>::infinity();
  static const double nan = std::numeric_limits<double>::quiet_NaN();
  if ((is_nan(z.real()) && is_inf(z.imag()))
      || (is_inf(z.real()) && is_nan(z.imag()))) {
    return {inf, nan};
  }
  V r = sqrt(norm(z));
  V theta = arg(z);
  return {log(r), theta};
}

/**
 * Return the base 10 logarithm of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return base 10 logarithm of the argument
 */
template <typename V>
inline std::complex<V> complex_log10(const std::complex<V>& z) {
  static const double inv_log_10 = 1 / log(10);
  return log(z) * V(inv_log_10);
}

/**
 * Return the first argument raised to the power of the second
 * argument.  At least one of the arguments must be a complex number.
 *
 * @tparam U type of first argument
 * @tparam V type of second argument
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return first argument raised to the power of the second argument
 */
template <typename U, typename V>
inline complex_return_t<U, V> complex_pow(const U& x, const V& y) {
  return exp(y * log(x));
}

/**
 * Return the square root of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return square root of the argument
 */
template <typename V>
inline std::complex<V> complex_sqrt(const std::complex<V>& z) {
  auto m = sqrt(hypot(z.real(), z.imag()));
  auto at = 0.5 * atan2(z.imag(), z.real());
  return {m * cos(at), m * sin(at)};
}

/**
 * Return the hyperbolic sine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic sine of the argument
 */
template <typename V>
inline std::complex<V> complex_sinh(const std::complex<V>& z) {
  return 0.5 * (exp(z) - exp(-z));
}

/**
 * Return the hyperbolic cosine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic cosine of the argument
 */
template <typename V>
inline std::complex<V> complex_cosh(const std::complex<V>& z) {
  return 0.5 * (exp(z) + exp(-z));
}

/**
 * Return the hyperbolic tangent of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic tangent of the argument
 */
template <typename V>
inline auto complex_tanh(const std::complex<V>& z) {
  auto exp_z = exp(z);
  auto exp_neg_z = exp(-z);
  return (exp_z - exp_neg_z) / (exp_z + exp_neg_z);
}

/**
 * Return the hyperbolic arc sine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic arc sine of the argument
 */
template <typename V>
inline std::complex<V> complex_asinh(const std::complex<V>& z) {
  std::complex<double> y_d = asinh(value_of_rec(z));
  auto y = log(z + sqrt(1 + z * z));
  return copysign(y, y_d);
}

/**
 * Return the hyperbolic arc cosine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic arc cosine of the argument
 */
template <typename V>
inline std::complex<V> complex_acosh(const std::complex<V>& z) {
  std::complex<double> y_d = acosh(value_of_rec(z));
  auto y = log(z + sqrt(z * z - 1));
  return copysign(y, y_d);
}

/**
 * Return the hyperbolic arc tangent of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic arc tangent of the argument
 */
template <typename V>
inline std::complex<V> complex_atanh(const std::complex<V>& z) {
  std::complex<double> y_d = atanh(value_of_rec(z));
  V one(1);
  auto y = 0.5 * (log(one + z) - log(one - z));
  return copysign(y, y_d);
}

/**
 * Return the sine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return sine of the argument
 */
template <typename V>
inline std::complex<V> complex_sin(const std::complex<V>& z) {
  return neg_i_times(sinh(i_times(z)));
}

/**
 * Return the cosine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return cosine of the argument
 */
template <typename V>
inline std::complex<V> complex_cos(const std::complex<V>& z) {
  return cosh(i_times(z));
}

/**
 * Return the tangent of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return tangent of the argument
 */
template <typename V>
inline std::complex<V> complex_tan(const std::complex<V>& z) {
  return neg_i_times(tanh(i_times(z)));
}

/**
 * Return the arc sine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return arc sine of the argument
 */
template <typename V>
inline std::complex<V> complex_asin(const std::complex<V>& z) {
  auto y_d = asin(value_of_rec(z));
  auto y = neg_i_times(asinh(i_times(z)));
  return copysign(y, y_d);
}

/**
 * Return the arc cosine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return arc cosine of the argument
 */
template <typename V>
inline std::complex<V> complex_acos(const std::complex<V>& z) {
  return V(0.5 * pi()) - asin(z);
}

/**
 * Return the arc tangent of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return arc tangent of the argument
 */
template <typename V>
inline std::complex<V> complex_atan(const std::complex<V>& z) {
  return neg_i_times(atanh(i_times(z)));
}

/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
template <typename U, typename V>
inline std::complex<return_type_t<U, V>> complex_polar(const U& r,
                                                       const V& theta) {
  using std::cos;
  using std::sin;
  if (!(r >= 0) || is_inf(theta)) {
    return {std::numeric_limits<double>::quiet_NaN()};
  }
  return {r * cos(theta), r * sin(theta)};
}

}  // namespace math
}  // namespace stan

namespace std {
// Function template specializations required to get around
// non-ADL-friendly coding in libstdc++ for g++.  For background, see:
// Walter E. Brown.  2017.  Thou Shalt Not Specialize std Function Templates!
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2017/p0551r1.pdf
//

/**
 * Return true if the respective parts of the two arguments are equal
 * and false otherwise.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second complex argument
 * @return if the arguments are equal
 */
template <>
inline bool operator==<stan::math::var>(const complex<stan::math::var>& lhs,
                                        const complex<stan::math::var>& rhs) {
  return equal_equal(lhs, rhs);
}

/**
 * Return true if the respective parts of the two arguments are equal
 * and false otherwise.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second scalar argument
 * @return if the arguments are equal
 */
template <>
inline bool operator==<stan::math::var>(const complex<stan::math::var>& lhs,
                                        const stan::math::var& rhs) {
  return equal_equal(lhs, rhs);
}

/**
 * Return true if the respective parts of the two arguments are equal
 * and false otherwise.
 *
 * @param[in] lhs first scalar argument
 * @parma[in] rhs second complex argument
 * @return if the arguments are equal
 */
template <>
inline bool operator==<stan::math::var>(const stan::math::var& lhs,
                                        const complex<stan::math::var>& rhs) {
  return equal_equal(lhs, rhs);
}

/**
 * Return true if the respective parts of the two arguments are not
 * equal and false otherwise.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second complex argument
 * @return if the arguments are not equal
 */
template <>
inline bool operator!=<stan::math::var>(const complex<stan::math::var>& lhs,
                                        const complex<stan::math::var>& rhs) {
  return not_equal(lhs, rhs);
}

/**
 * Return true if the respective parts of the two arguments are not
 * equal and false otherwise.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second scalar argument
 * @return if the arguments are not equal
 */
template <>
inline bool operator!=<stan::math::var>(const complex<stan::math::var>& lhs,
                                        const stan::math::var& rhs) {
  return not_equal(lhs, rhs);
}

/**
 * Return true if the respective parts of the two arguments are not
 * equal and false otherwise.
 *
 * @param[in] lhs first scalar argument
 * @parma[in] rhs second complex argument
 * @return if the arguments are not equal
 */
template <>
inline bool operator!=<stan::math::var>(const stan::math::var& lhs,
                                        const complex<stan::math::var>& rhs) {
  return not_equal(lhs, rhs);
}

/**
 * Return the complex arc hyperbolic sine of the specified complex
 * argument with branch cuts outside the interval `[-i, i]` along the
 * imaginary axis.
 *
 * @param[in] z complex argument
 * @return arc hyperbolic sine of argument
 */
template <>
inline complex<stan::math::var> asinh<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_asinh(z);
}

/**
 * Return the complex arc hyperbolic cosine of the specified complex
 * argument with branch cuts at values less than 1 along the real
 * axis.
 *
 * @param[in] z complex argument
 * @return arc hyperbolic cosine of argument
 */
template <>
inline complex<stan::math::var> acosh<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_acosh(z);
}

/**
 * Return the complex arc hyperbolic tangent of the specified complex
 * argument with branch cuts outside the interval `[-i, i]` along the
 * real axis.
 *
 * @param[in] z complex argument
 * @return arc hyperbolic tanget of argument
 */
template <>
inline complex<stan::math::var> atanh<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_atanh(z);
}

/**
 * Return the value of the specified argument.
 *
 * @param[in] val the complex number argument
 * @return a copy of the argument
 */
template <>
inline complex<stan::math::var> operator+<stan::math::var>(
    const complex<stan::math::var>& val) {
  return val;
}

/**
 * Return the negation of the specified argument.
 *
 * @param[in] val the complex number argument
 * @return negated argument
 */
template <>
inline complex<stan::math::var> operator-<stan::math::var>(
    const complex<stan::math::var>& val) {
  return stan::math::complex_negate(val);
}

/**
 * Returns the sum of its arguments.
 *
 * @param[in] lhs complex first argument
 * @param[in] rhs complex second argument
 * @return sum of the arguments
 */
template <>
inline complex<stan::math::var> operator+<stan::math::var>(
    const complex<stan::math::var>& lhs, const complex<stan::math::var>& rhs) {
  return add(lhs, rhs);
}

/**
 * Returns the sum of its arguments.
 *
 * @param[in] lhs complex first argument
 * @param[in] rhs scalar second argument
 * @return sum of the arguments
 */
template <>
inline complex<stan::math::var> operator+<stan::math::var>(
    const complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  return add(lhs, rhs);
}

/**
 * Return the sum of the two arguments.
 *
 * @param[in] lhs scalar first argument
 * @param[in] rhs complex second argument
 * @return sum of the arguments
 */
template <>
inline complex<stan::math::var> operator+<stan::math::var>(
    const stan::math::var& lhs, const complex<stan::math::var>& rhs) {
  return add(lhs, rhs);
}

/**
 * Returns the result of subtracting the second argument from the
 * first.
 *
 * @param[in] lhs complex first argument
 * @param[in] rhs complex second argument
 * @return difference between first and second argument
 */
template <>
inline complex<stan::math::var> operator-<stan::math::var>(
    const complex<stan::math::var>& lhs, const complex<stan::math::var>& rhs) {
  return subtract(lhs, rhs);
}

/**
 * Returns the result of subtracting the second argument from the
 * first.
 *
 * @param[in] lhs complex first argument
 * @param[in] rhs scalar second argument
 * @return difference between first and second argument
 */
template <>
inline complex<stan::math::var> operator-<stan::math::var>(
    const complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  return subtract(lhs, rhs);
}

/**
 * Returns the result of subtracting the second argument from the
 * first.
 *
 * @param[in] lhs scalar first argument
 * @param[in] rhs complex second argument
 * @return difference between first and second argument
 */
template <>
inline complex<stan::math::var> operator-<stan::math::var>(
    const stan::math::var& lhs, const complex<stan::math::var>& rhs) {
  return subtract(lhs, rhs);
}

/**
 * Return the product of the first and second argument.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second complex argument
 * @return product of the first and second argument
 */
template <>
inline complex<stan::math::var> operator*<stan::math::var>(
    const complex<stan::math::var>& lhs, const complex<stan::math::var>& rhs) {
  return multiply(lhs, rhs);
}

/**
 * Return the product of the first and second argument.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second scalar argument
 * @return product of the first and second argument
 */
template <>
inline complex<stan::math::var> operator*<stan::math::var>(
    const complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  return multiply(lhs, rhs);
}

/**
 * Return the product of the first and second argument.
 *
 * @param[in] lhs first scalar argument
 * @parma[in] rhs second complex argument
 * @return product of the first and second argument
 */
template <>
inline complex<stan::math::var> operator*<stan::math::var>(
    const stan::math::var& lhs, const complex<stan::math::var>& rhs) {
  return multiply(lhs, rhs);
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second complex argument
 * @return quotient of the first and second argument
 */
template <>
inline complex<stan::math::var> operator/<stan::math::var>(
    const complex<stan::math::var>& lhs, const complex<stan::math::var>& rhs) {
  return divide(lhs, rhs);
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * @param[in] lhs first complex argument
 * @parma[in] rhs second scalar argument
 * @return quotient of the first and second argument
 */
template <>
inline complex<stan::math::var> operator/<stan::math::var>(
    const complex<stan::math::var>& lhs, const stan::math::var& rhs) {
  return divide(lhs, rhs);
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * @param[in] lhs first scalar argument
 * @parma[in] rhs second complex argument
 * @return quotient of the first and second argument
 */
template <>
inline complex<stan::math::var> operator/(const stan::math::var& lhs,
                                          const complex<stan::math::var>& rhs) {
  return divide(lhs, rhs);
}

/**
 * Return the real component of the specified complex number.
 *
 * @param[in] z complex argument
 * @return real component of argment
 */
template <>
inline stan::math::var real<stan::math::var>(
    const complex<stan::math::var>& z) {
  return stan::math::real(z);
}

/**
 * Return the imaginary component of the specified complex number.
 *
 * @param[in] z complex argument
 * @return real component of argment
 */
template <>
inline stan::math::var imag<stan::math::var>(
    const complex<stan::math::var>& z) {
  return stan::math::imag(z);
}

/**
 * Return the magnitude of the specified complex number.
 *
 * @param[in] complex argument
 * @return magnitude of argument
 */
template <>
inline stan::math::var abs<stan::math::var>(const complex<stan::math::var>& z) {
  return complex_abs(z);
}

/**
 * Return the phase angle of the specified complex number.
 *
 * @param[in] z complex argument
 * @return phase angle of argument
 */
template <>
inline stan::math::var arg<stan::math::var>(const complex<stan::math::var>& z) {
  return complex_arg(z);
}

/**
 * Return the squared magnitude of the specified complex number.
 *
 * @param[in] z complex argument
 * @return squared magnitude of argument
 */
template <>
inline stan::math::var norm<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_norm(z);
}

/**
 * Return the complex conjugate of the specified complex number.
 *
 * @param[in] z complex argument
 * @return complex conjugate of argument
 */
template <>
inline complex<stan::math::var> conj<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_conj(z);
}

/**
 * Return the projection of the specified complex arguent onto the
 * Riemann sphere.
 *
 * @param[in] z complex argument
 * @return project of argument onto the Riemann sphere
 */
template <>
inline complex<stan::math::var> proj<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_proj(z);
}

/**
 * Returns complex number with specified magnitude and phase angle.
 *
 * @param[in] r magnitude
 * @param[in] theta phase angle
 * @return complex number with magnitude and phase angle
 */
template <>
inline complex<stan::math::var> polar<stan::math::var>(
    const stan::math::var& r, const stan::math::var& theta) {
  return complex_polar(r, theta);
}

/**
 * Return the natural exponent of the specified complex number.
 *
 * @param[in] complex argument
 * @return exponential of argument
 */
template <>
inline complex<stan::math::var> exp<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_exp(z);
}

/**
 * Return the natural logarithm of the specified complex value with a
 * branch cut along the negative real axis.
 *
 * @param[in] z complex argument
 * @return natural logarithm of argument
 */
template <>
inline complex<stan::math::var> log<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_log(z);
}

/**
 * Return the base 10 logarithm of the specified complex value with a
 * branch cut along the negative real axis.
 *
 * @param[in] z complex argument
 * @return base 10 logarithm of argument
 */
template <>
inline complex<stan::math::var> log10<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_log10(z);
}

template <>
inline complex<stan::math::var> pow<stan::math::var>(
    const complex<stan::math::var>& x, const complex<stan::math::var>& y) {
  return complex_pow(x, y);
}

// Following std:: templates cannot be specialized in both g++ and clang++

// can be specialized in g++, but not in clang++
// complex<var> pow(const complex<var>& x, const var& y);
// complex<var> pow(const var& x, const complex<var>& y);
// complex<var> pow(const complex<var>& x, int y);

// can't be specialized in g++ or clang++
// complex<var> pow(const complex<var>& x, const complex<double>& y);
// complex<var> pow(const complex<double>& x, const complex<var>& y);
// complex<var> pow(const complex<var>& x, const double& y);
// complex<var> pow(const double& x, const complex<var>& y);

/**
 * Return the square root of the specified complex number with a
 * branch cut along the negative real axis.
 *
 * @param[in] z complex number argument
 * @return square root of argument
 */
template <>
inline complex<stan::math::var> sqrt<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_sqrt(z);
}

/**
 * Return the complex hyperbolic sine of the specified complex
 * argument.
 *
 * @param[in] z complex argument
 * @return hyperbolic sine of argument
 */
template <>
inline complex<stan::math::var> sinh<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_sinh(z);
}

/**
 * Return the complex hyperbolic cosine of the specified complex
 * argument.
 *
 * @param[in] z complex argument
 * @return hyperbolic cosine of argument
 */
template <>
inline complex<stan::math::var> cosh<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_cosh(z);
}

/**
 * Return the complex hyperbolic tangent of the specified complex
 * argument.
 *
 * @param[in] z complex argument
 * @return hyperbolic tangent of argument
 */
template <>
inline complex<stan::math::var> tanh<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_tanh(z);
}

/**
 * Return the complex sine of the specified complex number.
 *
 * @param[in] z a complex argument
 * @return complex sine of argument
 */
template <>
inline complex<stan::math::var> sin<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_sin(z);
}

/**
 * Return the complex cosine of the specified complex number.
 *
 * @param[in] z a complex argument
 * @return complex cosine of argument
 */
template <>
inline complex<stan::math::var> cos<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_cos(z);
}

/**
 * Return the complex tangent of the specified complex number.
 *
 * @param[in] z a complex argument
 * @return complex tangent of argument
 */
template <>
inline complex<stan::math::var> tan<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_tan(z);
}

/**
 * Return the complex arc sine of the specified complex argument, with
 * branch cuts outside the interval `[-1, 1]` along the real axis.
 *
 * @param[in] z complex argument
 * @return complex arc sine of the argument
 */
template <>
inline complex<stan::math::var> asin<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_asin(z);
}

/**
 * Return the complex arc cosine of the specified complex argument, with
 * branch cuts outside the interval `[-1, 1]` along the real axis.
 *
 * @param[in] z complex argument
 * @return complex arc cosine of the argument
 */
template <>
inline complex<stan::math::var> acos<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_acos(z);
}

/**
 * Return the complex arc tangent of the specified complex argument,
 * with branch cuts outside the interval `[-i, i]` along the imaginary
 * axis.
 *
 * @param[in] z complex argument
 * @return complex arc tangent of the argument
 */
template <>
inline complex<stan::math::var> atan<stan::math::var>(
    const complex<stan::math::var>& z) {
  return complex_atan(z);
}
}  // namespace std

// OVERLOADS FOR ADL IN STAN
namespace stan {
namespace math {

/**
 * Writes specified complex number to specified ostream in
 * form `(real, imag)`.
 *
 * @tparam T complex number type
 * @tparam CharT character type for output stream
 * @tparam Traits traits for output stream
 * @param[in,out] o character output stream
 * @param[in] x complex number to be inserted
 * @return specified output stream
 */
template <typename T, class CharT, class Traits>
inline std::basic_ostream<CharT, Traits>& write_complex(
    std::basic_ostream<CharT, Traits>& o, const std::complex<T>& x) {
  std::basic_ostringstream<CharT, Traits> s;
  s.flags(o.flags());
  s.imbue(o.getloc());
  s.precision(o.precision());
  s << '(' << x.real() << "," << x.imag() << ')';
  return o << s.str();
}

/**
 * Reads complex number from specified input stream.
 * Supported formats are `real`, `(real)`, and `(real,imag)`.
 *
 * @tparam T complex number type
 * @tparam CharT character type for input stream
 * @tparam Traits traits for input stream
 * @param[in,out] is input stream from which to read
 * @param[in,out] x complex number to write
 * @return specified input stream
 */
template <typename T, class CharT, class Traits>
inline std::basic_istream<CharT, Traits>& read_complex(
    std::basic_istream<CharT, Traits>& is, std::complex<T>& x) {
  std::complex<double> y;
  is >> y;
  x = y;
  return is;
}

/**
 * Writes specified complex number to specified ostream in
 * form `(real, imag)`.
 *
 * @tparam CharT character type for output stream
 * @tparam Traits traits for output stream
 * @param[in,out] o character output stream
 * @param[in] xcomplex number to be inserted
 * @return specified output stream
 */
template <class CharT, class Traits>
inline std::basic_ostream<CharT, Traits>& operator<<(
    std::basic_ostream<CharT, Traits>& o, const std::complex<var>& x) {
  return write_complex(o, x);
}

/**
 * Read a complex number from the specified input stream.  The
 * supported formats are `real`, `(real)`, and `(real, imaginary)`.
 *
 * @tparam T forward mode value type
 * @tparam CharT character type for output stream
 * @tparam Traits traits for input stream
 * @param[in,out] is input stream from which to read
 * @param[in,out] x complex number to write
 * @return specified input stream
 */
template <class CharT, class Traits>
inline std::basic_istream<CharT, Traits>& operator>>(
    std::basic_istream<CharT, Traits>& is, std::complex<var>& x) {
  return read_complex(is, x);
}

/**
 * Writes specified complex number to specified ostream in
 * form `(real, imag)`.
 *
 * @tparam T forward mode value type
 * @tparam CharT character type for output stream
 * @tparam Traits traits for output stream
 * @param[in,out] o character output stream
 * @param[in] x complex number to be inserted
 * @return specified output stream
 */
template <typename T, class CharT, class Traits>
inline std::basic_ostream<CharT, Traits>& operator<<(
    std::basic_ostream<CharT, Traits>& o, const std::complex<fvar<T>>& x) {
  write_complex(o, x);
}

/**
 * Read a complex number from the specified input stream.  The
 * supported formats are `real`, `(real)`, and `(real, imaginary)`.
 *
 * @tparam T forward mode value type
 * @tparam CharT character type for output stream
 * @tparam Traits traits for input stream
 * @param[in,out] is input stream from which to read
 * @param[in,out] x complex number to write
 * @return specified input stream
 */
template <typename T, class CharT, class Traits>
inline std::basic_istream<CharT, Traits>& operator>>(
    std::basic_istream<CharT, Traits>& is, std::complex<fvar<T>>& x) {
  return read_complex(is, x);
}

// {int, double, var, complex<double>, complex<var>}
// with at least one var, at least one complex
//
// built-in template functions
//  A. complex<T>, complex<T>
//  B. complex<T>, T
//  C. T, complex<T>
//
// instantiations for autodiff type T
// never mix autodiff types other than with double
//  1. (int, complex<T>)
//  2. (double, complex<T>)
//  3. (T, complex<double>)
//  4. (T, complex<T>)
//  5. (complex<double>, T)
//  6. (complex<double>, complex<T>)
//  7. (complex<T>, int)
//  8. (complex<T>, double)
//  9. (complex<T>, T)
// 10. (complex<T>, complex<double>)
// 11. (complex<T>, complex<T>)

std::complex<var> operator+(const std::complex<var>& z) { return z; }
template <typename T>
inline std::complex<fvar<T>> operator+(const std::complex<var>& z) {
  return z;
}

std::complex<var> operator-(const std::complex<var>& z) {
  return complex_negate(z);
}
template <typename T>
inline std::complex<fvar<T>> operator-(const std::complex<var>& z) {
  return complex_negate(z);
}

/**
 * Return the sum of the two arguments.
 *
 * @tparam U value type of first argument
 * @tparam V value type of second argument
 * @param[in] x first complex argument
 * @param[in] y second complex argument
 * @return sum of the arguments
 */
template <typename U, typename V, require_any_complex_t<U, V>...>
inline auto operator+(U&& x, V&& y) {
  return add(std::forward<U>(x), std::forward<V>(y));
}

/**
 * Return the difference of the two arguments.
 *
 * @tparam U value type of first argument
 * @tparam V value type of second argument
 * @param[in] x first complex argument
 * @param[in] y second complex argument
 * @return difference between first and second argument
 */
template <typename U, typename V, require_any_complex_t<U, V>...>
inline auto operator-(U&& x, V&& y) {
  return subtract(std::forward<U>(x), std::forward<V>(y));
}

/**
 * Return the product of the two arguments.
 *
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return product of the arguments
 */
// var 1
std::complex<var> operator*(int lhs, const std::complex<var>& rhs) {
  return multiply(lhs, rhs);
}
// var 2
std::complex<var> operator*(double lhs, const std::complex<var>& rhs) {
  return multiply(lhs, rhs);
}
// var 3
std::complex<var> operator*(const var& lhs, const std::complex<double>& rhs) {
  return multiply(lhs, rhs);
}
// var 4
std::complex<var> operator*(const var& lhs, const std::complex<var>& rhs) {
  return multiply(lhs, rhs);
}
// var 5
std::complex<var> operator*(const std::complex<double>& lhs, const var& rhs) {
  return multiply(lhs, rhs);
}
// var 6
std::complex<var> operator*(const std::complex<double>& lhs,
                            const std::complex<var>& rhs) {
  return multiply(lhs, rhs);
}
// var 7
std::complex<var> operator*(const std::complex<var>& lhs, int rhs) {
  return multiply(lhs, rhs);
}
// var 8
std::complex<var> operator*(const std::complex<var>& lhs, double rhs) {
  return multiply(lhs, rhs);
}
// var 9
std::complex<var> operator*(const std::complex<var>& lhs, const var& rhs) {
  return multiply(lhs, rhs);
}
// var 10
std::complex<var> operator*(const std::complex<var>& lhs,
                            const std::complex<double>& rhs) {
  return multiply(lhs, rhs);
}
// var 11
std::complex<var> operator*(const std::complex<var>& lhs,
                            const std::complex<var>& rhs) {
  return multiply(lhs, rhs);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> operator*(int lhs,
                                       const std::complex<fvar<T>>& rhs) {
  return multiply(lhs, rhs);
}
// fvar 2
template <typename T>
inline std::complex<fvar<T>> operator*(double lhs,
                                       const std::complex<fvar<T>>& rhs) {
  return multiply(lhs, rhs);
}
// fvar 3
template <typename T>
inline std::complex<fvar<T>> operator*(const fvar<T>& lhs,
                                       const std::complex<double>& rhs) {
  return multiply(lhs, rhs);
}
// fvar 4
template <typename T>
inline std::complex<fvar<T>> operator*(const fvar<T>& lhs,
                                       const std::complex<fvar<T>>& rhs) {
  return multiply(lhs, rhs);
}
// fvar 5
template <typename T>
inline std::complex<fvar<T>> operator*(const std::complex<double>& lhs,
                                       const fvar<T>& rhs) {
  return multiply(lhs, rhs);
}
// fvar 6
template <typename T>
inline std::complex<fvar<T>> operator*(const std::complex<double>& lhs,
                                       const std::complex<fvar<T>>& rhs) {
  return multiply(lhs, rhs);
}
// fvar 7
template <typename T>
inline std::complex<fvar<T>> operator*(const std::complex<fvar<T>>& lhs,
                                       int rhs) {
  return multiply(lhs, rhs);
}
// fvar 8
template <typename T>
inline std::complex<fvar<T>> operator*(const std::complex<fvar<T>>& lhs,
                                       double rhs) {
  return multiply(lhs, rhs);
}
// fvar 9
template <typename T>
inline std::complex<fvar<T>> operator*(const std::complex<fvar<T>>& lhs,
                                       const fvar<T>& rhs) {
  return multiply(lhs, rhs);
}
// fvar 10
template <typename T>
inline std::complex<fvar<T>> operator*(const std::complex<fvar<T>>& lhs,
                                       const std::complex<double>& rhs) {
  return multiply(lhs, rhs);
}
// fvar 11
template <typename T>
inline std::complex<fvar<T>> operator*(const std::complex<fvar<T>>& lhs,
                                       const std::complex<fvar<T>>& rhs) {
  return multiply(lhs, rhs);
}

/**
 * Return the quotient of the two arguments.
 *
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return quotient of the arguments
 */
// 1
template <typename U, typename V, require_any_complex_t<U, V>...,
 require_any_autodiff_t<value_type_t<U>, value_type_t<V>>...>
inline auto operator/(U&& lhs, V&& rhs) {
  return divide(std::forward<U>(lhs), std::forward<V>(rhs));
}

// var 1
inline var abs(const std::complex<var>& z) { return complex_abs(z); }
// fvar 1
template <typename T>
inline fvar<T> abs(const std::complex<fvar<T>>& z) {
  return complex_abs(z);
}

// var 1
inline var arg(const std::complex<var>& z) { return complex_arg(z); }
// fvar 1
template <typename T>
inline fvar<T> arg(const std::complex<fvar<T>>& z) {
  return complex_arg(z);
}

// var 1
inline var norm(const std::complex<var>& z) { return complex_norm(z); }
// fvar 1
template <typename T>
inline fvar<T> norm(const std::complex<fvar<T>>& z) {
  return complex_norm(z);
}

// var 1
inline std::complex<var> conj(const std::complex<var>& z) {
  return complex_conj(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> conj(const std::complex<fvar<T>>& z) {
  return complex_conj(z);
}

// var 1
inline std::complex<var> proj(const std::complex<var>& z) {
  return complex_proj(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> proj(const std::complex<fvar<T>>& z) {
  return complex_proj(z);
}

// var 1
inline std::complex<var> exp(const std::complex<var>& z) {
  return complex_exp(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> exp(const std::complex<fvar<T>>& z) {
  return complex_exp(z);
}

// var 1
inline std::complex<var> log(const std::complex<var>& z) {
  return complex_log(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> log(const std::complex<fvar<T>>& z) {
  return complex_log(z);
}

// var 1
inline std::complex<var> log10(const std::complex<var>& z) {
  return complex_log10(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> log10(const std::complex<fvar<T>>& z) {
  return complex_log10(z);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @param[in] lhs first argument
 * @param[in] rhs second argument
 * @return first argument to the power of the second
 */
// var 1
inline std::complex<var> pow(int lhs, const std::complex<var>& rhs) {
  return complex_pow(lhs, rhs);
}
// var 2
inline std::complex<var> pow(double lhs, const std::complex<var>& rhs) {
  return complex_pow(lhs, rhs);
}
// var 3
inline std::complex<var> pow(const var& lhs, const std::complex<double>& rhs) {
  return complex_pow(lhs, rhs);
}
// var 4
inline std::complex<var> pow(const var& lhs, const std::complex<var>& rhs) {
  return complex_pow(lhs, rhs);
}
// var 5
inline std::complex<var> pow(const std::complex<double>& lhs, const var& rhs) {
  return complex_pow(lhs, rhs);
}
// var 6
inline std::complex<var> pow(const std::complex<double>& lhs,
                             const std::complex<var>& rhs) {
  return complex_pow(lhs, rhs);
}
// var 7
inline std::complex<var> pow(const std::complex<var>& lhs, int rhs) {
  return complex_pow(lhs, rhs);
}
// var 8
inline std::complex<var> pow(const std::complex<var>& lhs, double rhs) {
  return complex_pow(lhs, rhs);
}
// var 9
inline std::complex<var> pow(const std::complex<var>& lhs, const var& rhs) {
  return complex_pow(lhs, rhs);
}
// var 10
inline std::complex<var> pow(const std::complex<var>& lhs,
                             const std::complex<double>& rhs) {
  return complex_pow(lhs, rhs);
}
// var 11
inline std::complex<var> pow(const std::complex<var>& lhs,
                             const std::complex<var>& rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> pow(int lhs, const std::complex<fvar<T>>& rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 2
template <typename T>
inline std::complex<fvar<T>> pow(double lhs, const std::complex<fvar<T>>& rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 3
template <typename T>
inline std::complex<fvar<T>> pow(const fvar<T>& lhs,
                                 const std::complex<double>& rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 4
template <typename T>
inline std::complex<fvar<T>> pow(const fvar<T>& lhs,
                                 const std::complex<fvar<T>>& rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 5
template <typename T>
inline std::complex<fvar<T>> pow(const std::complex<double>& lhs,
                                 const fvar<T>& rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 6
template <typename T>
inline std::complex<fvar<T>> pow(const std::complex<double>& lhs,
                                 const std::complex<fvar<T>>& rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 7
template <typename T>
inline std::complex<fvar<T>> pow(const std::complex<fvar<T>>& lhs, int rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 8
template <typename T>
inline std::complex<fvar<T>> pow(const std::complex<fvar<T>>& lhs, double rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 9
template <typename T>
inline std::complex<fvar<T>> pow(const std::complex<fvar<T>>& lhs,
                                 const fvar<T>& rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 10
template <typename T>
inline std::complex<fvar<T>> pow(const std::complex<fvar<T>>& lhs,
                                 const std::complex<double>& rhs) {
  return complex_pow(lhs, rhs);
}
// fvar 11
template <typename T>
inline std::complex<fvar<T>> pow(const std::complex<fvar<T>>& lhs,
                                 const std::complex<fvar<T>>& rhs) {
  return complex_pow(lhs, rhs);
}

// var 1
inline std::complex<var> sqrt(const std::complex<var>& z) {
  return complex_sqrt(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> sqrt(const std::complex<fvar<T>>& z) {
  return complex_sqrt(z);
}

// var 1
inline std::complex<var> sin(const std::complex<var>& z) {
  return complex_sin(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> sin(const std::complex<fvar<T>>& z) {
  return complex_sin(z);
}

// var 1
inline std::complex<var> cos(const std::complex<var>& z) {
  return complex_cos(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> cos(const std::complex<fvar<T>>& z) {
  return complex_cos(z);
}

// var 1
inline std::complex<var> tan(const std::complex<var>& z) {
  return complex_tan(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> tan(const std::complex<fvar<T>>& z) {
  return complex_tan(z);
}

// var 1
inline std::complex<var> asin(const std::complex<var>& z) {
  return complex_asin(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> asin(const std::complex<fvar<T>>& z) {
  return complex_asin(z);
}

// var 1
inline std::complex<var> acos(const std::complex<var>& z) {
  return complex_acos(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> acos(const std::complex<fvar<T>>& z) {
  return complex_acos(z);
}

// var 1
inline std::complex<var> atan(const std::complex<var>& z) {
  return complex_atan(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> atan(const std::complex<fvar<T>>& z) {
  return complex_atan(z);
}

// var 1
inline std::complex<var> sinh(const std::complex<var>& z) {
  return complex_sinh(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> sinh(const std::complex<fvar<T>>& z) {
  return complex_sinh(z);
}

// var 1
inline std::complex<var> cosh(const std::complex<var>& z) {
  return complex_cosh(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> cosh(const std::complex<fvar<T>>& z) {
  return complex_cosh(z);
}

// var 1
inline std::complex<var> tanh(const std::complex<var>& z) {
  return complex_tanh(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> tanh(const std::complex<fvar<T>>& z) {
  return complex_tanh(z);
}

// var 1
inline std::complex<var> asinh(const std::complex<var>& z) {
  return complex_asinh(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> asinh(const std::complex<fvar<T>>& z) {
  return complex_asinh(z);
}

// var 1
inline std::complex<var> acosh(const std::complex<var>& z) {
  return complex_acosh(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> acosh(const std::complex<fvar<T>>& z) {
  return complex_acosh(z);
}

// var 1
inline std::complex<var> atanh(const std::complex<var>& z) {
  return complex_atanh(z);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> atanh(const std::complex<fvar<T>>& z) {
  return complex_atanh(z);
}

// var 1
inline std::complex<var> polar(const var& r, const var& theta) {
  return complex_polar(r, theta);
}
// var 2
inline std::complex<var> polar(const var& r, double theta) {
  return complex_polar(r, theta);
}
// var 3
inline std::complex<var> polar(double r, const var& theta) {
  return complex_polar(r, theta);
}
// var 4
inline std::complex<var> polar(const var& r, int theta) {
  return complex_polar(r, theta);
}
// var 5
inline std::complex<var> polar(int r, const var& theta) {
  return complex_polar(r, theta);
}
// fvar 1
template <typename T>
inline std::complex<fvar<T>> polar(const fvar<T>& r, const fvar<T>& theta) {
  return complex_polar(r, theta);
}
// fvar 2
template <typename T>
inline std::complex<fvar<T>> polar(const fvar<T>& r, double theta) {
  return complex_polar(r, theta);
}
// fvar 3
template <typename T>
inline std::complex<fvar<T>> polar(double r, const fvar<T>& theta) {
  return complex_polar(r, theta);
}
// fvar 4
template <typename T>
inline std::complex<fvar<T>> polar(const fvar<T>& r, int theta) {
  return complex_polar(r, theta);
}
// fvar 5
template <typename T>
inline std::complex<fvar<T>> polar(int r, const fvar<T>& theta) {
  return complex_polar(r, theta);
}

}  // namespace math
}  // namespace stan

#endif
