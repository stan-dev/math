#ifndef STAN_MATH_FWD_CORE_STD_COMPLEX_HPP
#define STAN_MATH_FWD_CORE_STD_COMPLEX_HPP

#include <stan/math/prim/core/complex_base.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <complex>
#include <type_traits>

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
   * Constructs complex number from real part or with default zero
   * value, setting imaginary part to zero.  This is also the nullary
   * constructor.
   *
   * @param[in] re the real part (default zero).
   */
  complex(const value_type& re = value_type(0))  // NOLINT(runtime/explicit)
      : stan::math::complex_base<stan::math::fvar<T>>(re) {}

  /**
   * Construct a complex number with the real and imaginary components
   * that are copies of those of the specified complex number.
   *
   * @tparam V value type of complex argument
   * @param[in] other another complex to use as source
   */
  template <typename V>
  complex(const complex<V>& other)
      : stan::math::complex_base<stan::math::fvar<T>>(other) {}

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
   * Construct a complex number with the specified real component and
   * zero imaginary component.
   *
   * @tparam V type of real component
   * @param[in] re real component
   */
  template <typename V,
            typename = std::enable_if_t<std::is_arithmetic<V>::value>>
  complex(V re)  // NOLINT(runtime/explicit)
      : stan::math::complex_base<stan::math::fvar<T>>(re) {}

  /**
   * Destroy this complex number.
   */
  ~complex() {}

  /**
   * Assign the specified value to the real part of this complex
   * number and set imaginary part to zero.
   *
   * @tparam V type of value
   * @param[in] x value to assign
   * @return this complex number after assigning other
   */
  template <typename V>
  complex_type& operator=(const V& x) {
    return base_t::operator=(x);
  }

  /**
   * Adds other to this and return this.
   *
   * @tparam V type of value
   * @param[in] other value to add
   * @return this complex number after adding argument
   */
  template <typename V>
  complex_type& operator+=(const V& other) {
    return base_t::operator+=(other);
  }

  /**
   * Subtracts other from this and return this.
   *
   * @tparam V type of value
   * @param[in] other value to subtract
   * @return this complex number after subtracting argument
   */
  template <typename V>
  complex_type& operator-=(const V& other) {
    return base_t::operator-=(other);
  }

  /**
   * Multiplies this by other and return this.
   *
   * @tparam V type of value
   * @param[in] other value to multiply
   * @return this complex number after multiplying by argument
   */
  template <typename V>
  complex_type& operator*=(const V& other) {
    return base_t::operator*=(other);
  }

  /**
   * Divides this by other and return this.
   *
   * @tparam V type of value
   * @param[in] other value to divide by
   * @return this complex number after dividing by argument
   */
  template <typename V>
  complex_type& operator/=(const V& other) {
    return base_t::operator/=(other);
  }
};

}  // namespace std

#endif
