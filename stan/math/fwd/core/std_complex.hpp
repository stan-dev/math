#ifndef STAN_MATH_FWD_CORE_STD_COMPLEX_HPP
#define STAN_MATH_FWD_CORE_STD_COMPLEX_HPP

#include <stan/math/prim/core/complex_base.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/prim/meta.hpp>
#include <complex>

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

   constexpr complex() = default;
   /**
    * Construct complex number from real and imaginary parts.
    *
    * @tparam U type of real part
    * @tparam V type of imaginary part
    * @param[in] re real part
    * @param[in] im imaginary part
    */
   template <typename U, typename V>
   constexpr complex(const U& re, const V& im) : base_t(re, im) {}

   /**
    * Construct a complex base with the specified real part and a zero
    * imaginary part.
    *
    * @tparam Scalar real type (must be assignable to `V`)
    * @param[in] re real part
    */
   template <typename Scalar, typename = stan::require_stan_scalar_t<T>>
   constexpr complex(Scalar&& re) : base_t(re) {} // NOLINT(runtime/explicit)

   template <typename Scalar, typename = stan::require_arithmetic_t<Scalar>>
   constexpr auto& operator=(const std::complex<Scalar>& x) {
     this->re_ = x.real();
     this->im_ = x.imag();
     return *this;
   }

};

}  // namespace std

#endif
