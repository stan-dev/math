#ifndef STAN_MATH_PRIM_CPLX_COMPLEX_HPP
#define STAN_MATH_PRIM_CPLX_COMPLEX_HPP

#include <stan/math/prim/cplx/zeroing.hpp>
#include <complex>

namespace stan {
namespace math {

/** 
 *  Complex class that forwards the interface of std::complex, brining
 *  stan::math's namespace into ADL for functions called on data types
 *  that inherit from this class, as well as allowing template
 *  specializations of std::complex to indirectly work with different
 *  underlying data types. For example, std::complex<var> specializes
 *  the std::complex<T> template to inherit from
 *  stan::math::complex<T>, which in turn inherits from 
 *  std::complex<zeroing<T>>. This causes stan::math to be pulled into 
 *  ADL on operations involving std::complex<var>, to handle, e.g.
 *  division operations that the base std::complex doesn't have
 *  defined usably for var, and, since zeroing<var> zero-initializes
 *  itself, will also correctly work with the remaining algorithms
 *  in std::complex<T> that require the expression T() to be 0.
 */
template <class T>
struct complex : std::complex<zeroing<T>> {
  using std::complex<zeroing<T>>::complex;
  // NOLINTNEXTLINE(runtime/explicit)
  complex(const std::complex<zeroing<T>>& t)
      : std::complex<zeroing<T>>(t) {}  ///< allow promotion
  /**
   * Without this converting ctor, copy initialization of
   * std::complex<zeroing<(f)var>>> from an (f)var fails, because the 
   * default implementation in the std::complex template requires the
   * argument type to match the template type. This is the same reason
   * why std::complex<double> can't be multiplied by an int.
   *
   * tparam R type of real argument. Won't fire on non-arithmetic types.
   * tparam I type of imaginary argument.
   * param[in] real, the pure real component of the complex number.
   * param[in] imag, the pure imaginary component of the complex number
   */
  template <class R = double, class I = double,
            std::enable_if_t<is_arith_like<R>::value>* = nullptr>
  // NOLINTNEXTLINE(runtime/explicit)
  complex(R const& real =0, I const& imag =0) 
      : std::complex<zeroing<T>>(real, imag) {}
  // downcast
  operator std::complex<T>() const {
   return *static_cast<std::complex<T>*>(this);
  }
};

}  // namespace math
}  // namespace stan
#endif
