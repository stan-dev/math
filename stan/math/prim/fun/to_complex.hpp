#ifndef STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP
#define STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return a complex value from a real component and an imaginary component.
 * Default values for both components is 0.
 *
 * @tparam T type of real component
 * @tparam S type of imaginary component
 * @param[in] re real component (default = 0)
 * @param[in] im imaginary component (default = 0)
 * @return complex value with specified real and imaginary components
 */
template <typename T = double, typename S = double,
          require_all_not_container_t<T, S>* = nullptr>
constexpr inline std::complex<stan::real_return_t<T, S>> to_complex(
    const T& re = 0, const S& im = 0) {
  return std::complex<stan::real_return_t<T, S>>(re, im);
}

/**
 * Return a complex valued container
 * from a real component and an imaginary component.
 *
 * @tparam T type of real component
 * @tparam S type of imaginary component
 * @param[in] re real component
 * @param[in] im imaginary component
 * @return complex valued container with specified real and imaginary components
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr,
          require_all_st_stan_scalar<T1, T2>* = nullptr>
inline auto to_complex(const T1& re, const T2& im) {
  return apply_scalar_binary(re, im, [&](const auto& c, const auto& d) {
    return stan::math::to_complex(c, d);
  });
}

}  // namespace math
}  // namespace stan

#endif
