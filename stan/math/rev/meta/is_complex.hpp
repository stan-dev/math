#ifndef STAN_MATH_PRIM_META_IS_COMPLEX_HPP
#define STAN_MATH_PRIM_META_IS_COMPLEX_HPP

#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/is_complex.hpp>

#include <complex>
#include <type_traits>

namespace stan {
namespace internal {

template <typename... Ts>
struct is_complex_impl<var_value<std::complex<Ts...>>> : std::true_type {};

}  // namespace internal

}

#endif
