#ifndef STAN_MATH_PRIM_SCAL_META_IS_DETECTED_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_DETECTED_HPP

#include <stan/math/prim/scal/meta/void_t.hpp>
#include <type_traits>

namespace stan {

template<typename, template <typename> class, typename = void>
struct is_detected : std::false_type {};

template<typename T, template <typename> class Op>
struct is_detected<T, Op, void_t<Op<T>>> : std::true_type {};


}  // namespace stan
#endif
