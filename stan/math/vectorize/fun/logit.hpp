
#ifndef STAN_MATH_VECTORIZED_FUN_LOGIT_HPP
#define STAN_MATH_VECTORIZED_FUN_LOGIT_HPP
#include <stan/math/prim/fun/logit.hpp>
namespace stan {
namespace math {

/**
 * Version of logit() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return the logit of each variable in the container.
 *
 * Note: The return must be evaluated otherwise the Ref object falls out
 * of scope
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto logit(const Container& x) {
  return make_holder(
      [](const auto& v_ref) {
        return apply_vector_unary<ref_type_t<Container>>::apply(
            v_ref,
            [](const auto& v) { return (v.array() / (1 - v.array())).log(); });
      },
      to_ref(x));
}


} // namespace math
} // namespace stan
#endif 

