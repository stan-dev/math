#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HELPER_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HELPER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/to_array_1d.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace stan {
namespace math {
namespace internal {
/**
 * Implementation for calculating the generalized hypergeometric function
 * \f$_pF_q(a_1,...,a_p;b_1,...,b_q;z)\f$.
 *
 * This is declared separatel to avoid circular dependencies between the
 * various hypergeometric functions.
 *
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @return Generalized hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_vector_st<std::is_arithmetic, Ta, Tb>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline double hypergeometric_pFq_helper(const Ta& a, const Tb& b, const Tz& z) {
  return boost::math::hypergeometric_pFq(to_array_1d(a), to_array_1d(b), z);
}
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
