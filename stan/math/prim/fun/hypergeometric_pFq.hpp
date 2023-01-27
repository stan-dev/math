#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/err/check_finite.hpp>
#include <stan/math/prim/fun/to_array_1d.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/to_row_vector.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments:
 * \f$_pF_q(a_1,...,a_p;b_1,...,b_q;z)\f$
 *
 * This function is not intended to be exposed to end users, only
 * used for p & q values that are stable with the grad_pFq
 * implementation.
 *
 * See 'grad_pFq.hpp' for the derivatives wrt each parameter
 *
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          typename ScalarT = return_type_t<Ta, Tb, Tz>,
          require_all_vector_st<std::is_arithmetic, Ta, Tb>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
return_type_t<Ta, Tb, Tz> hypergeometric_pFq(const Ta& a, const Tb& b,
                                             const Tz& z) {
  std::vector<ScalarT> a_ref = to_array_1d(a);
  std::vector<ScalarT> b_ref = to_array_1d(b);
  check_finite("hypergeometric_pFq", "a", a_ref);
  check_finite("hypergeometric_pFq", "b", b_ref);
  check_finite("hypergeometric_pFq", "z", z);

  check_not_nan("hypergeometric_pFq", "a", a_ref);
  check_not_nan("hypergeometric_pFq", "b", b_ref);
  check_not_nan("hypergeometric_pFq", "z", z);

  bool condition_1 = (a_ref.size() > (b_ref.size() + 1)) && (z != 0);
  bool condition_2 = (a_ref.size() == (b_ref.size() + 1)) && (std::fabs(z) > 1.0);

  if (condition_1 || condition_2) {
    std::stringstream msg;
    msg << "hypergeometric function pFq does not meet convergence "
        << "conditions with given arguments: \n"
        << "a: [" << to_row_vector(a) << "]\n"
        << "b: [" << to_row_vector(b) << "]\n"
        << "z: " << z;
    throw std::domain_error(msg.str());
  }

  return boost::math::hypergeometric_pFq(a_ref, b_ref, z);
}
}  // namespace math
}  // namespace stan
#endif
