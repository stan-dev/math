#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/err/check_finite.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments:
 * \f$_pF_q(a_1,...,a_p;b_1,...,b_q;z)\f$
 *
 * See 'grad_pFq.hpp' for the derivatives wrt each parameter
 *
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @return Generalised hypergeometric function
 */
double hypergeometric_pFq(const Eigen::VectorXd& a, const Eigen::VectorXd& b,
                          double z) {
  check_finite("hypergeometric_pFq", "a", a);
  check_finite("hypergeometric_pFq", "b", b);
  check_finite("hypergeometric_pFq", "z", z);

  check_not_nan("hypergeometric_pFq", "a", a);
  check_not_nan("hypergeometric_pFq", "b", b);
  check_not_nan("hypergeometric_pFq", "z", z);

  bool condition_1 = (a.size() > (b.size() + 1)) && (z != 0);
  bool condition_2 = (a.size() == (b.size() + 1)) && (std::fabs(z) > 1);

  if (condition_1 || condition_2) {
    std::stringstream msg;
    msg << "hypergeometric function pFq does not meet convergence "
        << "conditions with given arguments. "
        << "a: " << a << ", b: " << b << ", "
        << ", z: " << z;
    throw std::domain_error(msg.str());
  }

  return boost::math::hypergeometric_pFq(
      std::vector<double>(a.data(), a.data() + a.size()),
      std::vector<double>(b.data(), b.data() + b.size()), z);
}
}  // namespace math
}  // namespace stan
#endif
