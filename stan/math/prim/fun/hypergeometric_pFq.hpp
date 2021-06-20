#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments
 *
 * @param[in] p p
 * @param[in] q q
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
double hypergeometric_pFq(const Eigen::VectorXd& p, const Eigen::VectorXd& q,
                          double z) {
  using std::fabs;
  check_finite("hypergeometric_pFq", "p", p);
  check_finite("hypergeometric_pFq", "q", q);
  check_finite("hypergeometric_pFq", "z", z);

  check_not_nan("hypergeometric_pFq", "p", p);
  check_not_nan("hypergeometric_pFq", "q", q);
  check_not_nan("hypergeometric_pFq", "z", z);

  bool condition_1 = (p.size() > (q.size() + 1));
  bool condition_2 = (p.size() == (q.size() + 1)) && (fabs(z) > 1);
  bool condition_3 = (p.size() > (q.size() + 1)) && (z != 0);

  if (condition_1 || condition_2 || condition_3) {
    std::stringstream msg;
    msg << "hypergeometric function pFq does not meet convergence "
        << "conditions with given arguments. "
        << "p: " << p << ", q: " << q << ", "
        << ", z: " << z;
    throw std::domain_error(msg.str());
  }

  return boost::math::hypergeometric_pFq(
                std::vector<double>(p.data(), p.data() + p.size()),
                std::vector<double>(q.data(), q.data() + q.size()),
                z);
}
}  // namespace math
}  // namespace stan
#endif
