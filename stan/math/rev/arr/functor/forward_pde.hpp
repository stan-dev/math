#ifndef STAN_MATH_REV_ARR_FUNCTOR_FORWARD_PDE_HPP
#define STAN_MATH_REV_ARR_FUNCTOR_FORWARD_PDE_HPP

#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/rev/core/precomputed_gradients.hpp>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

#include <algorithm>
#include <vector>

namespace stan {
namespace math {

// #ifdef STAN_EXTERN_PDE

/**
 * Return the solutions for the quantities of interest(QoI)
 * of the specified PDE problem.
 *
 * This function is templated to allow various PDE library
 * interfaces and corresponding input deck.
 *
 * @tparam F_pde_qoi type of PDE system interface. The functor
 * signature should follow
 * operator()(const vector<double>&, // theta
 *            const int,             // calculate sensitivity?
 *            const vector<double>&, // x_r
 *            const vector<int>&,    // x_i
 *            std::ostream*) -> std::vector<std::vector<double> >
 * It returns a vector of vectors, with each member vector
 * in the form
 *
 * {QoI, d(QoI)/d(theta1), d(QoI)/d(theta2), ...},
 *
 * namely, the first component is the quantity of interest,
 * and the rest are its derivatives with respect to theta.
 * @tparam T type of parameter.
 * @param[in] pde_qoi functor for the partial differential
 * equation's QoI calculation.
 * @param[in] theta parameter vector for the PDE.
 * @param[in] x_r continuous data vector for the PDE.
 * @param[in] x_i integer data vector for the PDE.
 * @param[out] msgs the print stream for warning messages.
 * @return a vector of quantities of interest.
 */
template <typename F_pde_qoi, typename T>
inline std::vector<T> forward_pde(const F_pde_qoi& pde_qoi,
                                  const std::vector<T>& theta,
                                  const std::vector<double>& x_r,
                                  const std::vector<int>& x_i,
                                  std::ostream* msgs = nullptr) {
  stan::math::check_not_nan("forward_pde", "theta", theta);

  const int need_sens = 1;
  std::vector<double> theta_d = stan::math::value_of(theta);
  std::vector<std::vector<double> > raw
      = pde_qoi(theta_d, need_sens, x_r, x_i, msgs);
  std::vector<T> res(raw.size());
  std::transform(raw.begin(), raw.end(), res.begin(),
                 [&theta](std::vector<double>& qoi_grad) -> T {
                   double qoi = qoi_grad[0];
                   std::vector<double> g(qoi_grad.begin() + 1, qoi_grad.end());
                   return stan::math::precomputed_gradients(qoi, theta, g);
                 });
  return res;
}

// #endif
}  // namespace math
}  // namespace stan

#endif
