#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_SOLVE_PDE_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_SOLVE_PDE_HPP

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
 * interfaces and corresponding input deck. This is
 * data-only version of the function. Therefore there is no
 * sensitivity is requested.
 *
 * @tparam F_pde type of PDE system interface. The functor
 * signature should follow
 * operator()(const vector<double>&, // theta
 *            const int,             // calculate sensitivity?
 *            const vector<double>&, // x_r
 *            const vector<int>&,    // x_i
 *            std::ostream*) -> std::vector<std::vector<double> >
 * It returns a vector of vectors, with each member vector
 * in the form
 *
 * {QoI}
 *
 * namely, a single-element vector of the quantity of interest.
 *
 * @param[in] pde functor for the partial differential equation.
 * @param[in] theta parameter vector for the PDE.
 * @param[in] x_r continuous data vector for the PDE.
 * @param[in] x_i integer data vector for the PDE.
 * @param[out] msgs the print stream for warning messages.
 * @return a vector of quantities of interest.
 */
template <typename F_pde>
inline std::vector<double> solve_pde(const F_pde& pde,
                                     const std::vector<double>& theta,
                                     const std::vector<double>& x_r,
                                     const std::vector<int>& x_i,
                                     std::ostream* msgs = nullptr) {
  stan::math::check_not_nan("solve_pde", "theta", theta);
  const int need_sens = 0;
  std::vector<std::vector<double> > raw = pde(theta, need_sens, x_r, x_i, msgs);
  std::vector<double> res(raw.size());
  std::transform(raw.begin(), raw.end(), res.begin(),
                 [&theta](std::vector<double>& qoi_grad) -> double {
                   return qoi_grad[0];
                 });
  return res;
}

// #endif
}
}  // namespace stan

#endif
