#ifndef STAN_MATH_REV_ARR_FUNCTOR_FORWARD_PDE_HPP
#define STAN_MATH_REV_ARR_FUNCTOR_FORWARD_PDE_HPP

#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

#include <algorithm>
#include <vector>

namespace stan {
namespace math {

#ifdef STAN_EXTERN_PDE

/**
 * Return the solutions for the quantities of interest(QoI)
 * of the specified PDE problem.
 *
 * This function is templated to allow various PDE library
 * interfaces and corresponding input deck.
 *
 * @tparam F_pde type of PDE system interface.
 * @tparam T type of parameter.
 * @param[in] pde functor for the partial differential equation.
 * @param[in] theta parameter vector for the PDE.
 * @param[in] x_r continuous data vector for the PDE.
 * @param[in] x_i integer data vector for the PDE.
 * @param[out] msgs the print stream for warning messages.
 * @return a vector of quantities of interest.
 */
  template<typename F_pde, typename T>
  inline std::vector<T> forward_pde(const F_pde& pde,
                                         const std::vector<T>& theta,
                                         const std::vector<double>& x_r,
                                         const std::vector<int>& x_i,
                                         std::ostream* msgs = nullptr) {
    stan::math::check_not_nan("forward_pde", "theta", theta);

    const int need_sens = 1;
    std::vector<double> theta_d = stan::math::value_of(theta);
    std::vector<std::vector<double> > raw =
      pde(theta_d, need_sens, x_r, x_i, msgs);
    std::vector<T> res(raw.size());
    std::transform(raw.begin(), raw.end(),
                   res.begin(), [&theta](std::vector<double>& qoi_grad)
                   -> T {
                     double qoi = qoi_grad[0];
                     std::vector<double> g(qoi_grad.begin() + 1,
                                           qoi_grad.end());
                     return stan::math::precomputed_gradients(qoi, theta, g);
                   });
    return res;
  }

#endif
}
}

#endif
