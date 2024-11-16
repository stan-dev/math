#ifndef STAN_MATH_PRIM_FUN_ONE_CPT_HPP
#define STAN_MATH_PRIM_FUN_ONE_CPT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

  /**
   * Solve one-cpt model
   */
  template <typename T, typename Tr,
	    require_all_eigen_t<T, Tr>* = nullptr,
	    require_all_vt_arithmetic<T, Tr>* = nullptr>
  inline auto one_cpt(const T& y,
		      double dt, double ka, double k10,
		      const Tr& rate) {
    using stan::math::exp;

    check_size_match("one_cpt", "y", y.rows(), "rate", rate.rows());
    const int n = y.rows();
    Eigen::Matrix<double, -1, 1> yt = Eigen::Matrix<double, -1, 1>::Zero(n);

    double exp1 = exp(-k10 * dt);
    double exp2 = exp(-ka * dt);

    // contribution from cpt 1 bolus dose
    yt(1) += y(1) * exp1;

    // contribution from cpt 1 infusion dose
    yt(1) += rate(1) * (1 - exp1) / k10;

    if (ka > 0.0) {
      // contribution from cpt 0 bolus dose
      yt(0) += y(0) * exp2;
      yt(1) += y(0) * ka / (ka - k10) * (exp1 - exp2);

      // contribution from cpt 0 infusion dose
      yt(0) += rate(0) * (1 - exp2) / ka;
      yt(1) += rate(0) * ka / (ka - k10) * ((1 - exp1) / k10 - (1 - exp2) / ka);
    } else {
      // no absorption, GUT is accumulating dosages.
      yt(0) += y(0) + rate(0) * dt;
    }
    return yt;
  }

}  // namespace math
}  // namespace stan
#endif
