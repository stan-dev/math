#ifndef STAN_MATH_REV_FUN_ONE_CPT_HPP
#define STAN_MATH_REV_FUN_ONE_CPT_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/one_cpt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

  /**
   * Solve one-cpt model
   */
  template<typename Tt0, typename Tt1, typename T, typename T1>
  template <typename T, require_all_eigen_t<T, Tr>* = nullptr,
	    require_any_var_t<T, Tr, T1, T2>* = nullptr>
  inline auto one_cpt(const T& y,
		      const Tt& dt,
		      const T1& ka, const T2& k10,
		      const Tr& rate) {
    using stan::math::one_cpt;

    check_size_match("one_cpt", "y", y.rows(), "rate", rate.rows());
    const int n = y.rows();
    Eigen::Matrix<var, -1, 1> yt = Eigen::Matrix<var, -1, 1>::Zero(n);

    auto exp1 = exp(-k10_ * dt);
    auto exp2 = exp(-ka_ * dt);

    // contribution from cpt 1 bolus dose
    yt(1) += y(1) * exp1;

    // contribution from cpt 1 infusion dose
    yt(1) += rate[1] * (1 - exp1) / k10_;

    if (ka_ > 0.0) {
      // contribution from cpt 0 bolus dose
      yt(0) += y(0) * exp2;
      yt(1) += y(0) * ka_ / (ka_ - k10_) * (exp1 - exp2);

      // contribution from cpt 0 infusion dose
      yt(0) += rate[0] * (1 - exp2) / ka_;
      yt(1) += rate[0] * ka_ / (ka_ - k10_) * ((1 - exp1) / k10_ - (1 - exp2) / ka_);
    } else {
      // no absorption, GUT is accumulating dosages.
      yt(0) += y(0) + rate[0] * dt;
    }
    return yt;
  }

}  // namespace math
}  // namespace stan
#endif
