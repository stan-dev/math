#ifndef STAN_MATH_FWD_FUN_GRAD_INC_BETA_HPP
#define STAN_MATH_FWD_FUN_GRAD_INC_BETA_HPP

#include <stan/math/prim/fun/grad_2F1.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/beta.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/inc_beta.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/fwd/fun/log1m.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Gradient of the incomplete beta function beta(a, b, z) with
 * respect to the first two arguments.
 *
 * Uses the equivalence to a hypergeometric function. See
 * http://dlmf.nist.gov/8.17#ii
 *
 * @tparam T inner type of the fvar
 * @param[out] g1 d/da
 * @param[out] g2 d/db
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 */
template <typename T>
void grad_inc_beta(fvar<T>& g1, fvar<T>& g2, fvar<T> a, fvar<T> b, fvar<T> z) {
  fvar<T> c1 = log(z);
  fvar<T> c2 = log1m(z);
  fvar<T> c3 = beta(a, b) * inc_beta(a, b, z);

  fvar<T> C = exp(a * c1 + b * c2) / a;

  fvar<T> dF1 = 0;
  fvar<T> dF2 = 0;
  fvar<T> dF3 = 0;
  fvar<T> dFz = 0;

  if (value_of_rec(C)) {
    std::forward_as_tuple(dF1, dF2, dF3, dFz)
        = grad_2F1<true>(a + b, fvar<T>(1.0), a + 1, z);
  }

  g1 = (c1 - 1.0 / a) * c3 + C * (dF1 + dF3);
  g2 = c2 * c3 + C * dF1;
}

}  // namespace math
}  // namespace stan
#endif
