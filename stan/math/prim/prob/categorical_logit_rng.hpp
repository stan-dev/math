#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_RNG_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/cumulative_sum.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a draw from a Categorical distribution given a
 * a vector of unnormalized log probabilities and a psuedo-random
 * number generator.
 *
 * This is a convenience wrapper around
 * <code>categorical_rng(softmax(beta), rng)</code>.
 *
 * @tparam RNG Type of pseudo-random number generator.
 * @param beta Vector of unnormalized log probabilities.
 * @param rng Pseudo-random number generator.
 * @return Categorical random variate
 */
template <typename T_beta, class RNG>
inline auto categorical_logit_rng(const T_beta& beta, RNG& rng) {
  using boost::uniform_01;
  using boost::variate_generator;

  vector_seq_view<T_beta> beta_vec(beta);
  size_t N = size_mvt(beta);

  VectorBuilder<true, int, value_type_t<T_beta>> output(N);
  variate_generator<RNG&, uniform_01<>> uniform01_rng(rng, uniform_01<>());

  for (size_t n = 0; n < N; ++n) {
    check_finite("categorical_logit_rng", "Probabilities parameter",
                 beta_vec[n]);

    const auto& theta = softmax(beta_vec[n]);
    Eigen::VectorXd index = cumulative_sum(theta);

    double c = uniform01_rng();
    int b = 0;

    while (c > index(b)) {
      b++;
    }
    output[n] = b + 1;
  }

  return output.data();
}
}  // namespace math
}  // namespace stan
#endif
