#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/cumulative_sum.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

template <class RNG>
inline int categorical_rng(
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& theta, RNG& rng) {
  using boost::uniform_01;
  using boost::variate_generator;
  static constexpr const char* function = "categorical_rng";
  check_simplex(function, "Probabilities parameter", theta);

  variate_generator<RNG&, uniform_01<> > uniform01_rng(rng, uniform_01<>());

  Eigen::VectorXd index(theta.rows());
  index.setZero();

  index = cumulative_sum(theta);

  // index may sum to just below 1, so stop at the last nonzero entry
  int last = theta.size() - 1;
  while (theta(last) == 0) {
    last--;
  }
  double c = uniform01_rng();
  int b = 0;
  while (b < last && c >= index(b)) {
    b++;
  }
  return b + 1;
}
}  // namespace math
}  // namespace stan
#endif
