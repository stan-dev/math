#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/cumulative_sum.hpp>
#include <random>

namespace stan {
namespace math {

template <class RNG>
inline int categorical_rng(
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& theta, RNG& rng) {
  static constexpr const char* function = "categorical_rng";
  check_simplex(function, "Probabilities parameter", theta);

  std::uniform_real_distribution<> uniform01_rng(0, 1);

  Eigen::VectorXd index(theta.rows());
  index.setZero();

  index = cumulative_sum(theta);

  double c = uniform01_rng(rng);
  int b = 0;
  while (c > index(b, 0)) {
    b++;
  }
  return b + 1;
}
}  // namespace math
}  // namespace stan
#endif
