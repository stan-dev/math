#ifndef STAN_MATH_PRIM_PROB_ORDERED_PROBIT_RNG_HPP
#define STAN_MATH_PRIM_PROB_ORDERED_PROBIT_RNG_HPP

#include <stanh/prim/fun/Phi.hpp>
#include <stanh/prim/err/check_finite.hpp>
#include <stanh/prim/err/check_greater.hpp>
#include <stanh/prim/err/check_ordered.hpp>
#include <stanh/prim/prob/categorical_rng.hpp>

namespace stan {
namespace math {

template <class RNG>
inline int ordered_probit_rng(double eta, const Eigen::VectorXd& c, RNG& rng) {
  static const char* function = "ordered_probit";

  check_finite(function, "Location parameter", eta);
  check_greater(function, "Size of cut points parameter", c.size(), 0);
  check_ordered(function, "Cut points vector", c);
  check_finite(function, "Cut-points", c);

  Eigen::VectorXd cut(c.rows() + 1);
  cut(0) = 1 - Phi(eta - c(0));
  for (int j = 1; j < c.rows(); j++)
    cut(j) = Phi(eta - c(j - 1)) - Phi(eta - c(j));
  cut(c.rows()) = Phi(eta - c(c.rows() - 1));

  return categorical_rng(cut, rng);
}

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_PROB_ORDERED_PROBIT_RNG_HPP
#define STAN_MATH_PRIM_PROB_ORDERED_PROBIT_RNG_HPP

#include <stanh/prim/fun/Phi.hpp>
#include <stanh/prim/err/check_finite.hpp>
#include <stanh/prim/err/check_greater.hpp>
#include <stanh/prim/err/check_ordered.hpp>
#include <stanh/prim/prob/categorical_rng.hpp>

namespace stan {
namespace math {

template <class RNG>
inline int ordered_probit_rng(double eta, const Eigen::VectorXd& c, RNG& rng) {
  static const char* function = "ordered_probit";

  check_finite(function, "Location parameter", eta);
  check_greater(function, "Size of cut points parameter", c.size(), 0);
  check_ordered(function, "Cut points vector", c);
  check_finite(function, "Cut-points", c);

  Eigen::VectorXd cut(c.rows() + 1);
  cut(0) = 1 - Phi(eta - c(0));
  for (int j = 1; j < c.rows(); j++)
    cut(j) = Phi(eta - c(j - 1)) - Phi(eta - c(j));
  cut(c.rows()) = Phi(eta - c(c.rows() - 1));

  return categorical_rng(cut, rng);
}

}  // namespace math
}  // namespace stan
#endif
