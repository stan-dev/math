#ifndef STAN_MATH_PRIM_PROB_ORDERED_LOGISTIC_RNG_HPP
#define STAN_MATH_PRIM_PROB_ORDERED_LOGISTIC_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/prob/categorical_rng.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

template <class RNG>
inline int ordered_logistic_rng(
    double eta, const Eigen::Matrix<double, Eigen::Dynamic, 1>& c, RNG& rng) {
  using boost::variate_generator;
  static const char* function = "ordered_logistic";
  check_finite(function, "Location parameter", eta);
  check_greater(function, "Size of cut points parameter", c.size(), 0);
  check_ordered(function, "Cut points parameter", c);
  check_finite(function, "Cut points parameter", c(c.size() - 1));
  check_finite(function, "Cut points parameter", c(0));

  Eigen::VectorXd cut(c.rows() + 1);
  cut(0) = 1 - inv_logit(eta - c(0));
  for (int j = 1; j < c.rows(); j++) {
    cut(j) = inv_logit(eta - c(j - 1)) - inv_logit(eta - c(j));
  }
  cut(c.rows()) = inv_logit(eta - c(c.rows() - 1));

  return categorical_rng(cut, rng);
}

template <typename T_eta, typename T_c, class RNG,
          require_t<disjunction<is_vector<T_eta>, 
                                conjunction<is_std_vector<T_c>,
                                            is_vector<value_type_t<T_c>>
                                            >>>* = nullptr>
inline std::vector<int> ordered_logistic_rng(
  const T_eta& eta, const T_c& c, RNG& rng) {
  static const char* function = "ordered_logistic";
  check_greater(function, "Size of eta parameter", stan::math::size(eta), 0);
  check_greater(function, "Size of cut points parameter", stan::math::size(c), 0);
  scalar_seq_view<T_eta> eta_vec(eta);
  vector_seq_view<T_c> cut_vec(c);

  int ret_size = std::max(stan::math::size(eta), size_mvt(c));

  std::vector<int> rng_return(ret_size);
  for (int i = 0; i < ret_size; ++i) {
    rng_return[i] = ordered_logistic_rng(eta_vec[i], cut_vec[i], rng);
  }

  return rng_return;
}

}  // namespace math
}  // namespace stan
#endif
