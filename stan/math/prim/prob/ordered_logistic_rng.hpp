#ifndef STAN_MATH_PRIM_PROB_ORDERED_LOGISTIC_RNG_HPP
#define STAN_MATH_PRIM_PROB_ORDERED_LOGISTIC_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/prob/categorical_rng.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return an ordered-logistic random variate with specified location(s)
 * and cutpoint vector(s) using the specified random number generator.
 *
 * @tparam T_eta type of location parameter
 * @tparam T_c type of cutpoint parameter
 * @tparam RNG type of random number generator
 * @param eta (Sequence of) location parameters
 * @param c (Sequence of) cutpoint vectors
 * @param rng random number generator
 * @return (Sequence of) Ordered-logistic random variate(s)
 */
template <typename T_eta, typename T_c, class RNG>
inline auto ordered_logistic_rng(const T_eta& eta, const T_c& c, RNG& rng) {
  static const char* function = "ordered_logistic";
  check_greater(function, "Size of eta parameter", stan::math::size(eta), 0);
  check_greater(function, "Size of cut points parameter", stan::math::size(c),
                0);
  scalar_seq_view<T_eta> eta_vec(eta);
  vector_seq_view<T_c> cut_vec(c);

  if (stan::math::size_mvt(c) > 1) {
    check_consistent_sizes(function, "Locations", eta, "Cut points vector", c);
  }

  for (int c_i = 0; c_i < size_mvt(c); ++c_i) {
    check_greater(function, "Size of cut points parameter", cut_vec[c_i].size(),
                  0);
    check_ordered(function, "Cut points parameter", cut_vec[c_i]);
    check_finite(function, "Cut points parameter",
                 cut_vec[c_i][cut_vec[c_i].size() - 1]);
    check_finite(function, "Cut points parameter", cut_vec[c_i][0]);
  }

  for (int e_i = 0; e_i < stan::math::size(eta); ++e_i) {
    check_finite(function, "Location parameter", eta_vec[e_i]);
  }

  int ret_size = std::max(stan::math::size(eta), size_mvt(c));

  using return_t = std::conditional_t<
      disjunction<is_vector<T_eta>, is_std_vector<T_c>>::value,
      std::vector<int>, int>;

  VectorBuilder<true, int, return_t> output(ret_size);
  for (int i = 0; i < ret_size; ++i) {
    Eigen::VectorXd cut(cut_vec[i].rows() + 1);
    cut[0] = 1 - inv_logit(eta_vec[i] - cut_vec[i][0]);
    for (int j = 1; j < cut_vec[i].rows(); j++) {
      cut[j] = inv_logit(eta_vec[i] - cut_vec[i][j - 1])
               - inv_logit(eta_vec[i] - cut_vec[i][j]);
    }
    cut[cut_vec[i].rows()]
        = inv_logit(eta_vec[i] - cut_vec[i][cut_vec[i].rows() - 1]);

    output[i] = categorical_rng(cut, rng);
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
