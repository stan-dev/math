#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/err/check_finite.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/to_row_vector.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalized hypergeometric function applied to the
 * input arguments:
 * \f$_pF_q(a_1,...,a_p;b_1,...,b_q;z)\f$
 *
 * See 'grad_pFq.hpp' for the derivatives wrt each parameter
 *
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @return Generalized hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_vector_st<std::is_arithmetic, Ta, Tb>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_pFq(Ta&& a, Tb&& b, Tz&& z) {
  decltype(auto) a_ref = to_ref(std::forward<Ta>(a));
  decltype(auto) b_ref = to_ref(std::forward<Tb>(b));
  check_finite("hypergeometric_pFq", "a", a_ref);
  check_finite("hypergeometric_pFq", "b", b_ref);
  check_finite("hypergeometric_pFq", "z", z);
  check_not_nan("hypergeometric_pFq", "a", a_ref);
  check_not_nan("hypergeometric_pFq", "b", b_ref);
  check_not_nan("hypergeometric_pFq", "z", z);

  const bool condition_1 = (a_ref.size() > (b_ref.size() + 1)) && (z != 0);
  const bool condition_2
      = (a_ref.size() == (b_ref.size() + 1)) && (std::fabs(z) > 1);

  if (condition_1 || condition_2) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << "hypergeometric function pFq does not meet convergence "
             "conditions with given arguments. "
             "a: "
          << to_row_vector(a_ref) << ", "
          << "b: " << to_row_vector(b_ref) << ", "
          << "z: " << z;
      throw std::domain_error(msg.str());
    }();
  }
  // For plain vectors, we can use Eigen's Map to avoid unnecessary copies
  using a_ref_t = decltype(a_ref);
  using b_ref_t = decltype(b_ref);
  constexpr bool is_a_plain_vec = std::is_same_v<std::decay_t<a_ref_t>, plain_type_t<a_ref_t>>;
  constexpr bool is_b_plain_vec = std::is_same_v<std::decay_t<b_ref_t>, plain_type_t<b_ref_t>>;
  if constexpr (is_a_plain_vec && is_b_plain_vec) {
    // We use type erasure not do a hard copy here
    using map_t = Eigen::Map<Eigen::VectorXd>;
    auto map_a = map_t(const_cast<double*>(a_ref.data()), a_ref.size());
    auto map_b = map_t(const_cast<double*>(b_ref.data()), b_ref.size());
    return boost::math::hypergeometric_pFq(map_a, map_b, z);
  } else {
    // We need pointers to `a` and `b`'s data here so we hard evaluate.
    decltype(auto) a_eval = eval(a_ref);
    decltype(auto) b_eval = eval(b_ref);
    return boost::math::hypergeometric_pFq(
        std::vector<double>(a_eval.data(), a_eval.data() + a_eval.size()),
        std::vector<double>(b_eval.data(), b_eval.data() + b_eval.size()), z);
  }
}
}  // namespace math
}  // namespace stan
#endif
