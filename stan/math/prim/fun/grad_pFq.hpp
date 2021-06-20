#ifndef STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP
#define STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_rising_factorial.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/max.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {

/**
 * Returns the gradient of generalised hypergeometric function wrt to the
 * input arguments:
 * \f$ _pF_q(a_1,...,a_p;b_1,...,b_q;z) \f$
 * 
 * The derivatives wrt p and q are defined using a Kampé de Fériet function,
 *   (https://en.wikipedia.org/wiki/Kamp%C3%A9_de_F%C3%A9riet_function).
 *   This is implemented below as an infinite sum (on the log scale) until
 *   convergence.
 * 
 * \f$ \frac{\partial}{\partial a_1} = 
 *     \frac{z\prod_{j=2}^pa_j}{\prod_{j=1}^qb_j}
 *     F_{q+1\:0\;1}^{p\:1\:2}\left(\begin{array}&a_1+1,...,a_p+1;1;1,a_1\\
 *        2, b_1+1,...,b_1+1;;a_1+1\end{array};z,z\right) \f$
 * 
 * \f$ \frac{\partial}{\partial b_1}= 
 *     -\frac{z\prod_{j=1}^pa_j}{b_1\prod_{j=1}^qb_j}
 *     F_{q+1\:0\;1}^{p\:1\:2}\left(\begin{array}&a_1+1,...,a_p+1;1;1,b_1\\
 *        2, b_1+1,...,b_1+1;;b_1+1\end{array};z,z\right) \f$
 *
 * \f$ \frac{\partial}{\partial z}= \frac{\prod_{j=1}^pa_j}{\prod_{j=1}^qb_j}
 *      {}_pF_q(a_1 + 1,...,a_p + 1; b_1 +1,...,b_q+1;z) \f$
 *
 * @tparam calc_p Boolean for whether to calculate derivatives wrt to 'p'
 * @tparam calc_q Boolean for whether to calculate derivatives wrt to 'q'
 * @tparam calc_z Boolean for whether to calculate derivatives wrt to 'z'
 * @tparam TupleT TType of tuple containing references
 * @tparam Tp Type of p
 * @tparam Tq Type of q
 * @tparam Tz Type of z
 * @param[in] grad_tuple Tuple of references to evaluate gradients into
 * @param[in] p Vector of 'p' arguments to function
 * @param[in] q Vector of 'q' arguments to function
 * @param[in] z Scalar z argument
 * @param[in] precision Convergence criteria for infinite sum
 * @param[in] max_steps Maximum number of iterations for infinite sum
 * @return Generalised hypergeometric function
 */
template <bool calc_p, bool calc_q, bool calc_z,
          typename TupleT, typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_t<Tp, Tq>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
void grad_pFq_impl(TupleT&& grad_tuple, const Tp& p, const Tq& q, const Tz& z,
              double precision = 1e-14, int max_steps = 1e6) {
  using std::max;
  using scalar_t = return_type_t<Tp, Tq, Tz>;
  using Tp_plain = plain_type_t<Tp>;
  using Tq_plain = plain_type_t<Tq>;
  using T_vec = Eigen::Matrix<scalar_t, -1, 1>;

  Tp_plain pp1 = (p.array() + 1).matrix();
  Tq_plain qp1 = (q.array() + 1).matrix();
  Tp_plain log_p = log(p);
  Tq_plain log_q = log(q);
  scalar_type_t<Tz> log_z = log(z);
  scalar_type_t<Tp> sum_log_p = sum(log_p);
  scalar_type_t<Tq> sum_log_q = sum(log_q);

  // Only need the infinite sum for partials wrt p & q
  if (calc_p || calc_q) {
    double log_precision = log(precision);

    int m = 0;
    scalar_t outer_diff = 0;

    // Declare vectors to accumulate sums into
    // where NEGATIVE_INFTY is zero on the log scale
     T_vec dp_infsum =  T_vec::Constant(p.size(), NEGATIVE_INFTY);
     T_vec dq_infsum =  T_vec::Constant(q.size(), NEGATIVE_INFTY);

    while ((outer_diff > log_precision) && (m < max_steps)) {
      // Vectors to accumulate outer sum into
      T_vec dp_iter_m =  T_vec::Constant(p.size(), NEGATIVE_INFTY);
      T_vec dp_mn =  T_vec::Constant(p.size(), NEGATIVE_INFTY);
      T_vec dq_iter_m = T_vec::Constant(q.size(), NEGATIVE_INFTY);
      T_vec dq_mn = T_vec::Constant(q.size(), NEGATIVE_INFTY);

      double log_phammer_1m = log_rising_factorial(1, m);
      double lgamma_mp1 = lgamma(m + 1);

      int n = 0;
      scalar_t inner_diff = 0;

      while ((inner_diff > log_precision) & (n < max_steps)) {
        // Numerator term
        scalar_t term1_mn = (m + n) * log_z
                            + sum(log_rising_factorial(pp1, m + n))
                            + log_phammer_1m + log_rising_factorial(1, n);
        // Denominator term
        scalar_t term2_mn = lgamma_mp1 + lgamma(n + 1)
                            + sum(log_rising_factorial(qp1, m + n))
                            + log_rising_factorial(2, m + n);
        if (calc_p) {
          // Division (on log scale) for the p & q partials
          dp_mn = (term1_mn + log_rising_factorial(p, n).array())
                    - (term2_mn + log_rising_factorial(pp1, n).array());

          // Perform a row-wise log_sum_exp to accumulate current iteration
          dp_iter_m = dp_iter_m.binaryExpr(dp_mn, [&](auto& a, auto& b) {
                                            return log_sum_exp(a, b);
                                          });
        }
        if (calc_q) {
          dq_mn = (term1_mn + log_rising_factorial(q, n).array())
                    - (term2_mn + log_rising_factorial(qp1, n).array());
          dq_iter_m = dq_iter_m.binaryExpr(dq_mn, [&](auto& a, auto& b) {
                                            return log_sum_exp(a, b);
                                          });
        }

        // Series convergence assessed by whether the sum of all terms is
        //   smaller than the specified criteria (precision)
        inner_diff = max(log_sum_exp(dp_mn), log_sum_exp(dq_mn));
        n += 1;
      }
      if (calc_p) {
        // Accumulate sums once the inner loop for the current iteration has
        //   converged
        dp_infsum = dp_infsum.binaryExpr(dp_iter_m, [&](auto& a, auto& b) {
                                          return log_sum_exp(a, b);
                                        });
      }
      if (calc_q) {
        dq_infsum = dq_infsum.binaryExpr(dq_iter_m, [&](auto& a, auto& b) {
                                        return log_sum_exp(a, b);
                                      });
      }

      // Assess convergence of outer loop
      outer_diff = max(log_sum_exp(dp_iter_m), log_sum_exp(dq_iter_m));
      m += 1;
    }
    if (m == max_steps) {
      throw_domain_error("grad_pFq", "k (internal counter)", max_steps,
                         "exceeded ",
                         " iterations, hypergeometric function gradient "
                         "did not converge.");
    }

    if (calc_p) {
      // Workaround to construct vector where each element is the product of
      //   all other elements
      Eigen::VectorXi ind_vector = Eigen::VectorXi::LinSpaced(p.size(), 0,
                                                              p.size() - 1);

      Tp_plain prod_excl_curr = ind_vector.unaryExpr([&log_p,
                                                      &sum_log_p] (int i) {
                                                  return sum_log_p - log_p[i];
                                                });
      T_vec pre_mult_p = (log_z + prod_excl_curr.array() - sum_log_q).matrix();

      // Evaluate gradients into provided containers
      std::get<0>(grad_tuple) = exp(pre_mult_p + dp_infsum);
    }

    if (calc_q) {
      T_vec pre_mult_q = (log_z + sum_log_p) - (log_q.array() + sum_log_q);
      std::get<1>(grad_tuple) = -exp(pre_mult_q + dq_infsum);
    }
  }

  if (calc_z) {
    std::get<2>(grad_tuple) = exp(sum_log_p - sum_log_q)
                                * hypergeometric_pFq(pp1, qp1, z);
  }
}
}  // namespace internal

/**
 * Wrapper function for calculating gradients wrt all parameters
 * 
 * @tparam Tp Type of p
 * @tparam Tq Type of q
 * @tparam Tz Type of z
 * @param[in] grad_p Vector to evaluate gradients wrt p into
 * @param[in] grad_q Vector to evaluate gradients wrt q into
 * @param[in] grad_z Scalar to evaluate gradient wrt z into
 * @param[in] p Vector of 'p' arguments to function
 * @param[in] q Vector of 'q' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Tp, typename Tq, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Tp, Tq>* = nullptr>
void grad_pFq(plain_type_t<Tp>& grad_p, plain_type_t<Tq>& grad_q,
              plain_type_t<Tz>& grad_z, const Tp& p, const Tq& q, const Tz& z,
              double precision = 1e-14, int max_steps = 1e6) {
  internal::grad_pFq_impl<true, true, true>(
    std::forward_as_tuple(grad_p, grad_q, grad_z),
    p, q, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt only p
 * 
 * @tparam Tp Type of p
 * @tparam Tq Type of q
 * @tparam Tz Type of z
 * @param[in] grad_p Vector to evaluate gradients wrt p into
 * @param[in] p Vector of 'p' arguments to function
 * @param[in] q Vector of 'q' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Tp, typename Tq, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Tp, Tq>* = nullptr>
void grad_pFq_p(plain_type_t<Tp>& grad_p, const Tp& p, const Tq& q,
                const Tz& z, double precision = 1e-14,
                int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Tp>;
  internal::grad_pFq_impl<true, false, false>(std::forward_as_tuple(grad_p,
                                                     plain_type_t<Tp>{},
                                                     scalar_type_t<Tp>{}),
                               p, q, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt only q
 * 
 * @tparam Tp Type of p
 * @tparam Tq Type of q
 * @tparam Tz Type of z
 * @param[in] grad_q Vector to evaluate gradients wrt q into
 * @param[in] p Vector of 'p' arguments to function
 * @param[in] q Vector of 'q' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Tp, typename Tq, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Tp, Tq>* = nullptr>
void grad_pFq_q(plain_type_t<Tq>& grad_q, const Tp& p, const Tq& q,
                const Tz& z, double precision = 1e-14,
                int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Tq>;
  internal::grad_pFq_impl<false, true, false>(
    std::forward_as_tuple(plain_type_t<Tq>{}, grad_q, scalar_type_t<Tq>{}),
    p, q, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt only z
 * 
 * @tparam Tp Type of p
 * @tparam Tq Type of q
 * @tparam Tz Type of z
 * @param[in] grad_z Scalar to evaluate gradient wrt z into
 * @param[in] p Vector of 'p' arguments to function
 * @param[in] q Vector of 'q' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Tp, typename Tq, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Tp, Tq>* = nullptr>
void grad_pFq_z(plain_type_t<Tz>& grad_z, const Tp& p, const Tq& q,
                const Tz& z, double precision = 1e-14, int max_steps = 1e6) {
  internal::grad_pFq_impl<false, false, true>(
    std::forward_as_tuple(Eigen::Matrix<Tz, -1, 1>{},
                          Eigen::Matrix<Tz, -1, 1>{}, grad_z),
    p, q, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt p and q
 * 
 * @tparam Tp Type of p
 * @tparam Tq Type of q
 * @tparam Tz Type of z
 * @param[in] grad_p Vector to evaluate gradients wrt p into
 * @param[in] grad_q Vector to evaluate gradients wrt q into
 * @param[in] p Vector of 'p' arguments to function
 * @param[in] q Vector of 'q' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Tp, typename Tq, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Tp, Tq>* = nullptr>
void grad_pFq_pq(plain_type_t<Tp>& grad_p, plain_type_t<Tq>& grad_q,
                 const Tp& p, const Tq& q, const Tz& z,
                 double precision = 1e-14, int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Tp>;
  internal::grad_pFq_impl<true, true, false>(
    std::forward_as_tuple(grad_p, grad_q,  scalar_type_t<Tp>{}),
    p, q, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt p and z
 * 
 * @tparam Tp Type of p
 * @tparam Tq Type of q
 * @tparam Tz Type of z
 * @param[in] grad_p Vector to evaluate gradients wrt p into
 * @param[in] grad_z Scalar to evaluate gradient wrt z into
 * @param[in] p Vector of 'p' arguments to function
 * @param[in] q Vector of 'q' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Tp, typename Tq, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Tp, Tq>* = nullptr>
void grad_pFq_pz(plain_type_t<Tp>& grad_p, plain_type_t<Tz>& grad_z,
                 const Tp& p, const Tq& q, const Tz& z,
                 double precision = 1e-14, int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Tp>;
  internal::grad_pFq_impl<true, false, true>(
    std::forward_as_tuple(grad_p, plain_type_t<Tp>{}, grad_z),
    p, q, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt q and z
 * 
 * @tparam Tp Type of p
 * @tparam Tq Type of q
 * @tparam Tz Type of z
 * @param[in] grad_q Vector to evaluate gradients wrt q into
 * @param[in] grad_z Scalar to evaluate gradient wrt z into
 * @param[in] p Vector of 'p' arguments to function
 * @param[in] q Vector of 'q' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Tp, typename Tq, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Tp, Tq>* = nullptr>
void grad_pFq_qz(plain_type_t<Tq>& grad_q, plain_type_t<Tz>& grad_z,
                 const Tp& p, const Tq& q, const Tz& z,
                 double precision = 1e-14, int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Tq>;
  internal::grad_pFq_impl<false, true, true>(
    std::forward_as_tuple(plain_type_t<Tq>{}, grad_q, grad_z),
    p, q, z, precision, max_steps);
}

}  // namespace math
}  // namespace stan
#endif
