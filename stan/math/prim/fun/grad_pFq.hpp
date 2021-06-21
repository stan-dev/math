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
 * The derivatives wrt a and b are defined using a Kampé de Fériet function,
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
 * @tparam calc_a Boolean for whether to calculate derivatives wrt to 'a'
 * @tparam calc_b Boolean for whether to calculate derivatives wrt to 'b'
 * @tparam calc_z Boolean for whether to calculate derivatives wrt to 'z'
 * @tparam TupleT Type of tuple containing references
 * @tparam Ta Type of a
 * @tparam Tb Type of b
 * @tparam Tz Type of z
 * @param[in] grad_tuple Tuple of references to evaluate gradients into
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @param[in] precision Convergence criteria for infinite sum
 * @param[in] max_steps Maximum number of iterations for infinite sum
 * @return Generalised hypergeometric function
 */
template <bool calc_a, bool calc_b, bool calc_z,
          typename TupleT, typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
void grad_pFq_impl(TupleT&& grad_tuple, const Ta& a, const Tb& b, const Tz& z,
              double precision = 1e-14, int max_steps = 1e6) {
  using std::max;
  using scalar_t = return_type_t<Ta, Tb, Tz>;
  using Ta_plain = plain_type_t<Ta>;
  using Tb_plain = plain_type_t<Tb>;
  using T_vec = Eigen::Matrix<scalar_t, -1, 1>;

  Ta_plain ap1 = (a.array() + 1).matrix();
  Tb_plain bp1 = (b.array() + 1).matrix();
  Ta_plain log_a = log(a);
  Tb_plain log_b = log(b);
  scalar_type_t<Tz> log_z = log(z);
  scalar_type_t<Ta> sum_log_a = sum(log_a);
  scalar_type_t<Tb> sum_log_b = sum(log_b);

  // Only need the infinite sum for partials wrt p & b
  if (calc_a || calc_b) {
    double log_precision = log(precision);

    int m = 0;
    scalar_t outer_diff = 0;

    // Declare vectors to accumulate sums into
    // where NEGATIVE_INFTY is zero on the log scale
     T_vec da_infsum = T_vec::Constant(a.size(), NEGATIVE_INFTY);
     T_vec db_infsum = T_vec::Constant(b.size(), NEGATIVE_INFTY);

    while ((outer_diff > log_precision) && (m < max_steps)) {
      // Vectors to accumulate outer sum into
      T_vec da_iter_m = T_vec::Constant(a.size(), NEGATIVE_INFTY);
      T_vec da_mn = T_vec::Constant(a.size(), NEGATIVE_INFTY);
      T_vec db_iter_m = T_vec::Constant(b.size(), NEGATIVE_INFTY);
      T_vec db_mn = T_vec::Constant(b.size(), NEGATIVE_INFTY);

      double log_phammer_1m = log_rising_factorial(1, m);
      double lgamma_mp1 = lgamma(m + 1);

      int n = 0;
      scalar_t inner_diff = 0;

      while ((inner_diff > log_precision) & (n < max_steps)) {
        // Numerator term
        scalar_t term1_mn = (m + n) * log_z
                            + sum(log_rising_factorial(ap1, m + n))
                            + log_phammer_1m + log_rising_factorial(1, n);
        // Denominator term
        scalar_t term2_mn = lgamma_mp1 + lgamma(n + 1)
                            + sum(log_rising_factorial(bp1, m + n))
                            + log_rising_factorial(2, m + n);
        if (calc_a) {
          // Division (on log scale) for the a & b partials
          da_mn = (term1_mn + log_rising_factorial(a, n).array())
                    - (term2_mn + log_rising_factorial(ap1, n).array());

          // Perform a row-wise log_sum_exp to accumulate current iteration
          da_iter_m = da_iter_m.binaryExpr(da_mn, [&](auto& a, auto& b) {
                                            return log_sum_exp(a, b);
                                          });
        }
        if (calc_b) {
          db_mn = (term1_mn + log_rising_factorial(b, n).array())
                    - (term2_mn + log_rising_factorial(bp1, n).array());
          db_iter_m = db_iter_m.binaryExpr(db_mn, [&](auto& a, auto& b) {
                                            return log_sum_exp(a, b);
                                          });
        }

        // Series convergence assessed by whether the sum of all terms is
        //   smaller than the specified criteria (precision)
        inner_diff = max(log_sum_exp(da_mn), log_sum_exp(db_mn));
        n += 1;
      }
      if (calc_a) {
        // Accumulate sums once the inner loop for the current iteration has
        //   converged
        da_infsum = da_infsum.binaryExpr(da_iter_m, [&](auto& a, auto& b) {
                                          return log_sum_exp(a, b);
                                        });
      }
      if (calc_b) {
        db_infsum = db_infsum.binaryExpr(db_iter_m, [&](auto& a, auto& b) {
                                        return log_sum_exp(a, b);
                                      });
      }

      // Assess convergence of outer loop
      outer_diff = max(log_sum_exp(da_iter_m), log_sum_exp(db_iter_m));
      m += 1;
    }
    if (m == max_steps) {
      throw_domain_error("grad_pFq", "k (internal counter)", max_steps,
                         "exceeded ",
                         " iterations, hypergeometric function gradient "
                         "did not converge.");
    }

    if (calc_a) {
      // Workaround to construct vector where each element is the product of
      //   all other elements
      Eigen::VectorXi ind_vector = Eigen::VectorXi::LinSpaced(a.size(), 0,
                                                              a.size() - 1);

      Ta_plain prod_excl_curr = ind_vector.unaryExpr([&log_a,
                                                      &sum_log_a] (int i) {
                                                  return sum_log_a - log_a[i];
                                                });
      T_vec pre_mult_a = (log_z + prod_excl_curr.array() - sum_log_b).matrix();

      // Evaluate gradients into provided containers
      std::get<0>(grad_tuple) = exp(pre_mult_a + da_infsum);
    }

    if (calc_b) {
      T_vec pre_mult_b = (log_z + sum_log_a) - (log_b.array() + sum_log_b);
      std::get<1>(grad_tuple) = -exp(pre_mult_b + db_infsum);
    }
  }

  if (calc_z) {
    std::get<2>(grad_tuple) = exp(sum_log_a - sum_log_b)
                                * hypergeometric_pFq(ap1, bp1, z);
  }
}
}  // namespace internal

/**
 * Wrapper function for calculating gradients wrt all parameters
 * 
 * @tparam Ta Type of a
 * @tparam Tb Type of b
 * @tparam Tz Type of z
 * @param[in] grad_a Vector to evaluate gradients wrt a into
 * @param[in] grad_b Vector to evaluate gradients wrt b into
 * @param[in] grad_z Scalar to evaluate gradient wrt z into
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Ta, typename Tb, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr>
void grad_pFq(plain_type_t<Ta>& grad_a, plain_type_t<Tb>& grad_b,
              plain_type_t<Tz>& grad_z, const Ta& a, const Tb& b, const Tz& z,
              double precision = 1e-14, int max_steps = 1e6) {
  internal::grad_pFq_impl<true, true, true>(
    std::forward_as_tuple(grad_a, grad_b, grad_z),
    a, b, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt only a
 * 
 * @tparam Ta Type of a
 * @tparam Tb Type of b
 * @tparam Tz Type of z
 * @param[in] grad_a Vector to evaluate gradients wrt a into
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Ta, typename Tb, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr>
void grad_pFq_p(plain_type_t<Ta>& grad_a, const Ta& a, const Tb& b,
                const Tz& z, double precision = 1e-14,
                int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Ta>;
  internal::grad_pFq_impl<true, false, false>(
    std::forward_as_tuple(grad_a, plain_type_t<Ta>{}, scalar_type_t<Ta>{}),
    a, b, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt only b
 * 
 * @tparam Ta Type of a
 * @tparam Tb Type of b
 * @tparam Tz Type of z
 * @param[in] grad_b Vector to evaluate gradients wrt b into
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Ta, typename Tb, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr>
void grad_pFq_q(plain_type_t<Tb>& grad_b, const Ta& a, const Tb& b,
                const Tz& z, double precision = 1e-14,
                int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Tb>;
  internal::grad_pFq_impl<false, true, false>(
    std::forward_as_tuple(plain_type_t<Tb>{}, grad_b, scalar_type_t<Tb>{}),
    a, b, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt only z
 * 
 * @tparam Ta Type of a
 * @tparam Tb Type of b
 * @tparam Tz Type of z
 * @param[in] grad_z Scalar to evaluate gradient wrt z into
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Ta, typename Tb, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr>
void grad_pFq_z(plain_type_t<Tz>& grad_z, const Ta& a, const Tb& b,
                const Tz& z, double precision = 1e-14, int max_steps = 1e6) {
  internal::grad_pFq_impl<false, false, true>(
    std::forward_as_tuple(Eigen::Matrix<Tz, -1, 1>{},
                          Eigen::Matrix<Tz, -1, 1>{}, grad_z),
    a, b, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt a and b
 * 
 * @tparam Ta Type of a
 * @tparam Tb Type of b
 * @tparam Tz Type of z
 * @param[in] grad_a Vector to evaluate gradients wrt a into
 * @param[in] grad_b Vector to evaluate gradients wrt b into
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Ta, typename Tb, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr>
void grad_pFq_pq(plain_type_t<Ta>& grad_a, plain_type_t<Tb>& grad_b,
                 const Ta& a, const Tb& b, const Tz& z,
                 double precision = 1e-14, int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Ta>;
  internal::grad_pFq_impl<true, true, false>(
    std::forward_as_tuple(grad_a, grad_b, scalar_type_t<Ta>{}),
    a, b, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt a and z
 * 
 * @tparam Ta Type of a
 * @tparam Tb Type of b
 * @tparam Tz Type of z
 * @param[in] grad_a Vector to evaluate gradients wrt a into
 * @param[in] grad_z Scalar to evaluate gradient wrt z into
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Ta, typename Tb, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr>
void grad_pFq_pz(plain_type_t<Ta>& grad_a, plain_type_t<Tz>& grad_z,
                 const Ta& a, const Tb& b, const Tz& z,
                 double precision = 1e-14, int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Ta>;
  internal::grad_pFq_impl<true, false, true>(
    std::forward_as_tuple(grad_a, plain_type_t<Ta>{}, grad_z),
    a, b, z, precision, max_steps);
}

/**
 * Wrapper function for calculating gradients wrt b and z
 * 
 * @tparam Ta Type of a
 * @tparam Tb Type of b
 * @tparam Tz Type of z
 * @param[in] grad_b Vector to evaluate gradients wrt b into
 * @param[in] grad_z Scalar to evaluate gradient wrt z into
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 */
template <typename Ta, typename Tb, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr>
void grad_pFq_qz(plain_type_t<Tb>& grad_b, plain_type_t<Tz>& grad_z,
                 const Ta& a, const Tb& b, const Tz& z,
                 double precision = 1e-14, int max_steps = 1e6) {
  using scalar_t = scalar_type_t<Tb>;
  internal::grad_pFq_impl<false, true, true>(
    std::forward_as_tuple(plain_type_t<Tb>{}, grad_b, grad_z),
    a, b, z, precision, max_steps);
}

}  // namespace math
}  // namespace stan
#endif
