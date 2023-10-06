#ifndef STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP
#define STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/max.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/log_sum_exp_signed.hpp>

namespace stan {
namespace math {
namespace internal {
/**
 * Wrapper function for calculating gradients for the generalized
 * hypergeometric function. The function always returns a tuple with
 * three elements (gradients wrt a, b, and z, respectively), but the
 * elements will only be defined/calculated when the respective parameter
 * is not a primitive type.
 *
 * @tparam Ta Eigen type with either one row or column at compile time
 * @tparam Tb Eigen type with either one row or column at compile time
 * @tparam Tz Scalar type
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @param[in] precision Convergence criteria for infinite sum
 * @param[in] outer_steps Maximum number of iterations for outer_steps
 * @param[in] inner_steps Maximum number of iterations for inner_steps
 * @return Tuple of gradients
 */
template <bool CalcA, bool CalcB, bool CalcZ,
          typename TpFq, typename Ta, typename Tb, typename Tz>
auto grad_pFq_impl(const TpFq& pfq_val, const Ta& a_val, const Tb& b_val, const Tz& z_val,
                double precision = 1e-16,
                int max_steps = 1e8) {
  using std::max;
  using T_partials = return_type_t<Ta, Tb, Tz>;
  using Ta_partials = return_type_t<Ta>;
  using Tb_partials = return_type_t<Tb>;
  using Tz_partials = Tz;
  using T_ArrayPartials = promote_scalar_t<T_partials, Eigen::Array<double, -1, 1>>;
  using IntArray = Eigen::Array<int, -1, 1>;
  using Ta_PartialsArray = promote_scalar_t<Ta_partials, Eigen::Array<double, -1, 1>>;
  using Tb_PartialsArray = promote_scalar_t<Tb_partials, Eigen::Array<double, -1, 1>>;

  IntArray a_signs = sign(value_of_rec(a_val));
  IntArray b_signs = sign(value_of_rec(b_val));
  int z_sign = sign(value_of_rec(z_val));

  IntArray a_k_signs = a_signs;
  IntArray b_k_signs = b_signs;
  int z_k_sign = 1;

  Ta_PartialsArray a_k = a_val;
  Tb_PartialsArray b_k = b_val;
  Tz_partials log_z_pow_k = 0.0;
  Tz_partials log_z_k = log(abs(z_val));

  Ta_PartialsArray log_phammer_a_k = Ta_PartialsArray::Zero(a_val.size());
  Tb_PartialsArray log_phammer_b_k = Tb_PartialsArray::Zero(b_val.size());
  IntArray phammer_a_k_signs = IntArray::Ones(a_val.size());
  IntArray phammer_b_k_signs = IntArray::Ones(b_val.size());

  Ta_PartialsArray digamma_a = inv(a_k);
  Tb_PartialsArray digamma_b = inv(b_k);
  Ta_PartialsArray digamma_a_k = digamma_a;
  Tb_PartialsArray digamma_b_k = digamma_b;
  Ta_PartialsArray log_digamma_a_k = log(abs(digamma_a));
  Tb_PartialsArray log_digamma_b_k = log(abs(digamma_b));
  IntArray digamma_a_k_signs = a_k_signs;
  IntArray digamma_b_k_signs = b_k_signs;

  T_ArrayPartials a_log_infsum = T_ArrayPartials::Constant(a_val.size(), NEGATIVE_INFTY);
  T_ArrayPartials b_log_infsum = T_ArrayPartials::Constant(b_val.size(), NEGATIVE_INFTY);
  IntArray a_infsum_sign = IntArray::Ones(a_val.size());
  IntArray b_infsum_sign = IntArray::Ones(b_val.size());

  T_ArrayPartials a_grad = T_ArrayPartials::Constant(a_val.size(), NEGATIVE_INFTY);
  T_ArrayPartials b_grad = T_ArrayPartials::Constant(a_val.size(), NEGATIVE_INFTY);
  IntArray a_grad_signs = IntArray::Ones(a_val.size());
  IntArray b_grad_signs = IntArray::Ones(a_val.size());

  int k = 0;
  T_partials curr_log_prec = 1.0;
  double log_factorial_k = 0;
  while (k <= max_steps && curr_log_prec > log(precision)) {
    int base_sign = z_k_sign * phammer_a_k_signs.prod() * phammer_b_k_signs.prod();
    T_partials log_base = (sum(log_phammer_a_k) + log_z_pow_k)
                            - (log_factorial_k + sum(log_phammer_b_k));

    if (CalcA) {
      a_grad = log_digamma_a_k + log_base;
      a_grad_signs = base_sign * digamma_a_k_signs;

      for (int i = 0; i < a_val.size(); i++) {
        std::forward_as_tuple(a_log_infsum(i), a_infsum_sign(i))
          = log_sum_exp_signed(a_log_infsum(i), a_infsum_sign(i),
                                a_grad(i), a_grad_signs(i));
      }
    }

    if (CalcB) {
      b_grad = log_digamma_b_k + log_base;
      b_grad_signs = base_sign * digamma_b_k_signs;

      for (int i = 0; i < b_val.size(); i++) {
        std::forward_as_tuple(b_log_infsum(i), b_infsum_sign(i))
          = log_sum_exp_signed(b_log_infsum(i), b_infsum_sign(i),
                                b_grad(i), b_grad_signs(i));
      }
    }

    curr_log_prec = max(a_grad.maxCoeff(), b_grad.maxCoeff());
    log_phammer_a_k += log(abs(a_k));
    log_phammer_b_k += log(abs(b_k));
    phammer_a_k_signs *= a_k_signs;
    phammer_b_k_signs *= b_k_signs;

    log_z_pow_k += log_z_k;
    z_k_sign *= z_sign;

    digamma_a_k += select(a_k == 0.0, 0.0, inv(a_k));
    digamma_b_k += select(b_k == 0.0, 0.0, inv(b_k));
    log_digamma_a_k = log(abs(digamma_a_k));
    log_digamma_b_k = log(abs(digamma_b_k));
    digamma_a_k_signs = sign(value_of_rec(digamma_a_k));
    digamma_b_k_signs = sign(value_of_rec(digamma_b_k));

    a_k += 1.0;
    b_k += 1.0;
    a_k_signs = sign(value_of_rec(a_k));
    b_k_signs = sign(value_of_rec(b_k));

    k += 1;
    log_factorial_k += log(k);
  }
  std::tuple<Eigen::Matrix<T_partials, -1, 1>,
              Eigen::Matrix<T_partials, -1, 1>,
              T_partials> ret_tuple;

  if (CalcA) {
    std::get<0>(ret_tuple) = a_infsum_sign * exp(a_log_infsum) - digamma_a * pfq_val;
  }
  if (CalcB) {
    std::get<1>(ret_tuple) = digamma_b * pfq_val - b_infsum_sign * exp(b_log_infsum);
  }
  if (CalcZ) {
    T_partials pfq_p1_val = hypergeometric_pFq((a_val + 1).matrix(), (b_val + 1).matrix(), z_val);
    std::get<2>(ret_tuple) = prod(a_val) / prod(b_val) * pfq_p1_val;
  }
  return ret_tuple;
}
}


template <typename TpFq, typename Ta, typename Tb, typename Tz>
auto grad_pFq(const TpFq& pfq_val, const Ta& a, const Tb& b, const Tz& z,
                double precision = 1e-16,
                int max_steps = 1e8) {
  return internal::grad_pFq_impl<
    !is_constant_all<Ta>::value,
    !is_constant_all<Tb>::value,
    !is_constant_all<Tz>::value>(
      pfq_val,
      to_ref(as_value_column_array_or_scalar(a)),
      to_ref(as_value_column_array_or_scalar(b)),
      value_of(z)
    );
}

}  // namespace math
}  // namespace stan
#endif
