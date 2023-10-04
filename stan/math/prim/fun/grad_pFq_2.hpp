#ifndef STAN_MATH_PRIM_FUN_GRAD_PFQ_2_HPP
#define STAN_MATH_PRIM_FUN_GRAD_PFQ_2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>

namespace stan {
namespace math {
namespace internal {
  template <typename T>
  auto bin_sign(const T& x) {
    decltype(auto) sign1 = sign(value_of_rec(x));
    return select(sign1 == 0, 1, sign1);
  }
}

template <typename Ta, typename Tb, typename Tz>
auto grad_pFq_2(const Ta& a, const Tb& b, const Tz& z,
              double precision = 1e-16,
              int max_steps = 1e8) {
  using T_partials = partials_return_t<Ta, Tb, Tz>;
  using Ta_partials = promote_scalar_t<T_partials, plain_type_t<Ta>>;
  using Tb_partials = promote_scalar_t<T_partials, plain_type_t<Tb>>;
  using Tz_partials = promote_scalar_t<T_partials, plain_type_t<Tz>>;
  using T_ArrayPartials = promote_scalar_t<T_partials, Eigen::Array<double, -1, 1>>;
  using T_return = std::tuple<Ta_partials, Tb_partials, Tz_partials>;
  using IntArray = Eigen::Array<int, -1, 1>;

  T_ArrayPartials a_val = as_value_column_array_or_scalar(a);
  T_ArrayPartials b_val = as_value_column_array_or_scalar(b);
  Tz_partials z_val = value_of(z);

  IntArray a_signs = internal::bin_sign(a_val);
  IntArray b_signs = internal::bin_sign(b_val);
  int z_sign = internal::bin_sign(z_val);

  IntArray a_k_signs = a_signs;
  IntArray b_k_signs = b_signs;
  int z_k_sign = 1;

  T_ArrayPartials a_k = a_val;
  T_ArrayPartials b_k = b_val;
  Tz_partials log_z_pow_k = 0.0;
  Tz_partials log_z_k = log(abs(z_val));

  T_ArrayPartials log_phammer_a_k = T_ArrayPartials::Zero(a.size());
  T_ArrayPartials log_phammer_b_k = T_ArrayPartials::Zero(b.size());
  IntArray phammer_a_k_signs = IntArray::Ones(a.size());
  IntArray phammer_b_k_signs = IntArray::Ones(b.size());

  T_ArrayPartials digamma_a = inv(a_k);
  T_ArrayPartials digamma_b = inv(b_k);
  T_ArrayPartials digamma_a_k = digamma_a;
  T_ArrayPartials digamma_b_k = digamma_b;
  T_ArrayPartials log_digamma_a_k = log(abs(digamma_a));
  T_ArrayPartials log_digamma_b_k = log(abs(digamma_b));
  IntArray digamma_a_k_signs = a_k_signs;
  IntArray digamma_b_k_signs = b_k_signs;
  T_ArrayPartials a_infsum = T_ArrayPartials::Zero(a.size());
  T_ArrayPartials b_infsum = T_ArrayPartials::Zero(b.size());

  int k = 0;
  double curr_log_prec = 1.0;
  double curr_prec = 1.0;
  double log_factorial_k = 0;
  while (k <= max_steps && curr_log_prec > log(precision)) {
    int base_sign = z_k_sign * phammer_a_k_signs.prod() * phammer_b_k_signs.prod();
    T_partials log_base = (sum(log_phammer_a_k) + log_z_pow_k)
                          - (log_factorial_k + sum(log_phammer_b_k));
    T_ArrayPartials a_grad = log_digamma_a_k + log_base;
    T_ArrayPartials b_grad = log_digamma_b_k + log_base;

    curr_log_prec = max(a_grad.maxCoeff(), b_grad.maxCoeff());
    a_infsum += exp(a_grad) * base_sign * digamma_a_k_signs;
    b_infsum += exp(b_grad) * base_sign * digamma_b_k_signs;

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
    digamma_a_k_signs = internal::bin_sign(digamma_a_k);
    digamma_b_k_signs = internal::bin_sign(digamma_b_k);

    a_k += 1.0;
    b_k += 1.0;
    a_k_signs = internal::bin_sign(a_k);
    b_k_signs = internal::bin_sign(b_k);

    k += 1;
    log_factorial_k += log(k);
  }

  T_partials pfq_val = hypergeometric_pFq(a_val, b_val, z_val);
  T_partials pfq_p1_val = hypergeometric_pFq(a_val + 1, b_val + 1, z_val);

  T_return ret_tuple;
  std::get<0>(ret_tuple) = a_infsum - digamma_a * pfq_val;
  std::get<1>(ret_tuple) = digamma_b * pfq_val - b_infsum;
  std::get<2>(ret_tuple) = prod(a_val) / prod(b_val) * pfq_p1_val;
  return ret_tuple;
}

}  // namespace math
}  // namespace stan
#endif



/*     std::cout << "k: " << k << "\n"
              << "num: " << num << "\n"
              << "den: " << den << "\n"
              << "phammer_a_k: " << phammer_a_k.matrix().transpose() << "\n"
              << "phammer_b_k: " << phammer_b_k.matrix().transpose() << "\n"
              << "digamma_a_k: " << digamma_a_k.matrix().transpose() << "\n"
              << "digamma_b_k: " << digamma_b_k.matrix().transpose() << "\n"
              << "a: " << a_grad.matrix().transpose() << "\n"
              << "b: " << b_grad.matrix().transpose() << "\n"
              << std::endl; */
