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

  T_ArrayPartials a_val = as_value_column_array_or_scalar(a);
  T_ArrayPartials b_val = as_value_column_array_or_scalar(b);
  const auto& z_val = value_of(z);
  T_ArrayPartials a_k = a_val;
  T_ArrayPartials b_k = b_val;
  T_ArrayPartials phammer_a_k = T_ArrayPartials::Ones(a.size());
  T_ArrayPartials phammer_b_k = T_ArrayPartials::Ones(b.size());
  T_ArrayPartials digamma_a = inv(a_k);
  T_ArrayPartials digamma_b = inv(b_k);
  T_ArrayPartials digamma_a_k = digamma_a;
  T_ArrayPartials digamma_b_k = digamma_b;
  T_ArrayPartials a_infsum = T_ArrayPartials::Zero(a.size());
  T_ArrayPartials b_infsum = T_ArrayPartials::Zero(b.size());

  int k = 0;
  double curr_log_prec = 1.0;
  double curr_prec = 1.0;
  double factorial_k = 1;
  while (k <= max_steps && curr_prec > precision) {
    T_partials num = prod(phammer_a_k) * pow(z_val, k);
    T_partials den = factorial_k * prod(phammer_b_k);
    T_ArrayPartials a_grad = (digamma_a_k * num) / den;
    T_ArrayPartials b_grad = (digamma_b_k * num) / den;

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

    curr_prec = max(a_grad.abs().maxCoeff(), b_grad.abs().maxCoeff());
    a_infsum += a_grad;
    b_infsum += b_grad;


    phammer_a_k *= a_k;
    phammer_b_k *= b_k;
    digamma_a_k += select(a_k == 0.0, 0.0, inv(a_k));
    digamma_b_k += select(b_k == 0.0, 0.0, inv(b_k));
    a_k += 1.0;
    b_k += 1.0;
    k += 1;
    factorial_k *= k;
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
