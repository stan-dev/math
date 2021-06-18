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
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename Tp, typename Tq, typename Tz,
          require_stan_scalar_t<Tz>* = nullptr,
          require_all_eigen_vector_t<Tp, Tq>* = nullptr>
void grad_pFq(plain_type_t<Tp>& grad_p, plain_type_t<Tq>& grad_q,
              Tz& grad_z,
              const Tp& p, const Tq& q, const Tz& z,
              double precision = 1e-14, int max_steps = 1e6) {
    using std::max;
    using scalar_t = scalar_type_t<return_type_t<Tp, Tq, Tz>>;
    using Tp_plain = plain_type_t<Tp>;
    using Tq_plain = plain_type_t<Tq>;

    Tp_plain dp_infsum = Tp_plain::Constant(p.size(), NEGATIVE_INFTY);
    Tq_plain dq_infsum = Tq_plain::Constant(q.size(), NEGATIVE_INFTY);

    scalar_t log_z = log(z);
    double log_precision = log(precision);
  
    int m = 0;
    scalar_t outer_diff = 0;
    
    while((outer_diff > log_precision) & (m < max_steps)) {
      Tp_plain dp_iter_m = Tp_plain::Constant(p.size(),NEGATIVE_INFTY);
      Tq_plain dq_iter_m = Tq_plain::Constant(q.size(),NEGATIVE_INFTY);
      scalar_t log_phammer_1m = log_rising_factorial(1,m);
      scalar_t lgamma_mp1 = lgamma(m+1);
        int n = 0;
        scalar_t inner_diff = 0;
        while((inner_diff > log_precision) & (n < max_steps)) {
          scalar_t term1_mn = (m+n)*log_z + sum(log_rising_factorial(p.array()+1,m+n)) + log_phammer_1m + log_rising_factorial(1,n);
          scalar_t term2_mn = lgamma_mp1 + lgamma(n+1) + sum(log_rising_factorial(q.array()+1,m+n)) + log_rising_factorial(2,m+n);

          Tp_plain dp_mn = (term1_mn+log_rising_factorial(p, n).array()) - (term2_mn+log_rising_factorial(p.array()+1,n).array());
          Tq_plain dq_mn = (term1_mn+log_rising_factorial(q, n).array()) - (term2_mn+log_rising_factorial(q.array()+1,n).array());

          dp_iter_m = dp_iter_m.binaryExpr(dp_mn,[&](auto&& a, auto&& b){ return log_sum_exp(a, b); });
          dq_iter_m = dq_iter_m.binaryExpr(dq_mn,[&](auto&& a, auto&& b){ return log_sum_exp(a, b); });

          inner_diff = max(log_sum_exp(dp_mn), log_sum_exp(dq_mn));
           n += 1;
        }
        
      dp_infsum = dp_infsum.binaryExpr(dp_iter_m,[&](auto&& a, auto&& b){ return log_sum_exp(a, b); });
      dq_infsum = dq_infsum.binaryExpr(dq_iter_m,[&](auto&& a, auto&& b){ return log_sum_exp(a, b); });

      outer_diff = max(log_sum_exp(dp_iter_m), log_sum_exp(dq_iter_m));
      m += 1;
    }

    Eigen::VectorXi ind_vector = Eigen::VectorXi::LinSpaced(p.size(),0,p.size()-1);
    Tp_plain log_p = log(p);
    Tq_plain log_q = log(q);
    scalar_t sum_log_p = sum(log_p);
    scalar_t sum_log_q = sum(log_q);
    Tp_plain pre_mult_p = (log_z + ind_vector.unaryExpr([&log_p, &sum_log_p](int i){ return sum_log_p - log_p[i]; }).array() - sum_log_q).matrix();
    Tq_plain pre_mult_q = (log_z + sum_log_p) - (log_q.array() + sum_log_q);
  
    grad_p = exp(pre_mult_p.array() + dp_infsum.array()).matrix();
    grad_q = -exp(pre_mult_q.array() + dq_infsum.array()).matrix();
    grad_z = (prod(p) / prod(q)) * hypergeometric_pFq((p.array() + 1).matrix(), (q.array() + 1).matrix(), z);
}
}  // namespace math
}  // namespace stan
#endif
