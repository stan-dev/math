#ifndef STAN_MATH_REV_ARR_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_REV_ARR_FUN_LOG_SUM_EXP_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/calculate_chain.hpp>
#include <stan/math/prim/arr/fun/log_sum_exp.hpp>
#include <vector>
#include <limits>
#include <algorithm>

namespace stan {
namespace math {

namespace internal {
inline double log_sum_exp_as_double(const std::vector<var>& x) {
  using std::exp;
  using std::log;
  using std::numeric_limits;
  double max_val = std::max_element(x.begin(), x.end())->val();
  double sum = std::accumulate(x.begin(), x.end(), 0.0,
  [&max_val](auto& acc, auto&& x_i) {
    return acc + exp(x_i.val() - max_val);
  });
  return max_val + log(sum);
}

class log_sum_exp_vector_vari : public op_vector_vari {
 public:
  explicit log_sum_exp_vector_vari(const std::vector<var>& x)
      : op_vector_vari(log_sum_exp_as_double(x), x) {}
  void chain() {
    for (size_t i = 0; i < size_; ++i) {
      vis_[i]->adj_ += adj_ * calculate_chain(vis_[i]->val_, val_);
    }
  }
};
}  // namespace internal

/**
 * Returns the log sum of exponentials.
 */
inline var log_sum_exp(const std::vector<var>& x) {
  return var(new internal::log_sum_exp_vector_vari(x));
}

}  // namespace math
}  // namespace stan
#endif
