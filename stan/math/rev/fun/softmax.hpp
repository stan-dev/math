#ifndef STAN_MATH_REV_FUN_SOFTMAX_HPP
#define STAN_MATH_REV_FUN_SOFTMAX_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the softmax of the specified Eigen vector.  Softmax is
 * guaranteed to return a simplex.
 *
 * @param alpha Unconstrained input vector.
 * @return Softmax of the input.
 * @throw std::domain_error If the input vector is size 0.
 */
template <typename Mat, require_rev_matrix_t<Mat>* = nullptr>
inline plain_type_t<Mat> softmax(const Mat& alpha) {
  using mat_plain = plain_type_t<Mat>;
  if (alpha.size() == 0) {
    return alpha;
  }
  arena_t<mat_plain> alpha_arena = alpha;
  arena_t<Eigen::VectorXd> res_val = softmax(value_of(alpha_arena));
  arena_t<mat_plain> res = res_val;

  reverse_pass_callback([res_val, res, alpha_arena]() mutable {
    const auto& res_adj = to_ref(res.adj());
    alpha_arena.adj()
        += -res_val * res_adj.dot(res_val) + res_val.cwiseProduct(res_adj);
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
