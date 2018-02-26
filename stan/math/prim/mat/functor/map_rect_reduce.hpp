#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_REDUCE_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_REDUCE_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <vector>

namespace stan {
namespace math {
namespace internal {

template <typename F, typename T_shared_param, typename T_job_param>
class map_rect_reduce {};

template <typename F>
class map_rect_reduce<F, double, double> {
 public:
  matrix_d operator()(const vector_d& shared_params,
                      const vector_d& job_specific_params,
                      const std::vector<double>& x_r,
                      const std::vector<int>& x_i,
                      std::ostream* msgs = 0) const {
    const matrix_d out
        = F()(shared_params, job_specific_params, x_r, x_i, msgs).transpose();
    return out;
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
