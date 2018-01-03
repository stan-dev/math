#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_REDUCE_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_REDUCE_HPP

namespace stan {
namespace math {

template <typename F, typename T_shared_param, typename T_job_param>
class map_rect_reduce {};

template <typename F>
class map_rect_reduce<F, double, double> {
 public:
  matrix_d operator()(const vector_d& shared_params,
                      const vector_d& job_specific_params,
                      const std::vector<double>& x_r,
                      const std::vector<int>& x_i) const {
    const F f;
    const vector_d out = f(shared_params, job_specific_params, x_r, x_i, 0);
    return out.transpose();
  }
};

}  // namespace math
}  // namespace stan

#endif
