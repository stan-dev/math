#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_SERIAL_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_SERIAL_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/prim/mat/functor/map_rect_reduce.hpp>
#include <stan/math/prim/mat/functor/map_rect_combine.hpp>

#include <vector>

namespace stan {
namespace math {
namespace internal {

template <int call_id, typename F, typename T_shared_param,
          typename T_job_param>
Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type,
              Eigen::Dynamic, 1>
map_rect_serial(
    const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
    const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
        job_params,
    const std::vector<std::vector<double>>& x_r,
    const std::vector<std::vector<int>>& x_i, std::ostream* msgs = 0) {
  typedef map_rect_reduce<F, T_shared_param, T_job_param> ReduceF;
  typedef map_rect_combine<F, T_shared_param, T_job_param> CombineF;

  const int num_jobs = job_params.size();
  const vector_d shared_params_dbl = value_of(shared_params);

  matrix_d world_output(0, 0);
  std::vector<int> world_f_out(num_jobs, -1);

  for (int offset = 0, i = 0; i < num_jobs; offset += world_f_out[i], ++i) {
    const matrix_d job_output = ReduceF()(
        shared_params_dbl, value_of(job_params[i]), x_r[i], x_i[i], msgs);
    world_f_out[i] = job_output.cols();

    if (i == 0)
      world_output.resize(job_output.rows(), num_jobs * world_f_out[i]);

    if (world_output.cols() < offset + world_f_out[i])
      world_output.conservativeResize(Eigen::NoChange,
                                      2 * (offset + world_f_out[i]));

    world_output.block(0, offset, world_output.rows(), world_f_out[i])
        = job_output;
  }

  CombineF combine(shared_params, job_params);
  return combine(world_output, world_f_out);
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
