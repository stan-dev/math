#ifndef STAN_MATH_REV_MAT_FUNCTOR_MAP_RECT_CONCURRENT_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_MAP_RECT_CONCURRENT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/functor/map_rect_concurrent.hpp>
#include <stan/math/prim/mat/functor/map_rect_reduce.hpp>
#include <stan/math/prim/mat/functor/map_rect_combine.hpp>
#include <stan/math/rev/core/chainablestack.hpp>

#ifdef STAN_THREADS
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif

#include <algorithm>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <int call_id, typename F, typename T_shared_param,
          typename T_job_param>
Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type,
              Eigen::Dynamic, 1>
map_rect_concurrent(
    const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
    const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
        job_params,
    const std::vector<std::vector<double>>& x_r,
    const std::vector<std::vector<int>>& x_i, std::ostream* msgs) {
  using ReduceF = map_rect_reduce<F, T_shared_param, T_job_param>;
  using CombineF = map_rect_combine<F, T_shared_param, T_job_param>;

  const int num_jobs = job_params.size();
  const vector_d shared_params_dbl = value_of(shared_params);
  std::vector<matrix_d> job_output(num_jobs);
  std::vector<int> world_f_out(num_jobs, 0);

  auto execute_chunk = [&](std::size_t start, std::size_t end) -> void {
    for (std::size_t i = start; i != end; ++i) {
      job_output[i] = ReduceF()(shared_params_dbl, value_of(job_params[i]),
                                x_r[i], x_i[i], msgs);
      world_f_out[i] = job_output[i].cols();
    }
  };

#ifdef STAN_THREADS
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, num_jobs),
                    [&](const tbb::blocked_range<size_t>& r) {
                      execute_chunk(r.begin(), r.end());
                    });
#else
  execute_chunk(0, num_jobs);
#endif

  // collect results
  const int num_world_output
      = std::accumulate(world_f_out.begin(), world_f_out.end(), 0);
  matrix_d world_output(job_output[0].rows(), num_world_output);

  int offset = 0;
  for (const auto& job : job_output) {
    const int num_job_outputs = job.cols();

    world_output.block(0, offset, world_output.rows(), num_job_outputs) = job;

    offset += num_job_outputs;
  }
  CombineF combine(shared_params, job_params);
  return combine(world_output, world_f_out);
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
