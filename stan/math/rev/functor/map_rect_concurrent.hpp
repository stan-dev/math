#ifndef STAN_MATH_REV_FUNCTOR_MAP_RECT_CONCURRENT_HPP
#define STAN_MATH_REV_FUNCTOR_MAP_RECT_CONCURRENT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/functor/map_rect_concurrent.hpp>
#include <stan/math/prim/functor/map_rect_reduce.hpp>
#include <stan/math/prim/functor/map_rect_combine.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/team_thread_pool.hpp>

//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range.h>

#include <algorithm>
#include <numeric>
#include <thread>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <int call_id, typename F, typename T_shared_param,
          typename T_job_param, require_eigen_col_vector_t<T_shared_param>*>
inline Eigen::Matrix<return_type_t<T_shared_param, T_job_param>, Eigen::Dynamic,
                     1>
map_rect_concurrent(
    const T_shared_param& shared_params,
    const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
        job_params,
    const std::vector<std::vector<double>>& x_r,
    const std::vector<std::vector<int>>& x_i, std::ostream* msgs) {
  using ReduceF
      = map_rect_reduce<F, scalar_type_t<T_shared_param>, T_job_param>;
  using CombineF = map_rect_combine<F, T_shared_param, T_job_param>;

  const std::size_t num_jobs = job_params.size();
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
  auto& pool = stan::math::TeamThreadPool::instance();

  // Total participants includes caller (tid=0).
  const std::size_t max_team = pool.team_size();
  const std::size_t n
      = std::min<std::size_t>(max_team, num_jobs == 0 ? 1u : num_jobs);

  if (n <= 1 || num_jobs <= 1) {
    execute_chunk(0, num_jobs);
  } else {
    pool.parallel_region(n, [&](std::size_t tid) {
      const std::size_t nj = num_jobs;
      const std::size_t b0 = (nj * tid) / n;
      const std::size_t b1 = (nj * (tid + 1)) / n;
      if (b0 < b1) {
        execute_chunk(b0, b1);
      }
    });
  }
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
