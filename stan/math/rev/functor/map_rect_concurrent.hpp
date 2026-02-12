#ifndef STAN_MATH_REV_FUNCTOR_MAP_RECT_CONCURRENT_HPP
#define STAN_MATH_REV_FUNCTOR_MAP_RECT_CONCURRENT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/functor/map_rect_reduce.hpp>
#include <stan/math/prim/functor/map_rect_combine.hpp>
#include <stan/math/rev/core/team_thread_pool.hpp>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>

namespace stan {
namespace math {
namespace internal {

// NOTE: Do NOT redeclare the default template argument here.
// The prim header declares it already; rev header provides the definition.
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
  if (num_jobs == 0) {
    return Eigen::Matrix<return_type_t<T_shared_param, T_job_param>,
                         Eigen::Dynamic, 1>();
  }

  if (x_r.size() != num_jobs || x_i.size() != num_jobs) {
    throw std::invalid_argument(
        "map_rect_concurrent: x_r and x_i must have size == job_params.size()");
  }

  const vector_d shared_params_dbl = value_of(shared_params);

  std::vector<matrix_d> job_output(num_jobs);
  std::vector<int> world_f_out(num_jobs, 0);

  // Completion flags to avoid using an uncomputed job_output[k]
  // (turns Eigen assert into a clear runtime error).
  std::vector<unsigned char> done(num_jobs, 0);

  auto execute_chunk = [&](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      job_output[i] = ReduceF()(shared_params_dbl, value_of(job_params[i]),
                                x_r[i], x_i[i], msgs);
      world_f_out[i] = static_cast<int>(job_output[i].cols());
      done[i] = 1;
    }
  };
#ifdef STAN_THREADS
  auto& pool = stan::math::TeamThreadPool::instance();

  // If we're already on a pool worker thread, run serially.
  // This avoids nested team scheduling issues.
  if (stan::math::TeamThreadPool::in_worker_thread()) {
    execute_chunk(0, num_jobs);
  } else {
    const std::size_t max_team = pool.team_size();
    std::size_t n = std::min<std::size_t>(max_team, num_jobs);
    if (n < 1)
      n = 1;

    if (n <= 1 || num_jobs <= 1) {
      execute_chunk(0, num_jobs);
    } else {
      pool.parallel_region(n, [&](std::size_t tid) {
        const std::size_t b0 = (num_jobs * tid) / n;
        const std::size_t b1 = (num_jobs * (tid + 1)) / n;
        if (b0 < b1)
          execute_chunk(b0, b1);
      });
    }
  }
#else
  execute_chunk(0, num_jobs);
#endif

  // Verify all jobs computed before stitching.
  for (std::size_t i = 0; i < num_jobs; ++i) {
    if (!done[i]) {
      throw std::runtime_error("map_rect_concurrent: job " + std::to_string(i)
                               + " was not computed (race/scheduling bug)");
    }
  }

  const int num_world_output
      = std::accumulate(world_f_out.begin(), world_f_out.end(), 0);

  // Determine row count from the first non-empty job output.
  int out_rows = 0;
  for (std::size_t i = 0; i < num_jobs; ++i) {
    if (job_output[i].size() > 0) {
      out_rows = static_cast<int>(job_output[i].rows());
      break;
    }
  }

  matrix_d world_output(out_rows, num_world_output);

  int offset = 0;
  for (std::size_t i = 0; i < num_jobs; ++i) {
    const auto& job = job_output[i];
    const int c = static_cast<int>(job.cols());
    if (c == 0)
      continue;

    if (job.rows() != out_rows) {
      throw std::runtime_error(
          "map_rect_concurrent: inconsistent job output rows; job "
          + std::to_string(i) + " rows=" + std::to_string(job.rows())
          + " expected=" + std::to_string(out_rows));
    }

    world_output.block(0, offset, out_rows, c) = job;
    offset += c;
  }

  CombineF combine(shared_params, job_params);
  return combine(world_output, world_f_out);
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
