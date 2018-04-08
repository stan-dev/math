#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_CONCURRENT_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_CONCURRENT_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/prim/mat/functor/map_rect_reduce.hpp>
#include <stan/math/prim/mat/functor/map_rect_combine.hpp>

#include <vector>
#include <thread>
#include <future>
#include <cstdlib>

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
    const std::vector<std::vector<int>>& x_i, std::ostream* msgs = 0) {
  typedef map_rect_reduce<F, T_shared_param, T_job_param> ReduceF;
  typedef map_rect_combine<F, T_shared_param, T_job_param> CombineF;

  const int num_jobs = job_params.size();
  const vector_d shared_params_dbl = value_of(shared_params);
  std::vector<std::future<std::vector<matrix_d>>> futures;

  auto chunk_job = [&](int start, int end) {
    const int size = end - start;
    std::vector<matrix_d> chunk_f_out;
    chunk_f_out.reserve(size);
    for (int i = start; i != end; i++)
      chunk_f_out.push_back(ReduceF()(
          shared_params_dbl, value_of(job_params[i]), x_r[i], x_i[i], msgs));
    return chunk_f_out;
  };

  int num_threads = 1;

#ifdef STAN_THREADS
  const char* env_stan_num_threads = std::getenv("STAN_NUM_THREADS");
  if (env_stan_num_threads != nullptr) {
    const int env_num_threads = std::atoi(env_stan_num_threads);
    if (env_num_threads > 0)
      num_threads = env_num_threads;
    else if (env_num_threads == -1)
      num_threads = std::thread::hardware_concurrency();
    // anything else will use 1 thread.
  }
#endif

  const int num_jobs_per_thread = num_jobs / num_threads;
  int num_jobs_per_thread_remainder = num_jobs % num_threads;
  const int num_jobs_min = num_threads - num_jobs_per_thread_remainder;
  int job_start = 0;
  for (int j = 0; j < num_threads; j++) {
    int job_end = job_start + num_jobs_per_thread;
    // the excess jobs are assigned to the last processes
    if (j >= num_jobs_min && num_jobs_per_thread_remainder > 0) {
      job_end++;
      num_jobs_per_thread_remainder--;
    }
    if (j == num_threads - 1)
      job_end = num_jobs;
    // we only defer the first chunk such that the main thread is not
    // blocking if we use threading
    futures.emplace_back(std::async(j == 0 ? std::launch::deferred :
#ifndef STAN_THREADS
                                           std::launch::deferred,
#else
                                           std::launch::async,
#endif
                                    chunk_job, job_start, job_end));
    job_start = job_end;
  }

  // collect results
  std::vector<int> world_f_out;
  world_f_out.reserve(num_jobs);
  matrix_d world_output(0, 0);

  for (int j = 0, offset = 0; j < num_threads; j++) {
    const std::vector<matrix_d>& chunk_result = futures[j].get();
    if (j == 0)
      world_output.resize(chunk_result[0].rows(),
                          num_jobs * chunk_result[0].cols());

    for (const auto& job_result : chunk_result) {
      const int num_job_outputs = job_result.cols();
      world_f_out.push_back(num_job_outputs);

      if (world_output.cols() < offset + num_job_outputs)
        world_output.conservativeResize(Eigen::NoChange,
                                        2 * (offset + num_job_outputs));

      world_output.block(0, offset, world_output.rows(), num_job_outputs)
          = job_result;

      offset += num_job_outputs;
    }
  }

  CombineF combine(shared_params, job_params);
  return combine(world_output, world_f_out);
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
