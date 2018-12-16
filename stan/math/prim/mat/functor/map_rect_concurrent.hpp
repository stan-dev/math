#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_CONCURRENT_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_CONCURRENT_HPP

#include <stan/math/parallel/for_each.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/prim/mat/functor/map_rect_reduce.hpp>
#include <stan/math/prim/mat/functor/map_rect_combine.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iterator/counting_iterator.hpp>

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
    const std::vector<std::vector<int>>& x_i, std::ostream* msgs = nullptr) {
  typedef map_rect_reduce<F, T_shared_param, T_job_param> ReduceF;
  typedef map_rect_combine<F, T_shared_param, T_job_param> CombineF;
  typedef boost::counting_iterator<int> count_iter;

#ifdef STAN_THREADS
  constexpr std_par::execution::parallel_unsequenced_policy exec_policy
      = std_par::execution::par_unseq;
#else
  constexpr std_par::execution::sequenced_policy exec_policy
      = std_par::execution::seq;
#endif
  using std_par::for_each;

  const int num_jobs = job_params.size();
  const vector_d shared_params_dbl = value_of(shared_params);

  std::vector<int> world_f_out(num_jobs);
  std::vector<matrix_d> world_job_output(num_jobs);

  for_each(exec_policy, count_iter(0), count_iter(num_jobs), [&](int i) {
    world_job_output[i] = ReduceF()(shared_params_dbl, value_of(job_params[i]),
                                    x_r[i], x_i[i], msgs);
    world_f_out[i] = world_job_output[i].cols();
  });

  // collect results
  const int num_world_output
      = std::accumulate(world_f_out.begin(), world_f_out.end(), 0);
  matrix_d world_output(world_job_output[0].rows(), num_world_output);

  int offset = 0;
  for (const auto& job_result : world_job_output) {
    const int num_job_outputs = job_result.cols();

    world_output.block(0, offset, world_output.rows(), num_job_outputs)
        = job_result;

    offset += num_job_outputs;
  }
  CombineF combine(shared_params, job_params);
  return combine(world_output, world_f_out);
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
