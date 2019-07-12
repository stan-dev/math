#ifndef STAN_MATH_REV_MAT_FUNCTOR_MAP_RECT_CONCURRENT_TBB_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_MAP_RECT_CONCURRENT_TBB_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/functor/map_rect_concurrent.hpp>
#include <stan/math/prim/mat/functor/map_rect_reduce.hpp>
#include <stan/math/prim/mat/functor/map_rect_combine.hpp>
#include <stan/math/prim/scal/functor/parallel_for_each.hpp>
#include <stan/math/rev/core/chainablestack.hpp>

#include <vector>

namespace stan {
namespace math {
namespace internal {

template <int call_id, typename F, typename T_shared_param,
          typename T_job_param>
Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type,
              Eigen::Dynamic, 1>
map_rect_concurrent_tbb(
    const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
    const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
        job_params,
    const std::vector<std::vector<double>>& x_r,
    const std::vector<std::vector<int>>& x_i, std::ostream* msgs) {
  typedef map_rect_reduce<F, T_shared_param, T_job_param> ReduceF;
  typedef map_rect_combine<F, T_shared_param, T_job_param> CombineF;

  const int num_jobs = job_params.size();
  const vector_d shared_params_dbl = value_of(shared_params);
  typedef boost::counting_iterator<std::size_t> count_iter;

  std::vector<matrix_d> results = parallel_map(
      count_iter(0), count_iter(num_jobs), [&](std::size_t job) -> matrix_d {
        return ReduceF()(shared_params_dbl, value_of(job_params[job]), x_r[job],
                         x_i[job], msgs);
      });

  // collect results
  std::vector<int> world_f_out;
  std::size_t num_outputs = 0;
  for (std::size_t i = 0; i < num_jobs; ++i) {
    world_f_out.push_back(results[i].cols());
    num_outputs += results[i].cols();
  }

  world_f_out.reserve(num_jobs);

  matrix_d world_output(results[0].rows(), num_outputs);

  int offset = 0;
  for (const auto& job_result : results) {
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
