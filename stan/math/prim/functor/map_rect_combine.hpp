#ifndef STAN_MATH_PRIM_FUNCTOR_MAP_RECT_COMBINE_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAP_RECT_COMBINE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

#include <vector>

namespace stan {
namespace math {
namespace internal {

/* Template class for the combine step of map_rect which implements
 * the CombineF concept. The concept requires that
 *
 * - A nullary constructor creates a null combiner (used on the
 *    children)
 *
 * - A constructor which takes the shared and job-specific parameters
 *    (used on the root/main process)
 *
 * - Provides an operator() which takes as two arguments: (i) the
 *   function outputs as a ragged matrix and (ii) as second argument
 *   the output sizes of each function evaluation.
 *
 * This functor inserts the concatenated outputs of all reduce
 * operations into the autodiff stack. The concatenated results are
 * stored in a double only matrix and is ragged according to the
 * output sizes of each job.
 *
 * @tparam F type of user functor
 * @tparam T_shared_param type of shared parameters
 * @tparam T_job_param type of job specific parameters
 */
template <typename F, typename T_shared_param, typename T_job_param,
          require_eigen_col_vector_t<T_shared_param>* = nullptr>
class map_rect_combine {
  using ops_partials_t = internal::partials_propagator<
      return_type_t<T_shared_param, T_job_param>, void, T_shared_param,
      Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>;
  std::vector<ops_partials_t> ops_partials_;

  const std::size_t num_shared_operands_;
  const std::size_t num_job_operands_;

 public:
  using result_t = Eigen::Matrix<return_type_t<T_shared_param, T_job_param>,
                                 Eigen::Dynamic, 1>;

  map_rect_combine()
      : ops_partials_(), num_shared_operands_(0), num_job_operands_(0) {}
  map_rect_combine(
      const T_shared_param& shared_params,
      const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
          job_params)
      : ops_partials_(),
        num_shared_operands_(shared_params.rows()),
        num_job_operands_(dims(job_params)[1]) {
    ops_partials_.reserve(job_params.size());
    for (const auto& job_param : job_params) {
      ops_partials_.emplace_back(shared_params, job_param);
    }
  }

  result_t operator()(const matrix_d& world_result,
                      const std::vector<int>& world_f_out) {
    const std::size_t num_jobs = world_f_out.size();
    const std::size_t offset_job_params
        = is_constant_all<T_shared_param>::value ? 1 : 1 + num_shared_operands_;
    const std::size_t size_world_f_out = sum(world_f_out);

    result_t out(size_world_f_out);

    for (std::size_t i = 0, ij = 0; i != num_jobs; ++i) {
      for (int j = 0; j != world_f_out[i]; ++j, ++ij) {
        if (!is_constant_all<T_shared_param>::value) {
          edge<0>(ops_partials_[i]).partials_
              = world_result.block(1, ij, num_shared_operands_, 1);
        }

        if (!is_constant_all<T_job_param>::value) {
          edge<1>(ops_partials_[i]).partials_
              = world_result.block(offset_job_params, ij, num_job_operands_, 1);
        }

        out(ij) = ops_partials_[i].build(world_result(0, ij));
      }
    }

    return out;
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
