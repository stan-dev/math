#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MPI_MAP_RECT_COMBINE_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MPI_MAP_RECT_COMBINE_HPP

#include <boost/mpi.hpp>
#include <stan/math/prim/mat/functor/map_rect_combine.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename F, typename T_shared_param, typename T_job_param>
class mpi_map_rect_combine {
  boost::mpi::communicator world_;
  const std::size_t rank_ = world_.rank();
  const std::size_t world_size_ = world_.size();

  typedef map_rect_combine<F, T_shared_param, T_job_param> combine_t;
  combine_t combine_;

 public:
  typedef typename combine_t::result_type result_type;

  mpi_map_rect_combine() : combine_() {}
  mpi_map_rect_combine(
      const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
      const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>&
          job_params)
      : combine_(shared_params, job_params) {}

  result_type operator()(const matrix_d& local_result,
                         const std::vector<int>& world_f_out,
                         const std::vector<int>& job_chunks) {
    const std::size_t num_output_size_per_job = local_result.rows();
    const std::size_t size_world_f_out = sum(world_f_out);

    matrix_d world_result(num_output_size_per_job, size_world_f_out);

    std::vector<int> chunks_result(world_size_, 0);
    for (std::size_t i = 0, ij = 0; i != world_size_; ++i)
      for (std::size_t j = 0; j != job_chunks[i]; ++j, ++ij)
        chunks_result[i] += world_f_out[ij] * num_output_size_per_job;

    // collect results
    boost::mpi::gatherv(world_, local_result.data(), chunks_result[rank_],
                        world_result.data(), chunks_result, 0);

    // only add result to AD tree on root
    if (rank_ != 0)
      return result_type();

    // check if any of the workers flagged an error
    for (std::size_t i = 0, j = 0; i < world_f_out.size(); i++) {
      if (unlikely(world_result(0, j) == std::numeric_limits<double>::max()))
        throw std::domain_error("Error.");
      j += world_f_out[i];
    }

    return combine_(world_result, world_f_out);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
