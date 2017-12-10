#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_MPI_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_MPI_HPP


#include <stan/math/prim/mat/functor/map_rect.hpp>
#include <stan/math/prim/mat/functor/mpi_parallel_call.hpp>

#include <vector>

namespace stan {
  namespace math {

    template <typename F, typename T_shared_param, typename T_job_param>
    class map_rect_reduce;

    template <typename F, typename T_shared_param, typename T_job_param>
    class map_rect_combine;

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
      mpi_map_rect_combine(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                           const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>& job_params)
        : combine_(shared_params, job_params) {}

      result_type operator()(const matrix_d& local_result, const std::vector<int>& world_f_out, const std::vector<int>& job_chunks) {

        const std::size_t num_output_size_per_job = local_result.rows();
        const std::size_t size_world_f_out = sum(world_f_out);

        matrix_d world_result(num_output_size_per_job, size_world_f_out);

        std::vector<int> chunks_result(world_size_,0);
        for(std::size_t i=0, ij=0; i != world_size_; ++i)
          for(std::size_t j=0; j != job_chunks[i]; ++j, ++ij)
            chunks_result[i] += world_f_out[ij] * num_output_size_per_job;

        // collect results
        boost::mpi::gatherv(world_, local_result.data(), chunks_result[rank_], world_result.data(), chunks_result, 0);

        // only add result to AD tree on root
        if(rank_ != 0)
          return(result_type());

        return(combine_(world_result, world_f_out));
      }
    };

    template <int call_id, typename F, typename T_shared_param, typename T_job_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1>
    map_rect_mpi(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                 const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>& job_params,
                 const std::vector<std::vector<double>>& x_r,
                 const std::vector<std::vector<int>>& x_i,
                 std::ostream* msgs = 0) {
      typedef map_rect_reduce<F, T_shared_param, T_job_param> ReduceF;
      typedef mpi_map_rect_combine<F, T_shared_param, T_job_param> CombineF;
      
      mpi_parallel_call<call_id,ReduceF,CombineF> job_chunk(shared_params, job_params, x_r, x_i);

      return job_chunk.reduce();
    }
  }
}


#define STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, SHARED, JOB) \
  namespace stan { namespace math { namespace internal {                \
                       typedef FUNCTOR mpi_mr_ ## CALLID ## _ ## SHARED ## _ ## JOB ## _; \
                       typedef map_rect_reduce<mpi_mr_ ## CALLID ## _ ## SHARED ## _ ## JOB ## _, SHARED, JOB> mpi_mr_ ## CALLID ## _ ## SHARED ## _ ## JOB ## _red_ ; \
                       typedef mpi_map_rect_combine<mpi_mr_ ## CALLID ## _ ## SHARED ## _ ## JOB ## _, SHARED, JOB> mpi_mr_ ## CALLID ## _ ## SHARED ## _ ## JOB ## _comb_ ; \
                       typedef mpi_parallel_call<CALLID, mpi_mr_ ## CALLID ## _ ## SHARED ## _ ## JOB ## _red_, mpi_mr_ ## CALLID ## _ ## SHARED ## _ ## JOB ## _comb_> mpi_mr_ ## CALLID ## _ ## SHARED ## _ ## JOB ## _pcall_ ; \
      } } }                                                             \
  STAN_REGISTER_MPI_DISTRIBUTED_APPLY(stan::math::internal::mpi_mr_ ## CALLID ## _ ## SHARED ## _ ## JOB ## _pcall_)


#define STAN_REGISTER_MPI_MAP_RECT_ALL(CALLID, FUNCTOR)           \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, double, double)


// redefine register macro to use MPI variant
#undef STAN_REGISTER_MAP_RECT
#define STAN_REGISTER_MAP_RECT(CALLID, FUNCTOR) \
  STAN_REGISTER_MPI_MAP_RECT_ALL(CALLID, FUNCTOR)


#endif
