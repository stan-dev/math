#pragma once

#include <stan/math/prim/mat/functor/mpi_parallel_call.hpp>

#include <vector>
#include <map>
#include <type_traits>


namespace stan {
  namespace math {

    template <typename F, typename T_shared_param, typename T_job_param>
    struct map_rect_reduce;

    template <typename F>
    struct map_rect_reduce<F, double, double> {
      static std::size_t get_output_size(std::size_t num_shared_params, std::size_t num_job_specific_params) {
        return(1);
      }
      static matrix_d apply(const vector_d& shared_params, const vector_d& job_specific_params, const std::vector<double>& x_r, const std::vector<int>& x_i) {
        const F f;
        const vector_d out = f(shared_params, job_specific_params, x_r, x_i, 0);
        return( out.transpose() );
      }
    };
      
    template <typename F, typename T_shared_param, typename T_job_param>
    class map_rect_combine {
      boost::mpi::communicator world_;
      const std::size_t rank_ = world_.rank();
      const std::size_t world_size_ = world_.size();

      typedef operands_and_partials<Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>,
                                    Eigen::Matrix<T_job_param, Eigen::Dynamic, 1> > ops_partials_t;
      std::vector<ops_partials_t> ops_partials_;
      
      const std::size_t num_shared_operands_;
      const std::size_t num_job_operands_;

    public:

      typedef Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1> result_type;
      
      map_rect_combine() : ops_partials_(), num_shared_operands_(0), num_job_operands_(0) {}
      map_rect_combine(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                       const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>& job_params)
        : ops_partials_(), num_shared_operands_(shared_params.rows()), num_job_operands_(dims(job_params)[1]) {
        ops_partials_.reserve(job_params.size());
        for(std::size_t i = 0; i < job_params.size(); i++) {
          ops_partials_.push_back(ops_partials_t(shared_params, job_params[i]));
        }
      }

      result_type gather_outputs(const matrix_d& local_result, const std::vector<int>& world_f_out, const std::vector<int>& job_chunks) {

        const std::size_t num_jobs = world_f_out.size();
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

        result_type out(size_world_f_out);

        const std::size_t offset_job_params = is_constant_struct<T_shared_param>::value ? 1 : 1+num_shared_operands_ ;

        for(std::size_t i=0, ij=0; i != num_jobs; ++i) {
          for(std::size_t j=0; j != world_f_out[i]; ++j, ++ij) {
            // check if the outputs flags a failure
            if(unlikely(world_result(0,ij) == std::numeric_limits<double>::max())) {
              throw std::runtime_error("MPI error.");
            }

            if (!is_constant_struct<T_shared_param>::value) {
              ops_partials_[i].edge1_.partials_ = world_result.block(1,ij,num_shared_operands_,1);
            }
              
            if (!is_constant_struct<T_job_param>::value) {
              ops_partials_[i].edge2_.partials_ = world_result.block(offset_job_params,ij,num_job_operands_,1);
            }
            
            out(ij) = ops_partials_[i].build(world_result(0,ij));
          }
        }

        return(out);
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
      typedef map_rect_combine<F, T_shared_param, T_job_param> CombineF;
      
      mpi_parallel_call<call_id,ReduceF,CombineF> job_chunk(shared_params, job_params, x_r, x_i);

      return job_chunk.reduce();
    }
  }
}
