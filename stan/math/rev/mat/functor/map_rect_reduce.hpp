#pragma once

namespace stan {
  namespace math {

    template <typename F>
    struct map_rect_reduce<F, var, var> {
      static std::size_t get_output_size(std::size_t num_shared_params, std::size_t num_job_specific_params) {
        return(1+num_shared_params+num_job_specific_params);
      }
      static matrix_d apply(const vector_d& shared_params, const vector_d& job_specific_params, const std::vector<double>& x_r, const std::vector<int>& x_i) {
        const size_type num_shared_params = shared_params.rows();
        const size_type num_job_specific_params = job_specific_params.rows();
        const size_type num_params = num_shared_params  + num_job_specific_params;
        matrix_d out(1+num_params,0);

        try {
          start_nested();
          vector_v shared_params_v(num_shared_params);
          vector_v job_specific_params_v(num_job_specific_params);

          for(std::size_t i=0; i < num_shared_params; ++i)
            shared_params_v(i) = shared_params(i);
          for(std::size_t i=0; i < num_job_specific_params; ++i)
            job_specific_params_v(i) = job_specific_params(i);

          std::vector<var> z_vars(num_params);
          std::vector<double> z_grad(num_params);

          for(size_type i = 0; i < num_shared_params; ++i)
            z_vars[i] = shared_params_v(i);
          for(size_type i = 0; i < num_job_specific_params; ++i)
            z_vars[num_shared_params + i] = job_specific_params_v(i);

          const F f;
          vector_v fx_v = f(shared_params_v, job_specific_params_v, x_r, x_i, 0);

          const size_t size_f = fx_v.rows();

          out.resize(Eigen::NoChange, size_f);

          for(size_type i = 0; i < size_f; ++i) {
            out(0,i) = fx_v(i).val();            
            set_zero_all_adjoints_nested();
            fx_v(i).grad(z_vars, z_grad);
            out.block(1,i,num_params,1) = Eigen::Map<vector_d>(&z_grad[0], num_params);
          }
          recover_memory_nested();
        } catch(const std::exception& e) {
          recover_memory_nested();
          throw;
        }
        return( out );
      }
    };
    
    template <typename F>
    class map_rect_combine<F, var, var> {
      boost::mpi::communicator world_;
      const std::size_t rank_ = world_.rank();
      const std::size_t world_size_ = world_.size();

      typedef var T_shared_param;
      typedef var T_job_param;

      const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>* shared_params_operands_;
      const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>* job_params_operands_;
      
    public:

      typedef Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1> result_type;
      
      map_rect_combine() {}
      map_rect_combine(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                       const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>& job_params)
        : shared_params_operands_(&shared_params), job_params_operands_(&job_params) {}

      result_type gather_outputs(const matrix_d& local_result, const std::vector<int>& world_f_out, const std::vector<int>& job_chunks) const {

        const std::size_t num_jobs = world_f_out.size();
        const std::size_t num_output_size_per_job = local_result.rows();
        const std::size_t size_world_f_out = sum(world_f_out);

        // allocate memory on AD stack for final result. Note that the
        // gradients and the function results will land there
        double *world_result = 0;

        if(rank_ == 0)
          world_result = ChainableStack::memalloc_.alloc_array<double>( num_output_size_per_job * size_world_f_out );
            
        std::vector<int> chunks_result(world_size_,0);
        for(std::size_t i=0, ij=0; i != world_size_; ++i)
          for(std::size_t j=0; j != job_chunks[i]; ++j, ++ij)
            chunks_result[i] += world_f_out[ij] * num_output_size_per_job;

        // collect results
        boost::mpi::gatherv(world_, local_result.data(), chunks_result[rank_], world_result, chunks_result, 0);

        // only add result to AD tree on root
        if(rank_ != 0)
          return(result_type());

        result_type out(size_world_f_out);

        const std::size_t num_shared_operands = shared_params_operands_->size();
        const std::vector<int> dims_job_operands = dims(*job_params_operands_);
        const std::size_t num_job_operands = dims_job_operands[1]; //(*job_params_operands_)[0].size();
        const std::size_t num_operands = num_shared_operands + num_job_operands;

        vari** varis
          = ChainableStack::memalloc_.alloc_array<vari*>(num_jobs*(num_shared_operands+num_job_operands));

        for(std::size_t i=0, ik=0; i != num_jobs; i++) {
          // link the operands...
          for(std::size_t j=0; j != num_shared_operands; j++)
            varis[i * (num_operands) + j] = (*shared_params_operands_)(j).vi_;
          for(std::size_t j=0; j != num_job_operands; j++)
            varis[i * (num_operands) + num_shared_operands + j] = (*job_params_operands_)[i](j).vi_;
          
          // ...with partials of outputs
          for(std::size_t k=0; k != world_f_out[i]; k++, ik++) {
            const double val = *(world_result + (num_operands+1) * ik);
            if(unlikely(val == std::numeric_limits<double>::max())) {
              std::cout << "THROWING ON RANK " << rank_ << std::endl;
              throw std::runtime_error("MPI error.");
            }
            out(ik) = var(new precomputed_gradients_vari(val, num_operands, varis + i * (num_operands), world_result + (num_operands+1) * ik + 1));
          }
        }
      
        return(out);
      }
    };

  }
}

#undef STAN_REGISTER_MPI_MAP_RECT_ALL

#define STAN_REGISTER_MPI_MAP_RECT_ALL(CALLID, FUNCTOR)         \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, double, double) \
  STAN_REGISTER_MPI_MAP_RECT(CALLID, FUNCTOR, var, var)
