#pragma once

#include <vector>
#include <map>

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/prim/arr/functor/mpi_command.hpp>
#include <stan/math/prim/arr/functor/mpi_cluster.hpp>

namespace stan {
  namespace math {

    namespace internal {

      template <int member, typename T>
      class mpi_parallel_call_cache {
        // static members to hold locally cached data
        // of placing the cache inside mpi_parallel_call is that we will
        // cache the data multiple times (once for each ReduceF
        // type).
        typedef typename const T cached_t;
        static std::map<int, cached_t> local_;

        static bool is_cached(int callsite_id) { return(local_.count(callsite_id) == 1); }
        static cached_t& lookup(int callsite_id) { return(local_.find(callsite_id)->second); }
        static void cache(int callsite_id, cached_t element) { local_.insert(std::make_pair(callsite_id, element)); }
      };
      
    }

    // utility class to store and cache static data and output size
    // per job; manages memory allocation and cluster communication
    template <typename ReduceF, typename CombineF>
    class mpi_parallel_call {
      boost::mpi::communicator world_;
      const std::size_t rank_ = world_.rank();
      const std::size_t world_size_ = world_.size();

      // local caches which hold local slices of data
      typedef internal::mpi_parallel_call_cache<1, std::vector<std::vector<double>>> t_cache_x_r;
      typedef internal::mpi_parallel_call_cache<2, std::vector<std::vector<int>>> t_cache_x_i;
      typedef internal::mpi_parallel_call_cache<3, std::vector<int>> t_cache_f_out;
      typedef internal::mpi_parallel_call_cache<4, std::vector<int>> t_cache_meta;
      
      const int callsite_id_;

      const CombineF combine_;

      typedef typename CombineF::result_type result_type;

      vector_d local_shared_params_dbl_;
      matrix_d local_job_params_dbl_;

      int num_jobs_;
      int num_local_jobs_;

    public:
      // called on root
      template <typename T_shared_param, typename T_job_param>
      mpi_parallel_call(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                        const Eigen::Matrix<T_job_param, Eigen::Dynamic, Eigen::Dynamic>& job_params,
                        const std::vector<std::vector<double>>& x_r,
                        const std::vector<std::vector<int>>& x_i,
                        int callsite_id) : callsite_id_(callsite_id), combine_(shared_params, job_params) {
        if(rank_ != 0)
          throw std::runtime_error("problem sizes can only defined on the root.");

        // make childs aware of upcoming job
        mpi_broadcast_command<stan::math::mpi_distributed_apply<mpi_parallel_call<ReduceF,CombineF>>>();

        // send callsite id
        boost::mpi::broadcast(world_, callsite_id_, 0);

        vector_d shared_params_dbl = value_of(shared_params);
        matrix_d job_params_dbl = value_of(job_params);

        setup_call(shared_params_dbl, job_params_dbl, x_r, x_i);
      }

      // called on remote sites
      mpi_parallel_call(int callsite_id) : callsite_id_(callsite_id), combine_() {
        if(rank_ == 0)
          throw std::runtime_error("problem sizes must be defined on the root.");

        setup_call(vector_d(), matrix_d(), std::vector<std::vector<double>>(), std::vector<std::vector<int>>());
      }

      static void distributed_apply() {
        // entry point when called on remote
        int callsite_id;
        boost::mpi::broadcast(world_, callsite_id, 0);

        // call constructor for the remotes
        mpi_parallel_call<ReduceF,CombineF> job_chunk(callsite_id);

        job_chunk.reduce();
      }

      // all data is cached and local parameters are also available
      result_type reduce() {

        matrix_d local_output;
        std::vector<int> f_size(num_local_jobs_);

        t_cache_x_r::cache_t& local_x_r = t_cache_x_r::lookup(callsite_id_);
        t_cache_x_i::cache_t& local_x_i = t_cache_x_i::lookup(callsite_id_);

        int offset = 0;

        for(std::size_t i=0; i != num_local_jobs_; ++i) {
          const matrix_d job_output = ReduceF::apply(local_shared_params_dbl_, local_job_params_dbl_.col(i), local_x_r[i], local_x_i[i]);
          f_size[i] = job_output.cols();
          if(i==0)
            local_output.resize(job_output.rows(), num_local_jobs_);

          if(local_output.cols() < offset + f_size[i])
            local_output.conservativeResize(Eigen::NoChange, 2*(offset + f_size[i]));

          local_output.block(0, offset, local_output.rows(), f_size[i]);

          offset += f_size[i];
        }

        // we do not need to get rid of memory which we haven't used
        // as only upto offset will be send of MPI
        
        return combine_.gather_outputs(local_output_, f_size_);
      }

    private:

      template <typename T_cache>
      void cache_data(typename T_cache::cache_t& data) {
        // distribute data only if not in cache yet
        if(T_cache::is_cached(callsite_id_))
          return;

        // we have to transpose all data
        typename T_cache::cache_t trans_data = transponse(data);

        const auto flat_data = to_array_1d(trans_data);

        // transfer it and rebuild 2D array

        // finally
        T_cache::cache(callsite_id_, local_data);
      }

      void setup_call(const vector_d& shared_params, const matrix_d& job_params,
                      const std::vector<std::vectotr<double>>& x_r,
                      const std::vector<std::vectotr<int>>& x_i) {

        if(!t_cache_meta::is_cached(callsite_id_)) {
          std::vector<int> meta(3);
          meta[0] = shared_params.size();
          meta[1] = job_params.rows();
          meta[2] = job_params.cols();
          boost::mpi::broadcast(world_, meta.data(), 3, 0);
          t_cache_meta::cache(callsite_id_, meta);
        }

        const std::vector<int>& meta = t_cache_meta::lookup(callsite_id_);

        // broadcast shared parameters
        local_shared_params_dbl_ = shared_params;
        local_shared_params_dbl_.resize(meta[0]);
        boost::mpi::broadcast(world_, local_shared_params_dbl_.data(), meta[0], 0);

        // scatter job specific parameters
        const int num_job_params = meta[1];
        
        num_jobs_ = meta[2];
        
        const std::vector<int> job_chunks = mpi_map_chunks(num_jobs_, 1);

        num_local_jobs_ = job_chunks[rank_];
        local_job_params_dbl_.resize(num_job_params, num_local_jobs_);

        const std::vector<int> job_params_chunks = mpi_map_chunks(num_jobs_, num_job_params);
        boost::mpi::scatterv(world_, job_params.data(), job_params_chunks, local_job_params_.data(), 0);

        // distribute data if not yet cached
        cache_data<t_cache_x_r>(x_r);
        cache_data<t_cache_x_i>(x_i);
      }

    };

    template <typename F, typename T_shared_param, typename T_job_param>
    struct map_rect_mpi_traits;

    template<>
    template<typename F>
    struct map_rect_traits<F,double,double> {
      typedef simple_reducer<F> reduce_t;
      typedef simple_combine<F> combine_t;
    };

    template <typename F>
    struct simple_reducer {
      static matrix_d apply(const vector_d& shared_params, const vector_d& job_params, const std::vector<double>& x_r, const std::vector<int>& x_i) {
        const vector_d out = F::apply(shared_params, job_params, x_r, x_i);
        return( out.transpose() );
      }
    };

    // combine must be instantiated to be able to hold a reference to
    // the operands for the AD cases
    template <typename F>
    struct simple_combine {
      simple_combine() {}
      simple_combine(const vector_d& shared_params, const matrix_d& job_params) {}
      typedef vector_d result_type;
      vector_d gather_outputs(...) {
        // gather stuff into a single big object which we output as VectorXd.
      }
    };

    template <typename F, typename T_shared_param, typename T_job_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1>
    map_rect_mpi(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                 const Eigen::Matrix<T_job_param, Eigen::Dynamic, Eigen::Dynamic>& job_params,
                 const std::vector<std::vector<double>>& x_r,
                 const std::vector<std::vector<int>>& x_i,
                 const int callsite_id) {
      typedef typename map_rect_traits<F, T_shared_param, T_job_param>::reduce_t ReduceF;
      typedef typename map_rect_traits<F, T_shared_param, T_job_param>::combine_t CombineF;
      
      mpi_parallel_call<ReduceF,CombineF> job_chunk(shared_params, job_params, x_r, x_i, callsite_id);

      return job_chunk.reduce();
    }
  }
}
