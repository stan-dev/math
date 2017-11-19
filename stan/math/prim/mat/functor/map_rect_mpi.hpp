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
                        const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>& job_params,
                        const std::vector<std::vector<double>>& x_r,
                        const std::vector<std::vector<int>>& x_i,
                        int callsite_id) : callsite_id_(callsite_id), combine_(shared_params, job_params) {
        if(rank_ != 0)
          throw std::runtime_error("problem sizes can only defined on the root.");

        // make childs aware of upcoming job
        mpi_broadcast_command<stan::math::mpi_distributed_apply<mpi_parallel_call<ReduceF,CombineF>>>();

        // send callsite id
        boost::mpi::broadcast(world_, callsite_id_, 0);

        const std::size_t num_jobs = job_params.size();
        const std::size_t num_job_params = job_params[0].size();
        
        vector_d shared_params_dbl = value_of(shared_params);
        matrix_d job_params_dbl(num_job_params, num_jobs);

        for(std::size_t j=0; j != num_jobs; ++j)
          job_params_dbl.col(j) = value_of(job_params[j]);

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

        const std::vector<int> job_chunks = mpi_map_chunks(num_jobs_, 1);

        // id of first job out of all
        int start_job = 0;
        for(std::size_t n=0; n != rank_; ++n)
          start_job += job_chunks[n];

        matrix_d local_output(1, num_local_jobs_);
        std::vector<int> local_f_out(num_local_jobs_);

        t_cache_x_r::cache_t& local_x_r = t_cache_x_r::lookup(callsite_id_);
        t_cache_x_i::cache_t& local_x_i = t_cache_x_i::lookup(callsite_id_);

        // check if we know already output sizes
        if(t_cache_f_out::is_cached(callsite_id_)) {
          t_cache_f_out::cache_t& f_out = t_cache_f_out::lookup(callsite_id_);
          int num_outputs = 0;
          for(std::size_t j=start_job; j != start_job + num_local_jobs; ++j)
            num_outputs += f_out[j];
          local_output.resize(1, num_outputs);
        }

        int offset = 0;

        // TODO: add exception handling
        for(std::size_t i=0; i != num_local_jobs_; ++i) {
          const matrix_d job_output = ReduceF::apply(local_shared_params_dbl_, local_job_params_dbl_.col(i), local_x_r[i], local_x_i[i]);
          local_f_out[i] = job_output.cols();
          
          if(i==0)
            local_output.resize(job_output.rows(), Eigen::NoChange);

          if(local_output.cols() < offset + local_f_out[i])
            local_output.conservativeResize(Eigen::NoChange, 2*(offset + local_f_out[i]));

          local_output.block(0, offset, local_output.rows(), local_f_out[i]) = job_output;

          offset += local_f_out[i];
        }

        // during first execution we distribute the output sizes from
        // local jobs to the root
        if(!t_cache_f_out::is_cached(callsite_id_)) {
          std::vector<int> world_f_out(num_jobs_, 0);
          boost::mpi::gatherv(world_, local_f_out.data(), num_local_jobs_, world_f_out.data(), job_chunks, 0);
          // on the root we now have all sizes from all childs. Copy
          // over the local sizes to the world vector on each local
          // node in order to cache this information locally
          std::copy(local_f_out.begin(), local_f_out.end(), world_f_out.begin() + start_jobs);
          t_cache_f_out::cache(callsite_id_, world_f_out);
        }

        t_cache_f_out::cache_t& world_f_out = t_cache_f_out::lookup(callsite_id_);

        // check that cached sizes are the same as just collected from
        // this evaluation

        return combine_.gather_outputs(local_output_, world_f_out);
      }

    private:

      template <typename T_cache>
      void cache_data(typename T_cache::cache_t& data) {
        // distribute data only if not in cache yet
        if(T_cache::is_cached(callsite_id_))
          return;

        const std::size_t size_data = data[0].size();
        const std::vector<int> data_chunks = mpi_map_chunks(num_jobs_, size_data);
        const std::vector<int> job_chunks = mpi_map_chunks(num_jobs_, 1);

        const auto flat_data = to_array_1d(data);

        // transfer it ...
        decltype(flat_data) local_flat_data(data_chunks[rank_]);
        boost::mpi::scatterv(world_, flat_data.data(), data_chunks, local_flat_data.data(), 0);

        // ... and rebuild 2D array
        T_cache::cache_t local_data;
        auto local_iter = local_flat_data.begin();
        for(std::size_t i=0; i != job_chunks[rank_]; ++i) {
          typename T_cache::cache_t::value_type const data_elem(local_iter, local_iter + size_data);
          local_data.push_back(data_elem);
          local_iter += size_data;
        }

        // finally we cache it locally
        T_cache::cache(callsite_id_, local_data);
      }

      void setup_call(const vector_d& shared_params, const matrix_d& job_params,
                      const std::vector<std::vector<double>>& x_r,
                      const std::vector<std::vector<int>>& x_i) {

        if(!t_cache_meta::is_cached(callsite_id_)) {
          std::vector<int> meta(3);
          meta[0] = shared_params.size();
          meta[1] = job_params.rows();
          meta[2] = job_params.cols();
          boost::mpi::broadcast(world_, meta.data(), 3, 0);
          t_cache_meta::cache(callsite_id_, meta);
        }

        const std::vector<int>& meta = t_cache_meta::lookup(callsite_id_);
        const int num_shared_params = meta[0];
        const int num_job_params = meta[1];
        num_jobs_ = meta[2];

        // broadcast shared parameters
        local_shared_params_dbl_ = shared_params;
        local_shared_params_dbl_.resize(num_shared_params);
        boost::mpi::broadcast(world_, local_shared_params_dbl_.data(), num_shared_params, 0);

        // scatter job specific parameters
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
    struct map_rect_reduce;

    template <>
    template <typename F>
    struct map_rect_reduce<F, double, double> {
      static matrix_d apply(const vector_d& shared_params, const vector_d& job_specific_params, const std::vector<double>& x_r, const std::vector<int>& x_i) {
        const vector_d out = F::apply(shared_params, job_specific_params, x_r, x_i);
        return( out.transpose() );
      }
    };
      
    template <typename F, typename T_shared_param, typename T_job_param>
    class map_rect_combine {
      boost::mpi::communicator world_;
      const std::size_t rank_ = world_.rank();
      const std::size_t world_size_ = world_.size();

      const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>* shared_params_operands_;
      const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>* job_params_operands_;
      
    public:

      typedef Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1> result_type;
      
      map_rect_combine() {}
      map_rect_combine(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                       const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>& job_params)
        : shared_params_operands_(&shared_params), job_params_operands_(&job_params) {}

      result_type gather_outputs(const matrix_d& local_result, const std::vector<int>& world_f_out) {

        const std::size_t num_jobs = world_f_out.size();
        const std::size_t num_output_per_job = local_result.rows();
        const std::size_t size_world_f_out = sum(world_f_out);

        matrix_d world_result(num_output_per_job, size_world_f_out);

        const std::vector<int> job_chunks = mpi_map_chunks(num_jobs, 1);
        std::vector<int> chunks_result(world_size_,0);
        for(std::size_t i=0, ij=0; i != world_size_; ++i)
          for(std::size_t j=0; j != job_chunks[i]; ++j, ++ij)
            chunks_result[i] += world_f_out_[ij] * num_output_per_job;

        // collect results
        //std::cout << "gathering actual outputs..." << std::endl;
        boost::mpi::gatherv(world_, local_result.data(), local_result.size(), world_result, chunks_result, 0);

        // only add result to AD tree on root
        if(rank_ != 0)
          return(result_type());

        result_type out(size_world_f_out);

        const std::size_t num_shared_operands = *shared_params_operands_.size();
        const std::size_t num_job_operands = *job_params_operands_.size();

        const std::size_t offset_job_params = is_constant_struct<T_shared_param>::value ? 1 : 1+num_shared_operands ;

        for(std::size_t i=0, ij=0; i != num_jobs; ++i) {
          for(std::size_t j=0; j != world_f_out[i]; ++j, ++ij) {
            operands_and_partials<T_shared_param, T_job_param>
              ops_partials(*shared_params_operands_, *job_params_operands_[i]);

            if (!is_constant_struct<T_shared_param>::value) {
              for(std::size_t k=0; k != num_shared_operands; ++k)
                ops_partials.edge1_.partials[k] = world_result(1+k,ij);
            }
              
            if (!is_constant_struct<T_job_param>::value) {
              for(std::size_t k=0; k != num_job_operands; ++k)
                ops_partials.edge2_.partials[k] = world_result(offset_job_params+k,ij);
            }

            out(ij) = ops_partials.build(world_result(1,ij));
          }
        }

        return(out);
      }
    };

    template <typename F, typename T_shared_param, typename T_job_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1>
    map_rect_mpi(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                 const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>& job_params,
                 const std::vector<std::vector<double>>& x_r,
                 const std::vector<std::vector<int>>& x_i,
                 const int callsite_id) {
      typedef typename map_rect_reduce<F, T_shared_param, T_job_param> ReduceF;
      typedef typename map_rect_combine<F, T_shared_param, T_job_param> CombineF;
      
      mpi_parallel_call<ReduceF,CombineF> job_chunk(shared_params, job_params, x_r, x_i, callsite_id);

      return job_chunk.reduce();
    }
  }
}
