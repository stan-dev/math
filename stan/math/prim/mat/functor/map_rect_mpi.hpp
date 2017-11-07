#pragma once

#include <vector>
#include <map>

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/prim/arr/functor/mpi_command.hpp>
#include <stan/math/prim/arr/functor/mpi_cluster.hpp>

namespace stan {
  namespace math {

    namespace internal {
      class mpi_parallel_call_cache {
        // static members to hold locally cached data
        // of placing the cache inside mpi_parallel_call is that we will
        // cache the data multiple times (once for each ReduceF
        // type). We should probably should move this into a non-templated
        // class.
        static std::map<int, Eigen::MatrixXd> cache_x_r_;
        static std::map<int, std::vector<std::vector<int>>> cache_x_i_;
        static std::map<int, std::vector<int>> cache_f_size_;
        
      };
    }

    // utility class to store and cache static data and output size
    // per job; manages memory allocation and cluster communication
    template <typename T_shared_param, typename T_param, typename ReduceF, typename CombineF>
    class mpi_parallel_call {
      boost::mpi::communicator world_;
      const std::size_t rank_ = world_.rank();
      const std::size_t world_size_ = world_.size();

      // these point for a given instance to the data to be used
      std::map<int, Eigen::MatrixXd>::const_iterator local_x_r_;
      std::map<int, std::vector<std::vector<int>>>::const_iterator local_x_i_;
      std::map<int, std::vector<int>>::const_iterator local_f_size_;

      Eigen::MatrixXd local_output_;

      // initialized on root... note: we need references to the
      // operands only on the root. We could consider moving this into
      // an instance of CombineF since we only need these there
      const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>* eta_;
      const Eigen::Matrix<T_param, Eigen::Dynamic, Eigen::Dynamic>* theta_;
        
      Eigen::VectorXd eta_dbl_;
      Eigen::MatrixXd local_theta_dbl_;

      const int callsite_id_;

      typedef typename CombineF::result_type result_type;

    public:
      // called on root
      mpi_parallel_call(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& eta,
                        const Eigen::Matrix<T_param, Eigen::Dynamic, Eigen::Dynamic>& theta,
                        const Eigen::MatrixXd& x_r,
                        const std::vector<std::vector<int> >& x_i,
                        const int callsite_id) {
        if(rank_ != 0)
          throw std::runtime_error("problem sizes can only defined on the root.");

        // make childs aware of upcoming job
        mpi_broadcast_command<stan::math::mpi_distributed_apply<mpi_parallel_call<T_shared_param,T_param,ReduceF,CombineF>>>();

        // broadcast meta info like callsite

        // check of callsite has been seen already
        // if-first-call        
        // setup_local_data(x_r, x_i);
        // else
        // grab data from caches

        distribute_param();

        // if-first-call evaluate function with double args and record
        // how many outputs we get for each job. In case an exception
        // is thrown, we have to clear all the data cached so far!
      }

      // called on remote sites
      mpi_parallel_call() {
        if(rank_ == 0)
          throw std::runtime_error("problem sizes must be defined on the root.");

        // mirror constructor actions done on root
      }

      // mpi communication using caching... we return a const_iterator
      // to the slice of data which we process locally; is there a
      // better type to return? I only need a const reference to be
      // returned.
      std::map<int, std::vector<std::vector<int>>>::const_iterator scatter(const std::vector<int>& in_values);
      std::map<int, Eigen::MatrixXd>::const_iterator scatter(const Eigen::MatrixXd& in_values);

      static void distributed_apply() {
        // entry point when called on remote

        // call default constructor
        mpi_parallel_call<T_shared_param,T_param,ReduceF,CombineF> job_chunk;

        job_chunk.reduce();
      }

      // all data is cached and local parameters are also available
      result_type reduce() {
        // get output size for each element which includes all the
        // calculated gradients, etc.
        const std::size_t num_element_outputs = ReduceF::get_element_outputs(eta.rows(), theta.rows());
        // allocate final storage of all local outputs
        local_output_.resize(num_element_outputs, sum(local_f_size->second));

        const std::size_t num_local_jobs = local_theta_.cols();
        std::size_t offset = 0;

        // evaluate for each job ReduceF which stores results directly
        // in final location
        for(std::size_t i = 0; i != num_local_jobs; ++i) {
          const std::size_t f_size = local_f_size_[i];
          ReduceF::apply(eta_, local_theta_.col(i), local_x_r_.col(i), local_x_i_[i], local_output_.block(0, offset, num_element_outputs, f_size));
          offset += f_size;
        }

        return CombineF::gather_outputs(eta_, theta_, f_size_, local_output_);
      }

    };

    template <typename F, typename T_shared_param, typename T_param>
    struct map_rect_mpi_traits;

    template<>
    template<typename F>
    struct map_rect_traits<F,double,double> {
      typedef simple_reducer<F> reduce_t;
      typedef simple_combine<F> combine_t;
    };

    template <typename F>
    struct simple_reducer {
      static std::size_t get_elements_output(std::size_t size_eta, std::size_t size_theta) { return(1) };
      static void apply(const Eigen::VectorXd& eta, const Eigen::VectorXd& theta, const Eigen::VectorXd& x_r, const std::vector<int>& x_i, Eigen::MatrixXd& out) {
        out = F::apply(eta, theta, x_r, x_i);
      }
    };

    template <typename F>
    struct simple_combine {
      typedef Eigen::VectorXd result_type;
      static Eigen::VectorXd gather_outputs(...) {
        // gather stuff into a single big object which we output as VectorXd.
      }
    };

    template <typename F, typename T_shared_param, typename T_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_param>::type, Eigen::Dynamic, 1>
    map_rect_mpi(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& eta,
                 const Eigen::Matrix<T_param, Eigen::Dynamic, Eigen::Dynamic>& theta,
                 const Eigen::MatrixXd& x_r,
                 const std::vector<std::vector<int> >& x_i,
                 const int callsite_id) {
      typedef typename map_rect_traits<T_shared_param, T_param, F>::reduce_t ReduceF;
      typedef typename map_rect_traits<T_shared_param, T_param, F>::combine_t CombineF;
      
      mpi_parallel_call<T_shared_param,T_param,ReduceF,CombineF> job_chunk(eta, theta, x_r, x_i, callsite_id);

      return job_chunk.reduce();
    }
  }
}
