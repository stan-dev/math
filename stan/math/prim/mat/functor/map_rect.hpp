#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_MAP_RECT_HPP

#define STAN_REGISTER_MAP_RECT(CALLID, FUNCTOR) 

#include <stan/math/prim/scal/meta/operands_and_partials.hpp>

#ifdef STAN_HAS_MPI
#include <stan/math/prim/mat/functor/map_rect_mpi.hpp>
#endif

namespace stan {
  namespace math {

    template <typename F, typename T_shared_param, typename T_job_param>
    class map_rect_reduce {
    };

    template <typename F>
    class map_rect_reduce<F, double, double> {
    public:
      matrix_d operator()(const vector_d& shared_params, const vector_d& job_specific_params,
                          const std::vector<double>& x_r, const std::vector<int>& x_i) const {
        const F f;
        const vector_d out = f(shared_params, job_specific_params, x_r, x_i, 0);
        return( out.transpose() );
      }
    };
      
    template <typename F, typename T_shared_param, typename T_job_param>
    class map_rect_combine {

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

      result_type operator()(const matrix_d& world_result, const std::vector<int>& world_f_out) {

        const std::size_t num_jobs = world_f_out.size();
        const std::size_t offset_job_params = is_constant_struct<T_shared_param>::value ? 1 : 1+num_shared_operands_ ;
        const std::size_t size_world_f_out = sum(world_f_out);

        result_type out(size_world_f_out);

        for(std::size_t i=0, ij=0; i != num_jobs; ++i) {
          for(std::size_t j=0; j != world_f_out[i]; ++j, ++ij) {
            // check if the outputs flags a failure
            if(unlikely(world_result(0,ij) == std::numeric_limits<double>::max()))
              throw std::runtime_error("Error.");

            if (!is_constant_struct<T_shared_param>::value)
              ops_partials_[i].edge1_.partials_ = world_result.block(1,ij,num_shared_operands_,1);
              
            if (!is_constant_struct<T_job_param>::value)
              ops_partials_[i].edge2_.partials_ = world_result.block(offset_job_params,ij,num_job_operands_,1);
            
            out(ij) = ops_partials_[i].build(world_result(0,ij));
          }
        }

        return(out);
      }
    };
    
    template <int call_id, typename F, typename T_shared_param, typename T_job_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1>
    map_rect_serial(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
                    const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1> >& job_params,
                    const std::vector<std::vector<double> >& x_r,
                    const std::vector<std::vector<int> >& x_i,
                    std::ostream* msgs = 0) {
      check_nonzero_size("map_rect_serial", "job parameters", job_params);
      check_matching_sizes("map_rect_serial", "job parameters", job_params, "continuous data", x_r);
      check_matching_sizes("map_rect_serial", "job parameters", job_params, "integer data", x_i);

      typedef map_rect_reduce< F, T_shared_param, T_job_param> ReduceF;
      typedef map_rect_combine<F, T_shared_param, T_job_param> CombineF;

      const int num_jobs = job_params.size();
      const vector_d shared_params_dbl = value_of(shared_params);
      
      matrix_d world_output(0, 0);
      std::vector<int> world_f_out(num_jobs, -1);

      int offset = 0;
      for(std::size_t i=0; i < num_jobs; ++i) {
        const matrix_d job_output = ReduceF()(shared_params_dbl, value_of(job_params[i]), x_r[i], x_i[i]);
        world_f_out[i] = job_output.cols();

        if(i == 0) world_output.resize(job_output.rows(), num_jobs * world_f_out[i]);

        if(world_output.cols() < offset + world_f_out[i])
          world_output.conservativeResize(Eigen::NoChange, 2*(offset + world_f_out[i]));

        world_output.block(0, offset, world_output.rows(), world_f_out[i]) = job_output;

        offset += world_f_out[i];
      }

      CombineF combine(shared_params, job_params);
      return combine(world_output, world_f_out);

      /* obsolete version which builds an AD graph       
      typedef typename stan::return_type<T_shared_param, T_job_param>::type result_type;
      Eigen::Matrix<result_type, Eigen::Dynamic, 1> out;
      const std::size_t num_jobs = job_params.size();
      int out_size = 0;
      const F f;
      
      for(std::size_t i = 0; i != num_jobs; ++i) {
        const Eigen::Matrix<result_type, Eigen::Dynamic, 1> fx = f(shared_params, job_params[i], x_r[i], x_i[i], msgs);
        const int fx_size = fx.rows();
        out_size += fx_size;
        if(i == 0) out.resize(num_jobs * fx_size);
        if(out.rows() < out_size)
          out.conservativeResize(2*out_size);
        out.segment(out_size - fx_size, fx_size) = fx;
      }
      out.conservativeResize(out_size);
      return(out);
      */
    }

    template <int call_id, typename F, typename T_shared_param, typename T_job_param>
    Eigen::Matrix<typename stan::return_type<T_shared_param, T_job_param>::type, Eigen::Dynamic, 1>
    map_rect(const Eigen::Matrix<T_shared_param, Eigen::Dynamic, 1>& shared_params,
             const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1> >& job_params,
             const std::vector<std::vector<double> >& x_r,
             const std::vector<std::vector<int> >& x_i,
             std::ostream* msgs = 0) {
#ifdef STAN_HAS_MPI
      return(map_rect_mpi<call_id,F,T_shared_param,T_job_param>(shared_params, job_params, x_r, x_i, msgs));
#else
      return(map_rect_serial<call_id,F,T_shared_param,T_job_param>(shared_params, job_params, x_r, x_i, msgs));
#endif
    }

  }
}

#endif
