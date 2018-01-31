#ifndef STAN_MATH_REV_MAT_FUNCTOR_IDAS_RESIDUAL_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_IDAS_RESIDUAL_HPP

#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
// #include <stan/math/rev/mat/functor/cvodes_ode_data.hpp>
// #include <stan/math/rev/arr/fun/decouple_ode_states.hpp>
#include <idas/idas.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <ostream>
#include <vector>

#define CHECK_IDAS_CALL(call) idas_check(call, #call)

inline void idas_check(int flag, const char* func) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << func << " failed with error flag " << flag;
    throw std::runtime_error(ss.str());
  }
}

namespace stan {
namespace math {

template<typename F, typename TYY, typename TYP, typename TPAR>
class idas_system {
 protected:
  const F& f_;
  const std::vector<TYY>& yy_;
  const std::vector<TYP>& yp_;
  const std::vector<TPAR>& theta_;
  const std::vector<double>& x_r_;
  const std::vector<int>& x_i_;
  const size_t N_;
  const size_t M_;
  const size_t n_sens_;
  N_Vector nv_yy_;
  N_Vector nv_yp_;
  void* mem_;
  SUNMatrix A_;
  SUNLinearSolver LS_;
  std::ostream* msgs_;

 public:
  static constexpr bool is_var_yy0 = stan::is_var<TYY>::value;
  static constexpr bool is_var_yp0 = stan::is_var<TYP>::value;
  static constexpr bool is_var_par = stan::is_var<TPAR>::value;
  static constexpr bool need_sens = is_var_yy0 && is_var_yp0 && is_var_par;

  using return_type =  std::vector<
    std::vector<typename stan::return_type<TYY, TYP, TPAR>::type> >;

  idas_system(const F& f,
                        const std::vector<TYY>& yy0,
                        const std::vector<TYP>& yp0,
                        const std::vector<TPAR>& theta,
                        const std::vector<double>& x_r,
                        const std::vector<int>& x_i,
                        std::ostream* msgs) :
    f_(f), yy_(yy0), yp_(yp0), theta_(theta), x_r_(x_r), x_i_(x_i), 
    N_(yy0.size()),
    M_(theta.size()),
    n_sens_((is_var_yy0? 0 : N_) + (is_var_yp0? 0 : N_) + (is_var_par? 0 : M_)),
    mem_(IDACreate()),
    msgs_(msgs)
  {
    nv_yy_ = N_VMake_Serial(N_, value_of(yy_).data());
    nv_yp_ = N_VMake_Serial(N_, value_of(yp_).data());

    A_ = SUNDenseMatrix(N_, N_);
    LS_ = SUNDenseLinearSolver(nv_yy_, A_);

    if (mem_ == NULL)
      throw std::runtime_error("IDACreate failed to allocate memory");    
  }

  ~idas_system(){
    SUNLinSolFree(LS_);
    SUNMatDestroy(A_);
    N_VDestroy_Serial(nv_yy_);
    N_VDestroy_Serial(nv_yp_);
    IDAFree(&mem_);
  }

  N_Vector& nv_yy() {return nv_yy_;}
  N_Vector& nv_yp() {return nv_yp_;}
  SUNMatrix& jacobi() {return A_;}
  SUNLinearSolver& linsol() {return LS_;}

  std::vector<double> yy() {return value_of(yy_);}
  std::vector<double> yp() {return value_of(yp_);}

  const size_t& n() {return N_;}

  const size_t& n_sens() {return n_sens_;}

  const size_t n_sys() {return N_ * (n_sens_ + 1);}

  void* mem() {return mem_;}

  void set_consistent_ic(double t1){
    if(is_var_yp0 && !is_var_yy0) {
      CHECK_IDAS_CALL(IDACalcIC(mem_, IDA_YA_YDP_INIT, t1));
      CHECK_IDAS_CALL(IDAGetConsistentIC(mem_, nv_yy_, nv_yp_));
    } else {
      CHECK_IDAS_CALL(IDACalcIC(mem_, IDA_Y_INIT, t1));
      CHECK_IDAS_CALL(IDAGetConsistentIC(mem_, nv_yy_, NULL));
    }
  }

  IDAResFn residual() {         // return a non-capture lambda
    return [](double t, N_Vector yy, N_Vector yp,
              N_Vector rr, void *user_data) -> int {
      using DAE = idas_system<F, TYY, TYP, TPAR>;
      DAE* dae = static_cast<DAE*>(user_data);

      auto N = NV_LENGTH_S(yy);
      auto yy_val = N_VGetArrayPointer(yy);
      std::vector<double> yy_vec(yy_val, yy_val + N);
      auto yp_val = N_VGetArrayPointer(yp); 
      std::vector<double> yp_vec(yp_val, yp_val + N);
      auto res = dae -> f_(t, yy_vec, yp_vec,
                           dae -> theta_,
                           dae -> x_r_,
                           dae -> x_i_,
                           dae -> msgs_);
      NV_DATA_S(rr) = res.data();
      return 0;
    };    
  }

};

template<typename F, typename TYY, typename TYP, typename TPAR>
class idas_forward_system: public idas_system<F, TYY, TYP, TPAR> {
  N_Vector* nv_yys_;
  N_Vector* nv_yps_;

public:
  idas_forward_system(const F& f,
                      const std::vector<TYY>& yy0,
                      const std::vector<TYP>& yp0,
                      const std::vector<TPAR>& theta,
                      const std::vector<double>& x_r,
                      const std::vector<int>& x_i,
                      std::ostream* msgs) :
    idas_system<F, TYY, TYP, TPAR>(f, yy0, yp0, theta, x_r, x_i, msgs) 
  {
    if (this -> need_sens) {
      nv_yys_ = N_VCloneVectorArray(this->n_sens_, this->nv_yy_);
      nv_yps_ = N_VCloneVectorArray(this->n_sens_, this->nv_yp_);
      for (size_t is=0;is<this->n_sens_;is++) {
        N_VConst(RCONST(0.0), nv_yys_[is]);
        N_VConst(RCONST(0.0), nv_yps_[is]);
      }
    }
  }

  ~idas_forward_system() {
    if (this -> need_sens) {
      N_VDestroyVectorArray_Serial(this->nv_yys_, this->n_sens_);
      N_VDestroyVectorArray_Serial(this->nv_yps_, this->n_sens_);
    }
  }

  N_Vector* nv_yys() {return nv_yys_;}
  N_Vector* nv_yps() {return nv_yps_;}
 
  void set_consistent_sens_ic() {
    if(this->is_var_yp0 && !this->is_var_yy0) {
      CHECK_IDAS_CALL(IDAGetSensConsistentIC(this->mem_, nv_yys_, nv_yps_));
    } else {
      CHECK_IDAS_CALL(IDAGetSensConsistentIC(this->mem_, nv_yys_, NULL));
    }
  }

  void cast_to_user_data() {        // inject dae info
    void* user_data = static_cast<void*>(this);
    CHECK_IDAS_CALL(IDASetUserData(this->mem_, user_data));
  }

  IDASensResFn sensitivity_residual(){
    return [](int n_sens, double t,
              N_Vector yy, N_Vector yp, N_Vector res,
              N_Vector* yys, N_Vector* yps,
              N_Vector* ress, void *user_data,
              N_Vector temp1, N_Vector temp2, N_Vector temp3) {
      using Eigen::Matrix;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      using Eigen::Dynamic;

      using DAE = idas_forward_system<F, TYY, TYP, TPAR>;

      DAE* dae = static_cast<DAE*>(user_data);

      const size_t& N = dae->N_;
      const size_t& M = dae->M_;

      Eigen::Map<VectorXd> vec_yy(N_VGetArrayPointer(yy), N);
      Eigen::Map<VectorXd> vec_yp(N_VGetArrayPointer(yp), N);
      Eigen::Map<VectorXd> vec_par(value_of(dae->theta_).data(), N);
      std::vector<double> vyy(vec_yy.data(), vec_yy.data() + N);
      std::vector<double> vyp(vec_yp.data(), vec_yp.data() + N);
        
      auto nv_to_mat = [](N_Vector* nv, size_t m, size_t n) {
        matrix_d mat;
        for(size_t j=0; j<m; ++j) {
          for(size_t i=0; i<n; ++i) {
            auto nvp = N_VGetArrayPointer(nv[i]);
            mat(j, i) = nvp[j];
          }
        }
        return mat;
      };
    
      auto yys_mat = nv_to_mat(yys, N, n_sens);
      auto yps_mat = nv_to_mat(yps, N, n_sens);

      try {
        stan::math::start_nested();

        // vecvar yy_var(yy_val, yy_val+N);
        // vecvar yp_var(yp_val, yp_val+N);

        MatrixXd J, r;
        VectorXd f_val;

        auto fyy = [&t, &vyp, &N, &dae](const matrix_v& x) {
          std::vector<var> yy(x.data(), x.data() + N);
          auto eval = dae -> f_(t, yy, vyp,
                           dae->theta_,
                           dae->x_r_,
                           dae->x_i_,
                           dae->msgs_);
          Eigen::Map<vector_v> res(eval.data(), N);
          return res;
        };
        stan::math::jacobian(fyy, vec_yy, f_val, J);
        r = J*yys_mat;

        auto fyp = [&t, &vyy, &N, &dae](const matrix_v& x) {
          std::vector<var> yp(x.data(), x.data() + N);
          auto eval = dae -> f_(t, vyy, yp,
                           dae->theta_,
                           dae->x_r_,
                           dae->x_i_,
                           dae->msgs_);
          Eigen::Map<vector_v> res(eval.data(), N);
          return res;
        };
        stan::math::jacobian(fyp, vec_yp, f_val, J);
        r += J*yps_mat;

        if (dae -> is_var_par) {
          auto fpar = [&t, &vyy, &vyp, &N, &M, &dae](const matrix_v& x) {
            std::vector<var> par(x.data(), x.data() + M);
            auto eval = dae -> f_(t, vyy, vyp,
                                  par,
                                  dae->x_r_,
                                  dae->x_i_,
                                  dae->msgs_);
            Eigen::Map<vector_v> res(eval.data(), N);
            return res;
          };
          stan::math::jacobian(fpar, vec_par, f_val, J);
          r += J;
        }

        for(size_t j=0; j<N; ++j) {
          for(size_t i=0; i<n_sens; ++i) {
            auto nvp = N_VGetArrayPointer(ress[i]);
            nvp[j] = r(j, i);
          }
        }
      
      } catch (const std::exception& e) {
        stan::math::recover_memory_nested();
        throw;
      }
      stan::math::recover_memory_nested();

      return 0;
    };
  }

  // inline std::vector<
  //   std::vector<typename stan::return_type<TYY, TYP, TPAR>::type> >
  // solution(const std::vector<std::vector<double> >& res_yy,
  //          const std::vector<MatrixXd>& res_ys) {

  //   constexpr bool has_par = stan::is_var<TPAR>::type;

  //   std::vector<typename stan::return_type<TYY, TYP, TPAR>::type> vars;
  //   vars

  //   precomputed_gradients(y[i][j], vars, temp_gradients);

  // }

  // clang++ -Wall -I . -isystem lib/eigen_3.3.3 -isystem lib/boost_1.65.1 -isystem lib/idas-2.1.0/include -std=c++1y -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION -Wno-unused-function -Wno-uninitialized -stdlib=libc++ -Wno-unknown-warning-option -Wno-tautological-compare -Wsign-compare -DNO_FPRINTF_OUTPUT -pipe lib/idas-2.1.0/lib/libsundials_idas.a lib/idas-2.1.0/lib/libsundials_nvecserial.a sandbox/idas_test.cpp   -o sandbox/idas_test

};


// TODO
template<typename F, typename TYY, typename TYP, typename TPAR>
class idas_adjoint_system: public idas_system<F, TYY, TYP, TPAR> {
};

}
}

#endif
