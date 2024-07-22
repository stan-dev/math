#ifndef STAN_MATH_REV_FUNCTOR_DAE_RESIDUAL_HPP
#define STAN_MATH_REV_FUNCTOR_DAE_RESIDUAL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/nested_rev_autodiff.hpp>
#include <stan/math/rev/core/zero_adjoints.hpp>
#include <stan/math/rev/core/accumulate_adjoints.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <ostream>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/**
 * IDAS DAE system that contains information on residual
 * equation functor, sensitivity residual equation functor,
 * as well as initial conditions. This is a base type that
 * is intended to contain common values used by forward
 * sensitivity system.
 *
 * @tparam F type of functor for DAE residual
 * @tparam Tyy scalar type of initial unknown values
 * @tparam Typ scalar type of initial unknown's derivative values
 * @tparam Tpar scalar type of parameters
 */
template <typename F, typename Tyy, typename Typ, typename... T_par>
class dae_system {
  using dae_type = dae_system<F, Tyy, Typ, T_par...>;

 protected:
  const F& f_; /**< residual functor */
  std::tuple<decltype(deep_copy_vars(std::declval<const T_par&>()))...>
      local_args_tuple_;
  std::tuple<plain_type_t<decltype(value_of(std::declval<const T_par&>()))>...>
      dbl_args_tuple_;
  std::tuple<const T_par&...> args_tuple_;
  std::ostream* msgs_;

 public:
  Tyy const& yy; /**< state variable y */
  Typ const& yp; /**< time derivatives of y */
  Eigen::VectorXd dbl_yy;
  Eigen::VectorXd dbl_yp;
  const size_t N;
  const size_t M;
  const size_t ns;
  vari** varis;
  std::vector<stan::math::var> all_vars;
  Eigen::VectorXd dbl_rr;

  static constexpr bool is_var_yy0 = stan::is_var<return_type_t<Tyy>>::value;
  static constexpr bool is_var_yp0 = stan::is_var<return_type_t<Typ>>::value;
  static constexpr bool is_var_par
      = stan::is_var<return_type_t<T_par...>>::value;
  static constexpr bool use_fwd_sens = is_var_yy0 || is_var_yp0 || is_var_par;

  using scalar_t = stan::return_type_t<Tyy, Typ, T_par...>;
  using return_t = std::vector<Eigen::Matrix<scalar_t, -1, 1>>;

  /**
   * Construct IDAS DAE system from initial condition and parameters
   *
   * @param[in] f DAE residual functor
   * @param[in] eq_id array for DAE's variable ID(1 for *
   *                  derivative variables, 0 for algebraic variables).
   * @param[in] yy0 initial condition
   * @param[in] yp0 initial condition for derivatives
   * @param[in, out] msgs the print stream for warning messages
   * @param args Extra arguments passed unmodified through to DAE right hand
   * side
   */
  dae_system(const F& f, const Tyy& yy0, const Typ& yp0, std::ostream* msgs,
             const T_par&... args)
      : f_(f),
        local_args_tuple_(deep_copy_vars(args)...),
        dbl_args_tuple_(value_of(args)...),
        args_tuple_(std::forward_as_tuple(args...)),
        msgs_(msgs),
        yy(yy0),
        yp(yp0),
        dbl_yy(stan::math::value_of(yy0)),
        dbl_yp(stan::math::value_of(yp0)),
        N(yy0.size()),
        M(count_vars(args...)),
        ns((is_var_yy0 ? N : 0) + (is_var_yp0 ? N : 0) + M),
        varis(ChainableStack::instance_->memalloc_.alloc_array<vari*>(ns)),
        all_vars(ns),
        dbl_rr(N) {
    save_varis(varis, yy0, yp0, args...);
    for (int i = 0; i < ns; ++i) {
      all_vars[i].vi_ = *(varis + i);
    }
  }

  void eval_residual(double t) {
    dbl_rr = math::apply(
        [&](auto&&... args) { return f_(t, dbl_yy, dbl_yp, msgs_, args...); },
        dbl_args_tuple_);
  }

  /**
   * Evaluate DAE residual according to IDAS signature
   */
  static int idas_res(double t, N_Vector yy, N_Vector yp, N_Vector rr,
                      void* user_data) {
    dae_type* dae = static_cast<dae_type*>(user_data);
    dae->dbl_yy = Eigen::Map<Eigen::VectorXd>(NV_DATA_S(yy), dae->N);
    dae->dbl_yp = Eigen::Map<Eigen::VectorXd>(NV_DATA_S(yp), dae->N);
    dae->eval_residual(t);
    Eigen::Map<Eigen::VectorXd>(NV_DATA_S(rr), dae->N) = dae->dbl_rr;

    return 0;
  }

  /**
   * Evaluate DAE sensitivity residual according to IDAS signature
   */
  static int idas_sens_res(int ns, double t, N_Vector yy, N_Vector yp,
                           N_Vector res, N_Vector* yys, N_Vector* yps,
                           N_Vector* ress, void* user_data, N_Vector temp1,
                           N_Vector temp2, N_Vector temp3) {
    dae_type* dae = static_cast<dae_type*>(user_data);

    const int n = dae->N;
    const int m = dae->M;

    for (int i = 0; i < ns; ++i) {
      N_VConst(0.0, ress[i]);
    }

    // Run nested autodiff in this scope
    stan::math::nested_rev_autodiff nested;

    Eigen::Matrix<var, -1, 1> yy_var(n), yp_var(n);
    for (int i = 0; i < n; ++i) {
      yy_var.coeffRef(i) = NV_Ith_S(yy, i);
      yp_var.coeffRef(i) = NV_Ith_S(yp, i);
    }

    Eigen::VectorXd g(m);
    Eigen::Matrix<var, -1, 1> fy = math::apply(
        [&](auto&&... args) {
          return dae->f_(t, yy_var, yp_var, dae->msgs_, args...);
        },
        dae->local_args_tuple_);

    for (int i = 0; i < n; ++i) {
      if (i > 0) {
        nested.set_zero_all_adjoints();
      }
      fy[i].grad();

      for (int j = 0; j < ns; ++j) {
        auto yysp = N_VGetArrayPointer(yys[j]);
        auto ypsp = N_VGetArrayPointer(yps[j]);
        auto ressp = N_VGetArrayPointer(ress[j]);
        for (int k = 0; k < n; ++k) {
          ressp[i] += yy_var[k].adj() * yysp[k] + yp_var[k].adj() * ypsp[k];
        }
      }

      if (is_var_par) {
        g.fill(0);
        stan::math::apply(
            [&](auto&&... args) {
              stan::math::accumulate_adjoints(g.data(), args...);
            },
            dae->local_args_tuple_);
        for (int j = 0; j < m; ++j) {
          auto ressp = N_VGetArrayPointer(ress[ns - m + j]);
          ressp[i] += g[j];
        }

        stan::math::for_each([](auto&& arg) { stan::math::zero_adjoints(arg); },
                             dae->local_args_tuple_);
      }
    }
    return 0;
  }
};

}  // namespace math
}  // namespace stan

#endif
