
#ifndef STAN_MATH_MANUAL_FORWARD_DECL_HPP
#define STAN_MATH_MANUAL_FORWARD_DECL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <type_traits>
#include <vector>
#include <boost/optional.hpp>
#include <unsupported/Eigen/FFT>
#ifndef TBB_INTERFACE_NEW
#include <tbb/tbb_stddef.h>

#if TBB_VERSION_MAJOR >= 2020
#define TBB_INTERFACE_NEW
#endif
#endif

#ifdef TBB_INTERFACE_NEW
#include <tbb/global_control.h>
#include <tbb/task_arena.h>
#else
#include <tbb/task_scheduler_init.h>
#endif


namespace stan {
namespace math {
namespace internal {
struct nonexisting_adjoint;
template <typename T, typename S, typename Enable>
class empty_broadcast_array;
template <typename T, typename F>
struct callback_vari;
template <typename ReturnType, typename Enable, typename... Ops>
class partials_propagator;
template <typename T, typename, typename>
struct arena_type_impl;
} // namespace internal

template <typename T>
struct fvar;

template <typename T, typename = void>
class vari_value;

using vari = vari_value<double, void>;

// ./stan/math/rev/core/var_value_fwd_declare.hpp:7
template <typename T, typename = void>
class var_value;
using var = var_value<double, void>;

template <typename F>
void finite_diff_grad_hessian(const F& f, const Eigen::VectorXd& x, double& fx,
                              Eigen::MatrixXd& hess,
                              std::vector<Eigen::MatrixXd>& grad_hess_fx,
                              double epsilon = 1e-04);

// ./stan/math/rev/core/var_value_fwd_declare.hpp:17
/**
 * Equivalent to `Eigen::Matrix`, except that the data is stored on AD stack.
 * That makes these objects trivially destructible and usable in `vari`s.
 *
 * @tparam MatrixType Eigen matrix type this works as (`MatrixXd`, `VectorXd`,
 * ...)
 */
template <typename MatrixType, typename>
class arena_matrix;

template <typename T, typename Enable = void>
class LDLT_factor;

template <typename T, typename = void>
class accumulator;
#ifdef TBB_INTERFACE_NEW
inline auto& init_threadpool_tbb(int n_threads = 0);
#else
inline auto& init_threadpool_tbb(int n_threads = 0);
#endif

// ./stan/math/prim/fun/grad_reg_inc_gamma.hpp:51
template <typename T1, typename T2>
return_type_t<T1, T2> grad_reg_inc_gamma(T1 a, T2 z, T1 g, T1 dig,
                                         double precision = 1e-6,
                                         int max_steps = 1e5);

// ./stan/math/prim/fun/grad_reg_lower_inc_gamma.hpp:109
template <typename T1, typename T2>
return_type_t<T1, T2> grad_reg_lower_inc_gamma(const T1& a, const T2& z,
                                               double precision = 1e-10,
                                               int max_steps = 1e5);

// ./stan/math/prim/fun/grad_pFq.hpp:87
template <bool calc_a, bool calc_b, bool calc_z,
          typename TpFq, typename Ta, typename Tb, typename Tz,
          typename T_Rtn,
          typename Ta_Rtn,
          typename Tb_Rtn>
std::tuple<Ta_Rtn, Tb_Rtn, T_Rtn> grad_pFq(const TpFq& pfq_val, const Ta& a,
                                           const Tb& b, const Tz& z,
                                           double precision = 1e-14,
                                           int max_steps = 1e6);
template <typename T1, typename T2, typename T3, typename T_z>
auto grad_2F1(const T1& a1, const T2& a2, const T3& b1, const T_z& z,
              double precision = 1e-14, int max_steps = 1e6);
}
}




#endif
