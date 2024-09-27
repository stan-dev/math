
#ifndef STAN_MATH_FORWARD_DECL_HPP
#define STAN_MATH_FORWARD_DECL_HPP
#include <type_traits>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
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
} // namespace internal

// ./stan/math/prim/meta/ad_promotable.hpp:18
template <typename From, typename To, typename = void>
struct ad_promotable;

// ./stan/math/prim/meta/conjunction.hpp:13
template <typename... T>
struct conjunction;

// ./stan/math/prim/meta/disjunction.hpp:13
template <typename... Conds>
struct disjunction;

// ./stan/math/prim/meta/append_return_type.hpp:21
template <typename T1, typename T2>
struct append_return_type;

// ./stan/math/prim/meta/append_return_type.hpp:35
template <>
struct append_return_type<int, int>;

// ./stan/math/prim/meta/holder.hpp:70
template <class ArgType, typename... Ptrs>
class Holder;

// ./stan/math/prim/meta/holder.hpp:114
template <typename ArgType, typename... Ptrs>
class Holder;

// ./stan/math/prim/meta/include_summand.hpp:36
template <bool propto, typename T = double, typename... T_pack>
struct include_summand;

// ./stan/math/prim/meta/index_type.hpp:24
template <typename T, typename = void>
struct index_type;

// ./stan/math/prim/meta/is_matrix_cl.hpp:10
template <typename T>
class arena_matrix_cl;

// ./stan/math/prim/meta/is_matrix_cl.hpp:17
class matrix_cl_base;

// ./stan/math/prim/meta/is_tuple.hpp:19
template <typename T>
struct is_tuple;

// ./stan/math/prim/meta/promote_scalar_type.hpp:21
template <typename T, typename S, typename Enable = void>
struct promote_scalar_type;

// ./stan/math/prim/meta/seq_view.hpp:10
template <typename T>
struct store_type;

// ./stan/math/prim/meta/seq_view.hpp:14
template <>
struct store_type<double>;

// ./stan/math/prim/meta/seq_view.hpp:18
template <>
struct store_type<int>;

// ./stan/math/prim/meta/seq_view.hpp:23
template <typename T>
struct pass_type;

// ./stan/math/prim/meta/seq_view.hpp:27
template <>
struct pass_type<double>;

// ./stan/math/prim/meta/seq_view.hpp:31
template <>
struct pass_type<int>;

// ./stan/math/prim/meta/seq_view.hpp:37
template <typename T, typename S>
class seq_view;

// ./stan/math/fwd/core/fvar.hpp:39
template <typename T>
struct fvar;

// ./stan/math/prim/functor/apply_scalar_unary.hpp:40
template <typename F, typename T, typename Enable = void>
struct apply_scalar_unary;

// ./stan/math/prim/functor/apply_vector_unary.hpp:16
template <typename T, typename Enable = void>
struct apply_vector_unary;

// ./stan/math/prim/fun/square.hpp:35
struct square_fun;

// ./stan/math/prim/core/complex_base.hpp:17
template <typename ValueType>
class complex_base;

// ./stan/math/memory/stack_alloc.hpp:72
class stack_alloc;

// ./stan/math/rev/core/autodiffstackstorage.hpp:88
template <typename ChainableT, typename ChainableAllocT>
struct AutodiffStackSingleton;

// ./stan/math/rev/core/chainablestack.hpp:7
class chainable_alloc;

// ./stan/math/rev/core/chainablestack.hpp:8
class vari_base;

// ./stan/math/rev/core/arena_allocator.hpp:22
template <typename T>
struct arena_allocator;

// ./stan/math/rev/core/chainable_alloc.hpp:16
class chainable_alloc;

// ./stan/math/rev/core/var_value_fwd_declare.hpp:7
template <typename T, typename = void>
class var_value;

// ./stan/math/rev/core/var_value_fwd_declare.hpp:17
template <typename MatrixType, typename = void>
class arena_matrix;

// ./stan/math/rev/core/chainable_object.hpp:24
template <typename T>
class chainable_object;

// ./stan/math/rev/core/chainable_object.hpp:77
template <typename T>
class unsafe_chainable_object;

// ./stan/math/rev/core/vari.hpp:16
template <typename T, typename = void>
class vari_value;

// ./stan/math/rev/core/vari.hpp:28
class vari_base;

// ./stan/math/rev/core/vari.hpp:205
template <typename T, typename = void>
class vari_view;

// ./stan/math/rev/core/vari.hpp:214
template <typename Derived>
class vari_view_eigen;

// ./stan/math/rev/core/init_chainablestack.hpp:27
class ad_tape_observer final;

// ./stan/math/rev/core/ddv_vari.hpp:9
class op_ddv_vari;

// ./stan/math/rev/core/dv_vari.hpp:9
class op_dv_vari;

// ./stan/math/rev/core/dvd_vari.hpp:9
class op_dvd_vari;

// ./stan/math/rev/core/dvv_vari.hpp:9
class op_dvv_vari;

// ./stan/math/rev/core/gevv_vvv_vari.hpp:11
class gevv_vvv_vari;

// ./stan/math/rev/core/read_var.hpp:17
template <typename EigRev, typename EigVari, typename EigDbl>
class vi_val_adj_functor;

// ./stan/math/rev/core/read_var.hpp:44
template <typename EigRev, typename EigDbl>
class val_adj_functor;

// ./stan/math/rev/core/read_var.hpp:82
template <typename EigVar, typename EigVari>
class vi_val_functor;

// ./stan/math/rev/core/read_var.hpp:106
template <typename EigVar, typename EigVari>
class vi_adj_functor;

// ./stan/math/rev/core/nested_rev_autodiff.hpp:27
class nested_rev_autodiff;

// ./stan/math/rev/core/matrix_vari.hpp:15
class op_matrix_vari;

// ./stan/math/prim/fun/inv.hpp:19
struct inv_fun;

// ./stan/math/rev/core/vv_vari.hpp:9
class op_vv_vari;

// ./stan/math/rev/core/vd_vari.hpp:9
class op_vd_vari;

// ./stan/math/rev/core/precomp_vv_vari.hpp:11
class precomp_vv_vari final;

// ./stan/math/rev/core/vvv_vari.hpp:9
class op_vvv_vari;

// ./stan/math/rev/core/precomp_vvv_vari.hpp:11
class precomp_vvv_vari final;

// ./stan/math/rev/core/precomputed_gradients.hpp:27
template <typename ContainerOperands = std::tuple<>,
          typename ContainerGradients = std::tuple<>>
class precomputed_gradients_vari_template;

// ./stan/math/prim/fun/fabs.hpp:30
struct fabs_fun;

// ./stan/math/prim/fun/floor.hpp:20
struct floor_fun;

// ./stan/math/prim/fun/abs.hpp:50
struct abs_fun;

// ./stan/math/prim/fun/log.hpp:22
struct log_fun;

// ./stan/math/prim/fun/LDLT_factor.hpp:20
template <typename T, typename Enable = void>
class LDLT_factor;

// ./stan/math/rev/core/profiling.hpp:21
class profile_info;

// ./stan/math/rev/core/profiling.hpp:150
template <typename T>
class profile;

// ./stan/math/rev/core/scoped_chainablestack.hpp:31
class ScopedChainableStack;

// ./stan/math/rev/core/stored_gradient_vari.hpp:17
class stored_gradient_vari;

// ./stan/math/rev/core/vdd_vari.hpp:9
class op_vdd_vari;

// ./stan/math/rev/core/vdv_vari.hpp:9
class op_vdv_vari;

// ./stan/math/rev/core/vector_vari.hpp:12
class op_vector_vari;

// ./stan/math/rev/core/vvd_vari.hpp:9
class op_vvd_vari;

// ./stan/math/fwd/fun/read_fvar.hpp:13
template <typename EigFvar, typename EigOut>
class read_fvar_functor;

// ./stan/math/prim/fun/sqrt.hpp:21
struct sqrt_fun;

// ./stan/math/prim/fun/inv_sqrt.hpp:27
struct inv_sqrt_fun;

// ./stan/math/prim/fun/accumulator.hpp:23
template <typename T, typename = void>
class accumulator;

// ./stan/math/prim/fun/exp.hpp:18
struct exp_fun;

// ./stan/math/prim/fun/cosh.hpp:22
struct cosh_fun;

// ./stan/math/prim/fun/cos.hpp:23
struct cos_fun;

// ./stan/math/prim/fun/asinh.hpp:30
struct asinh_fun;

// ./stan/math/prim/fun/asin.hpp:26
struct asin_fun;

// ./stan/math/prim/fun/acos.hpp:28
struct acos_fun;

// ./stan/math/prim/fun/acosh.hpp:62
struct acosh_fun;

// ./stan/math/prim/fun/atanh.hpp:51
struct atanh_fun;

// ./stan/math/prim/fun/atan.hpp:24
struct atan_fun;

// ./stan/math/prim/fun/lgamma.hpp:101
struct lgamma_fun;

// ./stan/math/prim/fun/digamma.hpp:59
struct digamma_fun;

// ./stan/math/prim/fun/log1p.hpp:55
struct log1p_fun;

// ./stan/math/prim/fun/log1m.hpp:56
struct log1m_fun;

// ./stan/math/prim/fun/cbrt.hpp:18
struct cbrt_fun;

// ./stan/math/prim/fun/sinh.hpp:21
struct sinh_fun;

// ./stan/math/prim/fun/sin.hpp:23
struct sin_fun;

// ./stan/math/prim/fun/trigamma.hpp:134
struct trigamma_fun;

// ./stan/math/prim/fun/erf.hpp:18
struct erf_fun;

// ./stan/math/prim/fun/erfc.hpp:19
struct erfc_fun;

// ./stan/math/prim/fun/exp2.hpp:14
struct exp2_fun;

// ./stan/math/prim/fun/expm1.hpp:18
struct expm1_fun;

// ./stan/math/prim/fun/tgamma.hpp:34
struct tgamma_fun;

// ./stan/math/prim/fun/sign.hpp:19
struct sign_fun;

// ./stan/math/prim/fun/inv_erfc.hpp:31
struct inv_erfc_fun;

// ./stan/math/prim/fun/Phi.hpp:52
struct Phi_fun;

// ./stan/math/prim/fun/inv_Phi.hpp:161
struct inv_Phi_fun;

// ./stan/math/prim/fun/inv_cloglog.hpp:60
struct inv_cloglog_fun;

// ./stan/math/prim/fun/log1m_exp.hpp:66
struct log1m_exp_fun;

// ./stan/math/prim/fun/inv_logit.hpp:70
struct inv_logit_fun;

// ./stan/math/prim/fun/log10.hpp:22
struct log10_fun;

// ./stan/math/prim/fun/log1p_exp.hpp:61
struct log1p_exp_fun;

// ./stan/math/prim/fun/log1m_inv_logit.hpp:58
struct log1m_inv_logit_fun;

// ./stan/math/prim/fun/log2.hpp:22
struct log2_fun;

// ./stan/math/prim/fun/log_inv_logit.hpp:56
struct log_inv_logit_fun;

// ./stan/math/prim/functor/operands_and_partials.hpp:15
template <typename Op1 = double, typename Op2 = double, typename Op3 = double,
          typename Op4 = double, typename Op5 = double, typename Op6 = double,
          typename Op7 = double, typename Op8 = double,
          typename T_return_type
          = return_type_t<Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8>>
class operands_and_partials;

// ./stan/math/prim/functor/operands_and_partials.hpp:151
template <typename Op1, typename Op2, typename Op3, typename Op4, typename Op5,
          typename Op6, typename Op7, typename Op8, typename T_return_type>
class operands_and_partials;

// ./stan/math/prim/fun/logit.hpp:62
struct logit_fun;

// ./stan/math/prim/fun/round.hpp:20
struct round_fun;

// ./stan/math/prim/fun/tanh.hpp:23
struct tanh_fun;

// ./stan/math/prim/fun/tan.hpp:23
struct tan_fun;

// ./stan/math/prim/fun/trunc.hpp:14
struct trunc_fun;

// ./stan/math/prim/fun/serializer.hpp:21
template <typename T>
struct deserializer;

// ./stan/math/prim/fun/serializer.hpp:180
template <typename T>
struct serializer;

// ./stan/math/prim/prob/std_normal_log_qf.hpp:133
struct std_normal_log_qf_fun;

// ./stan/math/prim/fun/promote_elements.hpp:20
template <typename T, typename S>
struct promote_elements;

// ./stan/math/prim/fun/array_builder.hpp:19
template <typename T>
class array_builder;

// ./stan/math/prim/fun/ceil.hpp:20
struct ceil_fun;

// ./stan/math/prim/fun/matrix_exp_action_handler.hpp:25
class matrix_exp_action_handler;

// ./stan/math/prim/fun/Phi_approx.hpp:40
struct Phi_approx_fun;

// ./stan/math/prim/fun/step.hpp:38
struct step_fun;

// ./stan/math/prim/fun/to_int.hpp:57
struct to_int_fun;

// ./stan/math/prim/fun/welford_covar_estimator.hpp:10
class welford_covar_estimator;

// ./stan/math/prim/fun/welford_var_estimator.hpp:10
class welford_var_estimator;

// ./stan/math/rev/fun/cov_exp_quad.hpp:21
template <typename T_x, typename T_sigma, typename T_l>
class cov_exp_quad_vari;

// ./stan/math/prim/functor/coupled_ode_system.hpp:14
template <bool ArithmeticArguments, typename F, typename T_y0, typename... Args>
struct coupled_ode_system_impl;

// ./stan/math/prim/functor/coupled_ode_system.hpp:104
template <typename F, typename T_y0, typename... Args>
struct coupled_ode_system;

// ./stan/math/rev/functor/algebra_system.hpp:23
template <typename T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct nlo_functor;

// ./stan/math/rev/functor/algebra_system.hpp:44
template <typename S>
struct hybrj_functor_solver;

// ./stan/math/rev/functor/kinsol_data.hpp:28
template <typename F1, typename... Args>
class kinsol_system_data;

// ./stan/math/rev/functor/algebra_solver_fp.hpp:34
template <typename F>
struct KinsolFixedPointEnv;

// ./stan/math/rev/functor/algebra_solver_fp.hpp:153
struct FixedPointADJac;

// ./stan/math/rev/functor/algebra_solver_fp.hpp:220
template <typename fp_env_type, typename fp_jac_type>
struct FixedPointSolver;

// ./stan/math/rev/functor/cvodes_integrator.hpp:33
template <int Lmm, typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args>
class cvodes_integrator;

// ./stan/math/rev/functor/idas_service.hpp:31
template <typename dae_type>
struct idas_service;

// ./stan/math/rev/functor/idas_integrator.hpp:22
class idas_integrator;

// ./stan/math/rev/functor/dae_system.hpp:35
template <typename F, typename Tyy, typename Typ, typename... T_par>
class dae_system;

// ./stan/math/rev/functor/cvodes_integrator_adjoint.hpp:35
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args>
class cvodes_integrator_adjoint_vari;

// ./stan/math/prim/meta/holder.hpp:238
template <typename T,
  typename... Ptrs,
  std::enable_if_t<sizeof...(Ptrs)>=1 >*>
Holder <T ,Ptrs ...>holder (T &&arg , Ptrs *...pointers );

// ./stan/math/prim/meta/holder.hpp:244
template <typename T>
T holder (T &&arg );

// ./stan/math/prim/meta/holder.hpp:349
template <typename F,
  typename... Args,
  require_not_plain_type_t<decltype (std::declval<F>()(std::declval<Args&>()...))>*>
auto make_holder (const F &func , Args &&...args );

// ./stan/math/prim/meta/holder.hpp:368
template <typename F,
  typename... Args,
  require_plain_type_t<decltype (std::declval<F>()(std::declval<Args&>()...))>*>
auto make_holder (const F &func , Args &&...args );

// ./stan/math/prim/meta/index_apply.hpp:25
template <std::size_tN,
  class F>
inline constexpr auto index_apply (F &&f );

// ./stan/math/prim/fun/sum.hpp:20
template <typename T,
  require_stan_scalar_t<T>*>
inline T sum (T &&m );

// ./stan/math/prim/fun/sum.hpp:32
template <typename T,
  require_not_var_t<T>*>
inline T sum (const std ::vector <T >&m );

// ./stan/math/prim/fun/sum.hpp:45
template <typename T,
  require_eigen_vt<std::is_arithmetic,
  T>*>
inline value_type_t <T >sum (const T &m );

// ./stan/math/prim/fun/sum.hpp:58
template <typename T,
  require_eigen_vt<is_complex,
  T>*>
inline value_type_t <T >sum (const T &m );

// ./stan/math/prim/meta/possibly_sum.hpp:17
template <typename CondSum,
  typename T,
  require_t<CondSum>*>
inline auto possibly_sum (T &&x );

// ./stan/math/prim/meta/possibly_sum.hpp:29
template <typename CondSum,
  typename T1,
  require_not_t<CondSum>*>
inline auto possibly_sum (T1 &&x );

// ./stan/math/prim/meta/static_select.hpp:22
template <bool Condition,
  typename T1,
  typename T2,
  std::enable_if_t<Condition>*>
T1 static_select (T1 &&a , T2 &&b );

// ./stan/math/prim/meta/static_select.hpp:28
template <bool Condition,
  typename T1,
  typename T2,
  std::enable_if_t<!Condition>*>
T2 static_select (T1 &&a , T2 &&b );

// ./stan/math/prim/fun/is_nan.hpp:19
template <typename T,
  require_integral_t<T>*>
inline bool is_nan (T x );

// ./stan/math/prim/fun/is_nan.hpp:33
template <typename T,
  require_floating_point_t<T>*>
inline bool is_nan (T x );

// ./stan/math/prim/fun/is_nan.hpp:38
template <typename T,
  require_eigen_t<T>*>
inline bool is_nan (const T &x );

// ./stan/math/fwd/core/operator_addition.hpp:17
template <typename T>
inline fvar <T >operator +(const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/core/operator_addition.hpp:30
template <typename T>
inline fvar <T >operator +(double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/core/operator_addition.hpp:43
template <typename T>
inline fvar <T >operator +(const fvar <T >&x1 , double x2 );

// ./stan/math/prim/core/operator_division.hpp:39
template <typename U,
  typename V,
  require_all_stan_scalar_t<U,
  V>*>
inline complex_return_t <U ,V >operator /(const std ::complex <U >&x , const std ::complex <V >&y );

// ./stan/math/prim/core/operator_division.hpp:54
template <typename U,
  typename V,
  require_all_stan_scalar_t<U,
  V>*>
inline complex_return_t <U ,V >operator /(const std ::complex <U >&x , const V &y );

// ./stan/math/prim/core/operator_division.hpp:68
template <typename U,
  typename V,
  require_all_stan_scalar_t<U,
  V>*>
inline complex_return_t <U ,V >operator /(const U &x , const std ::complex <V >&y );

// ./stan/math/fwd/core/operator_division.hpp:20
template <typename T>
inline fvar <T >operator /(const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/core/operator_division.hpp:34
template <typename T,
  typename U,
  require_arithmetic_t<U>*>
inline fvar <T >operator /(const fvar <T >&x1 , U x2 );

