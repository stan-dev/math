
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
    typename ...Ptrs,
    std::enable_if_t<sizeof ...(Ptrs)>=1 >*>
Holder <T ,Ptrs ...>holder (T &&arg , Ptrs *...pointers );

// ./stan/math/prim/meta/holder.hpp:244
template <typename T>
T holder (T &&arg );

// ./stan/math/prim/meta/holder.hpp:349
template <typename F,
    typename ...Args,
    require_not_plain_type_t<decltype (std::declval<F>()(std::declval<Args&>()...))>*>
auto make_holder (const F &func , Args &&...args );

// ./stan/math/prim/meta/holder.hpp:368
template <typename F,
    typename ...Args,
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

// ./stan/math/fwd/core/operator_division.hpp:48
template <typename T,
    typename U,
    require_arithmetic_t<U>*>
inline fvar <T >operator /(U x1 , const fvar <T >&x2 );

// ./stan/math/fwd/core/operator_division.hpp:54
template <typename T>
inline std ::complex <fvar <T >>operator /(const std ::complex <fvar <T >>&x1 , const std ::complex <fvar <T >>&x2 );

// ./stan/math/fwd/core/operator_division.hpp:59
template <typename T,
    typename U,
    require_arithmetic_t<U>*>
inline std ::complex <fvar <T >>operator /(const std ::complex <fvar <T >>&x1 , const std ::complex <U >&x2 );

// ./stan/math/fwd/core/operator_division.hpp:64
template <typename T>
inline std ::complex <fvar <T >>operator /(const std ::complex <fvar <T >>&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/core/operator_division.hpp:69
template <typename T,
    typename U,
    require_arithmetic_t<U>*>
inline std ::complex <fvar <T >>operator /(const std ::complex <fvar <T >>&x1 , U x2 );

// ./stan/math/fwd/core/operator_division.hpp:74
template <typename T,
    typename U,
    require_arithmetic_t<U>*>
inline std ::complex <fvar <T >>operator /(const std ::complex <U >&x1 , const std ::complex <fvar <T >>&x2 );

// ./stan/math/fwd/core/operator_division.hpp:79
template <typename T,
    typename U,
    require_arithmetic_t<U>*>
inline std ::complex <fvar <T >>operator /(const std ::complex <U >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/core/operator_division.hpp:85
template <typename T>
inline std ::complex <fvar <T >>operator /(const fvar <T >&x1 , const std ::complex <fvar <T >>&x2 );

// ./stan/math/fwd/core/operator_division.hpp:97
template <typename T,
    typename U,
    require_arithmetic_t<U>*>
inline std ::complex <fvar <T >>operator /(U x1 , const std ::complex <fvar <T >>&x2 );

// ./stan/math/fwd/core/operator_equal.hpp:18
template <typename T>
inline bool operator ==(const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_equal.hpp:33
template <typename T>
inline bool operator ==(const fvar <T >&x , double y );

// ./stan/math/fwd/core/operator_equal.hpp:47
template <typename T>
inline bool operator ==(double x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_greater_than.hpp:19
template <typename T>
inline bool operator >(const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_greater_than.hpp:34
template <typename T>
inline bool operator >(const fvar <T >&x , double y );

// ./stan/math/fwd/core/operator_greater_than.hpp:49
template <typename T>
inline bool operator >(double x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_greater_than_or_equal.hpp:19
template <typename T>
inline bool operator >=(const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_greater_than_or_equal.hpp:35
template <typename T>
inline bool operator >=(const fvar <T >&x , double y );

// ./stan/math/fwd/core/operator_greater_than_or_equal.hpp:51
template <typename T>
inline bool operator >=(double x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_less_than.hpp:19
template <typename T>
inline bool operator <(const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_less_than.hpp:34
template <typename T>
inline bool operator <(double x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_less_than.hpp:50
template <typename T>
inline bool operator <(const fvar <T >&x , double y );

// ./stan/math/fwd/core/operator_less_than_or_equal.hpp:20
template <typename T>
inline bool operator <=(const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_less_than_or_equal.hpp:35
template <typename T>
inline bool operator <=(const fvar <T >&x , double y );

// ./stan/math/fwd/core/operator_less_than_or_equal.hpp:50
template <typename T>
inline bool operator <=(double x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_logical_and.hpp:18
template <typename T>
inline bool operator &&(const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_logical_and.hpp:33
template <typename T>
inline bool operator &&(const fvar <T >&x , double y );

// ./stan/math/fwd/core/operator_logical_and.hpp:48
template <typename T>
inline bool operator &&(double x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_logical_or.hpp:18
template <typename T>
inline bool operator ||(const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_logical_or.hpp:33
template <typename T>
inline bool operator ||(const fvar <T >&x , double y );

// ./stan/math/fwd/core/operator_logical_or.hpp:48
template <typename T>
inline bool operator ||(double x , const fvar <T >&y );

// ./stan/math/prim/fun/as_column_vector_or_scalar.hpp:21
template <typename T,
    require_stan_scalar_t<T>*>
inline T as_column_vector_or_scalar (const T &a );

// ./stan/math/prim/fun/as_column_vector_or_scalar.hpp:31
template <typename T,
    typename S>
internal ::empty_broadcast_array <T ,S ,void >&as_column_vector_or_scalar (internal ::empty_broadcast_array <T , S , void >&a );

// ./stan/math/prim/fun/as_column_vector_or_scalar.hpp:43
template <typename T,
    require_eigen_col_vector_t<T>*>
inline T &&as_column_vector_or_scalar (T &&a );

// ./stan/math/prim/fun/as_column_vector_or_scalar.hpp:57
template <typename T,
    require_eigen_row_vector_t<T>*>
inline auto as_column_vector_or_scalar (T &&a );

// ./stan/math/prim/fun/as_column_vector_or_scalar.hpp:70
template <typename T,
    require_std_vector_t<T>*>
inline auto as_column_vector_or_scalar (T &&a );

// ./stan/math/prim/fun/square.hpp:26

inline double square (double x );

// ./stan/math/prim/fun/square.hpp:49
template <typename Container,
    require_not_stan_scalar_t<Container>*>
inline auto square (const Container &x );

// ./stan/math/prim/fun/square.hpp:66
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto square (const Container &x );

// ./stan/math/prim/core/operator_multiplication.hpp:40
template <typename U,
    typename V,
    require_all_stan_scalar_t<U,
    V>*>
inline complex_return_t <U ,V >operator *(const std ::complex <U >&x , const std ::complex <V >&y );

// ./stan/math/prim/core/operator_multiplication.hpp:55
template <typename U,
    typename V,
    require_stan_scalar_t<V>*>
inline complex_return_t <U ,V >operator *(const std ::complex <U >&x , const V &y );

// ./stan/math/prim/core/operator_multiplication.hpp:69
template <typename U,
    typename V,
    require_stan_scalar_t<U>*>
inline complex_return_t <U ,V >operator *(const U &x , const std ::complex <V >&y );

// ./stan/math/fwd/core/operator_multiplication.hpp:19
template <typename T>
inline fvar <T >operator *(const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_multiplication.hpp:32
template <typename T>
inline fvar <T >operator *(double x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_multiplication.hpp:45
template <typename T>
inline fvar <T >operator *(const fvar <T >&x , double y );

// ./stan/math/fwd/core/operator_multiplication.hpp:58
template <typename T>
inline std ::complex <stan ::math ::fvar <T >>operator *(const std ::complex <stan ::math ::fvar <T >>&x , const std ::complex <stan ::math ::fvar <T >>&y );

// ./stan/math/fwd/core/operator_multiplication.hpp:74
template <typename T>
inline std ::complex <stan ::math ::fvar <T >>operator *(const std ::complex <double >&x , const std ::complex <stan ::math ::fvar <T >>&y );

// ./stan/math/fwd/core/operator_multiplication.hpp:89
template <typename T>
inline std ::complex <stan ::math ::fvar <T >>operator *(const std ::complex <stan ::math ::fvar <T >>&x , const std ::complex <double >&y );

// ./stan/math/fwd/core/operator_not_equal.hpp:18
template <typename T>
inline bool operator !=(const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_not_equal.hpp:33
template <typename T>
inline bool operator !=(const fvar <T >&x , double y );

// ./stan/math/fwd/core/operator_not_equal.hpp:48
template <typename T>
inline bool operator !=(double x , const fvar <T >&y );

// ./stan/math/fwd/core/operator_subtraction.hpp:17
template <typename T>
inline fvar <T >operator -(const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/core/operator_subtraction.hpp:30
template <typename T>
inline fvar <T >operator -(double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/core/operator_subtraction.hpp:43
template <typename T>
inline fvar <T >operator -(const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/core/operator_unary_minus.hpp:16
template <typename T>
inline fvar <T >operator -(const fvar <T >&x );

// ./stan/math/fwd/core/operator_unary_not.hpp:17
template <typename T>
inline bool operator !(const fvar <T >&x );

// ./stan/math/fwd/core/operator_unary_plus.hpp:18
template <typename T>
inline fvar <T >operator +(const fvar <T >&x );

// ./stan/math/memory/stack_alloc.hpp:28
template <typename T>
bool is_aligned (T *ptr , unsigned int bytes_aligned );

// ./stan/math/rev/core/chainable_object.hpp:57
template <typename T>
auto make_chainable_ptr (T &&obj );

// ./stan/math/rev/core/chainable_object.hpp:112
template <typename T>
auto make_unsafe_chainable_ptr (T &&obj );

// ./stan/math/prim/fun/to_ref.hpp:16
template <typename T>
inline ref_type_t <T &&>to_ref (T &&a );

// ./stan/math/prim/fun/to_ref.hpp:28
template <bool Cond,
    typename T,
    std::enable_if_t<!Cond>*>
inline T to_ref_if (T &&a );

// ./stan/math/prim/fun/to_ref.hpp:42
template <bool Cond,
    typename T,
    std::enable_if_t<Cond>*>
inline ref_type_t <T &&>to_ref_if (T &&a );

// ./stan/math/rev/core/empty_nested.hpp:12

static inline bool empty_nested ();

// ./stan/math/rev/core/nested_size.hpp:10

static inline size_t nested_size ();

// ./stan/math/rev/core/grad.hpp:26

static inline void grad ();

// ./stan/math/rev/core/grad.hpp:50
template <typename Vari>
static void grad (Vari *vi );

// ./stan/math/rev/core/reverse_pass_callback.hpp:37
template <typename F>
inline void reverse_pass_callback (F &&functor );

// ./stan/math/rev/core/var.hpp:22
template <typename Vari>
static void grad (Vari *vi );

// ./stan/math/rev/core/accumulate_adjoints.hpp:14
template <typename ...Pargs>
inline double *accumulate_adjoints (double *dest , const var &x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:17
template <typename VarVec,
    require_std_vector_vt<is_var,
    VarVec>*>
inline double *accumulate_adjoints (double *dest , VarVec &&x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:21
template <typename VecContainer,
    require_std_vector_st<is_var,
    VecContainer>*>
inline double *accumulate_adjoints (double *dest , VecContainer &&x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:28
template <typename EigT,
    require_eigen_vt<is_var,
    EigT>*>
inline double *accumulate_adjoints (double *dest , EigT &&x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:32
template <typename Arith,
    require_st_arithmetic<Arith>*>
inline double *accumulate_adjoints (double *dest , Arith &&x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:36

inline double *accumulate_adjoints (double *dest );

// ./stan/math/rev/core/accumulate_adjoints.hpp:50
template <typename ...Pargs>
inline double *accumulate_adjoints (double *dest , const var &x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:69
template <typename VarVec,
    require_std_vector_vt<is_var,
    VarVec>*,
    typename ...Pargs>
inline double *accumulate_adjoints (double *dest , VarVec &&x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:94
template <typename VecContainer,
    require_std_vector_st<is_var,
    VecContainer>*,
    require_std_vector_vt<is_container,
    VecContainer>*,
    typename ...Pargs>
inline double *accumulate_adjoints (double *dest , VecContainer &&x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:118
template <typename EigT,
    require_eigen_vt<is_var,
    EigT>*,
    typename ...Pargs>
inline double *accumulate_adjoints (double *dest , EigT &&x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:138
template <typename Arith,
    require_st_arithmetic<Arith>*,
    typename ...Pargs>
inline double *accumulate_adjoints (double *dest , Arith &&x , Pargs &&...args );

// ./stan/math/rev/core/accumulate_adjoints.hpp:148

inline double *accumulate_adjoints (double *dest );

// ./stan/math/rev/core/build_vari_array.hpp:20
template <int R,
    int C>
vari **build_vari_array (const Eigen ::Matrix <var , R , C >&x );

// ./stan/math/rev/core/count_vars.hpp:145
template <typename ...Pargs>
inline size_t count_vars (Pargs &&...args );

// ./stan/math/rev/core/callback_vari.hpp:40
template <typename T,
    typename F>
internal ::callback_vari <plain_type_t <T >,F >*make_callback_vari (T &&value , F &&functor );

// ./stan/math/rev/core/callback_vari.hpp:60
template <typename T,
    typename F>
var_value <plain_type_t <T >>make_callback_var (T &&value , F &&functor );

// ./stan/math/rev/core/deep_copy_vars.hpp:22
template <typename Arith,
    typename 
;

// ./stan/math/rev/core/deep_copy_vars.hpp:33

inline auto deep_copy_vars (const var &arg );

// ./stan/math/rev/core/deep_copy_vars.hpp:44
template <typename VarVec,
    require_std_vector_vt<is_var,
    VarVec>*>
inline auto deep_copy_vars (VarVec &&arg );

// ./stan/math/rev/core/deep_copy_vars.hpp:60
template <typename VecContainer,
    require_std_vector_st<is_var,
    VecContainer>*>
inline auto deep_copy_vars (VecContainer &&arg );

// ./stan/math/rev/core/deep_copy_vars.hpp:78
template <typename EigT,
    require_eigen_vt<is_var,
    EigT>*>
inline auto deep_copy_vars (EigT &&arg );

// ./stan/math/prim/err/invalid_argument.hpp:26
template <typename T>
inline void invalid_argument (const char *function , const char *name , const T &y , const char *msg1 , const char *msg2 );

// ./stan/math/prim/err/invalid_argument.hpp:48
template <typename T>
inline void invalid_argument (const char *function , const char *name , const T &y , const char *msg1 );

// ./stan/math/prim/core/init_threadpool_tbb.hpp:100

inline tbb ::task_arena &init_threadpool_tbb (int n_threads =0);

// ./stan/math/prim/core/operator_addition.hpp:39
template <typename U,
    typename V,
    require_all_stan_scalar_t<U,
    V>*>
inline complex_return_t <U ,V >operator +(const std ::complex <U >&x , const std ::complex <V >&y );

// ./stan/math/prim/core/operator_addition.hpp:54
template <typename U,
    typename V,
    require_all_stan_scalar_t<U,
    V>*>
inline complex_return_t <U ,V >operator +(const std ::complex <U >&x , const V &y );

// ./stan/math/prim/core/operator_addition.hpp:68
template <typename U,
    typename V,
    require_all_stan_scalar_t<U,
    V>*>
inline complex_return_t <U ,V >operator +(const U &x , const std ::complex <V >&y );

// ./stan/math/prim/core/operator_equal_equal.hpp:20
template <typename U,
    typename V,
    typename >
inline bool operator ==(const std ::complex <U >&x , const std ::complex <V >&y );

// ./stan/math/prim/core/operator_equal_equal.hpp:35
template <typename U,
    typename V,
    typename >
inline bool operator ==(const std ::complex <U >&x , const V &y );

// ./stan/math/prim/core/operator_equal_equal.hpp:51
template <typename U,
    typename V,
    typename >
inline bool operator ==(const U &x , const std ::complex <V >&y );

// ./stan/math/prim/core/operator_minus.hpp:17
template <typename U,
    require_autodiff_t<U>>
inline std ::complex <U >operator -(const std ::complex <U >&x );

// ./stan/math/prim/core/operator_not_equal.hpp:20
template <typename U,
    typename V,
    typename >
inline bool operator !=(const std ::complex <U >&x , const std ::complex <V >&y );

// ./stan/math/prim/core/operator_not_equal.hpp:36
template <typename U,
    typename V,
    typename >
inline bool operator !=(const std ::complex <U >&x , const V &y );

// ./stan/math/prim/core/operator_not_equal.hpp:52
template <typename U,
    typename V,
    typename >
inline bool operator !=(const U &x , const std ::complex <V >&y );

// ./stan/math/prim/core/operator_plus.hpp:17
template <typename U,
    require_autodiff_t<U>>
inline std ::complex <U >operator +(const std ::complex <U >&x );

// ./stan/math/prim/core/operator_subtraction.hpp:39
template <typename U,
    typename V,
    require_all_stan_scalar_t<U,
    V>*>
inline complex_return_t <U ,V >operator -(const std ::complex <U >&x , const std ::complex <V >&y );

// ./stan/math/prim/core/operator_subtraction.hpp:54
template <typename U,
    typename V,
    require_all_stan_scalar_t<U,
    V>*>
inline complex_return_t <U ,V >operator -(const std ::complex <U >&x , const V &y );

// ./stan/math/prim/core/operator_subtraction.hpp:68
template <typename U,
    typename V,
    require_all_stan_scalar_t<U,
    V>*>
inline complex_return_t <U ,V >operator -(const U &x , const std ::complex <V >&y );

// ./stan/math/rev/core/read_var.hpp:138
template <typename EigVar,
    typename EigVari,
    typename EigDbl>
inline void read_vi_val_adj (const EigVar &VarMat , EigVari &VariMat , EigDbl &ValMat , EigDbl &AdjMat );

// ./stan/math/rev/core/read_var.hpp:157
template <typename EigRev,
    typename EigDbl>
inline void read_val_adj (const EigRev &VarMat , EigDbl &ValMat , EigDbl &AdjMat );

// ./stan/math/rev/core/read_var.hpp:174
template <typename EigVar,
    typename EigVari,
    typename EigDbl>
inline void read_vi_val (const EigVar &VarMat , EigVari &VariMat , EigDbl &ValMat );

// ./stan/math/rev/core/read_var.hpp:192
template <typename EigVar,
    typename EigVari,
    typename EigDbl>
inline void read_vi_adj (const EigVar &VarMat , EigVari &VariMat , EigDbl &AdjMat );

// ./stan/math/rev/core/recover_memory_nested.hpp:23

static inline void recover_memory_nested ();

// ./stan/math/rev/core/set_zero_all_adjoints_nested.hpp:20

static inline void set_zero_all_adjoints_nested ();

// ./stan/math/rev/core/start_nested.hpp:16

static inline void start_nested ();

// ./stan/math/prim/fun/dims.hpp:19
template <typename T,
    require_stan_scalar_t<T>*>
inline void dims (const T &x , std ::vector <int >&result );

// ./stan/math/prim/fun/dims.hpp:29
template <typename T,
    require_matrix_t<T>*>
inline void dims (const T &x , std ::vector <int >&result );

// ./stan/math/prim/fun/dims.hpp:44
template <typename T,
    typename Alloc>
inline void dims (const std ::vector <T , Alloc >&x , std ::vector <int >&result );

// ./stan/math/prim/fun/dims.hpp:57
template <typename T>
inline std ::vector <int >dims (const T &x );

// ./stan/math/prim/err/check_matching_dims.hpp:26
template <typename T1,
    typename T2,
    require_all_not_matrix_t<T1,
    T2>*>
inline void check_matching_dims (const char *function , const char *name1 , const T1 &y1 , const char *name2 , const T2 &y2 );

// ./stan/math/prim/err/check_matching_dims.hpp:76
template <typename T1,
    typename T2,
    require_any_t<conjunction<is_matrix<T1>,
    is_matrix<T2>>,
    conjunction<is_prim_or_rev_kernel_expression<T1>,
    is_prim_or_rev_kernel_expression<T2>>>*>
inline void check_matching_dims (const char *function , const char *name1 , const T1 &y1 , const char *name2 , const T2 &y2 );

// ./stan/math/prim/err/check_matching_dims.hpp:112
template <typename T1,
    typename T2,
    require_any_matrix_t<T1,
    T2>*>
inline void check_matching_dims (const char *function , const char *name1 , const T1 &y1 , const char *name2 , const T2 &y2 );

// ./stan/math/prim/err/check_matching_dims.hpp:137
template <bool check_compile,
    typename Mat1,
    typename Mat2,
    typename >
inline void check_matching_dims (const char *function , const char *name1 , const Mat1 &y1 , const char *name2 , const Mat2 &y2 );

// ./stan/math/prim/fun/as_array_or_scalar.hpp:18
template <typename T,
    require_stan_scalar_t<T>*>
inline T as_array_or_scalar (T &&v );

// ./stan/math/prim/fun/as_array_or_scalar.hpp:30
template <typename T,
    require_stan_scalar_t<T>*>
inline T &as_array_or_scalar (T &v );

// ./stan/math/prim/fun/as_array_or_scalar.hpp:42
template <typename T,
    require_eigen_array_t<T>*>
inline T as_array_or_scalar (T &&v );

// ./stan/math/prim/fun/as_array_or_scalar.hpp:54
template <typename T,
    typename >
inline auto as_array_or_scalar (T &&v );

// ./stan/math/prim/fun/as_array_or_scalar.hpp:67
template <typename T,
    require_std_vector_t<T>*
;

// ./stan/math/prim/fun/as_array_or_scalar.hpp:83
template <typename T,
    require_std_vector_vt<is_std_vector,
    T>*
;

// ./stan/math/prim/fun/inv.hpp:34
template <typename T,
    require_not_container_st<std::is_arithmetic,
    T>*>
inline auto inv (const T &x );

// ./stan/math/prim/fun/inv.hpp:49
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto inv (const Container &x );

// ./stan/math/prim/fun/constants.hpp:20

static constexpr double e ();

// ./stan/math/prim/fun/constants.hpp:27

static constexpr double egamma ();

// ./stan/math/prim/fun/constants.hpp:36

static constexpr double pi ();

// ./stan/math/prim/fun/constants.hpp:199

static constexpr inline double positive_infinity ();

// ./stan/math/prim/fun/constants.hpp:206

static constexpr inline double negative_infinity ();

// ./stan/math/prim/fun/constants.hpp:213

static constexpr inline double not_a_number ();

// ./stan/math/prim/fun/constants.hpp:221

static constexpr inline double machine_precision ();

// ./stan/math/prim/fun/constants.hpp:228

static constexpr inline double log10 ();

// ./stan/math/prim/fun/constants.hpp:235

static constexpr inline double sqrt2 ();

// ./stan/math/rev/core/operator_addition.hpp:53

inline var operator +(const var &a , const var &b );

// ./stan/math/rev/core/operator_addition.hpp:73
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline var operator +(const var &a , Arith b );

// ./stan/math/rev/core/operator_addition.hpp:95
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline var operator +(Arith a , const var &b );

// ./stan/math/rev/core/operator_addition.hpp:109
template <typename VarMat1,
    typename VarMat2,
    require_all_rev_matrix_t<VarMat1,
    VarMat2>*>
inline auto add (VarMat1 &&a , VarMat2 &&b );

// ./stan/math/rev/core/operator_addition.hpp:139
template <typename Arith,
    typename VarMat,
    require_st_arithmetic<Arith>*>
inline auto add (VarMat &&a , const Arith &b );

// ./stan/math/rev/core/operator_addition.hpp:165
template <typename Arith,
    typename VarMat,
    require_st_arithmetic<Arith>*>
inline auto add (const Arith &a , VarMat &&b );

// ./stan/math/rev/core/operator_addition.hpp:181
template <typename Var,
    typename EigMat,
    require_var_vt<std::is_arithmetic,
    Var>*>
inline auto add (const Var &a , const EigMat &b );

// ./stan/math/rev/core/operator_addition.hpp:200
template <typename EigMat,
    typename Var,
    require_eigen_vt<std::is_arithmetic,
    EigMat>*>
inline auto add (const EigMat &a , const Var &b );

// ./stan/math/rev/core/operator_addition.hpp:217
template <typename Var,
    typename VarMat,
    require_var_vt<std::is_arithmetic,
    Var>*>
inline auto add (const Var &a , VarMat &&b );

// ./stan/math/rev/core/operator_addition.hpp:246
template <typename Var,
    typename VarMat,
    require_var_vt<std::is_arithmetic,
    Var>*>
inline auto add (VarMat &&a , const Var &b );

// ./stan/math/rev/core/operator_addition.hpp:253
template <typename T1,
    typename T2,
    require_any_var_vt<std::is_arithmetic,
    T1,
    T2>*>
inline auto add (const T1 &a , const T2 &b );

// ./stan/math/rev/core/operator_addition.hpp:260
template <typename T1,
    typename T2,
    require_all_var_vt<std::is_arithmetic,
    T1,
    T2>*>
inline auto add (const T1 &a , const T2 &b );

// ./stan/math/rev/core/operator_addition.hpp:275
template <typename VarMat1,
    typename VarMat2,
    require_any_var_matrix_t<VarMat1,
    VarMat2>*>
inline auto operator +(VarMat1 &&a , VarMat2 &&b );

// ./stan/math/prim/fun/is_any_nan.hpp:20
template <typename T>
inline bool is_any_nan (const T &x );

// ./stan/math/prim/fun/is_any_nan.hpp:34
template <typename T,
    typename ...Ts>
inline bool is_any_nan (const T &x , const Ts &...xs );

// ./stan/math/prim/fun/value_of.hpp:20
template <typename T,
    require_st_arithmetic<T>*>
inline T value_of (T &&x );

// ./stan/math/prim/fun/value_of.hpp:48
template <typename T,
    require_std_vector_t<T>*>
inline auto value_of (const T &x );

// ./stan/math/prim/fun/value_of.hpp:70
template <typename EigMat,
    require_eigen_dense_base_t<EigMat>*>
inline auto value_of (EigMat &&M );

// ./stan/math/prim/fun/value_of.hpp:80
template <typename EigMat,
    require_eigen_sparse_base_t<EigMat>*>
inline auto value_of (EigMat &&M );

// ./stan/math/prim/fun/value_of.hpp:96
template <typename EigMat,
    require_eigen_sparse_base_t<EigMat>*>
inline auto value_of (EigMat &&M );

// ./stan/math/rev/core/operator_subtraction.hpp:56

inline var operator -(const var &a , const var &b );

// ./stan/math/rev/core/operator_subtraction.hpp:77
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline var operator -(const var &a , Arith b );

// ./stan/math/rev/core/operator_subtraction.hpp:100
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline var operator -(Arith a , const var &b );

// ./stan/math/rev/core/operator_subtraction.hpp:116
template <typename VarMat1,
    typename VarMat2,
    require_all_rev_matrix_t<VarMat1,
    VarMat2>*>
inline auto subtract (const VarMat1 &a , const VarMat2 &b );

// ./stan/math/rev/core/operator_subtraction.hpp:146
template <typename Arith,
    typename VarMat,
    require_st_arithmetic<Arith>*>
inline auto subtract (const VarMat &a , const Arith &b );

// ./stan/math/rev/core/operator_subtraction.hpp:172
template <typename Arith,
    typename VarMat,
    require_st_arithmetic<Arith>*>
inline auto subtract (const Arith &a , const VarMat &b );

// ./stan/math/rev/core/operator_subtraction.hpp:198
template <typename Var,
    typename EigMat,
    require_var_vt<std::is_arithmetic,
    Var>*>
inline auto subtract (const Var &a , const EigMat &b );

// ./stan/math/rev/core/operator_subtraction.hpp:217
template <typename EigMat,
    typename Var,
    require_eigen_vt<std::is_arithmetic,
    EigMat>*>
inline auto subtract (const EigMat &a , const Var &b );

// ./stan/math/rev/core/operator_subtraction.hpp:237
template <typename Var,
    typename VarMat,
    require_var_vt<std::is_arithmetic,
    Var>*>
inline auto subtract (const Var &a , const VarMat &b );

// ./stan/math/rev/core/operator_subtraction.hpp:266
template <typename Var,
    typename VarMat,
    require_rev_matrix_t<VarMat>*>
inline auto subtract (const VarMat &a , const Var &b );

// ./stan/math/rev/core/operator_subtraction.hpp:285
template <typename T1,
    typename T2,
    require_any_var_vt<std::is_arithmetic,
    T1,
    T2>*>
inline auto subtract (const T1 &a , const T2 &b );

// ./stan/math/rev/core/operator_subtraction.hpp:292
template <typename T1,
    typename T2,
    require_all_var_vt<std::is_arithmetic,
    T1,
    T2>*>
inline auto subtract (const T1 &a , const T2 &b );

// ./stan/math/rev/core/operator_subtraction.hpp:307
template <typename VarMat1,
    typename VarMat2,
    require_any_var_matrix_t<VarMat1,
    VarMat2>*>
inline auto operator -(const VarMat1 &a , const VarMat2 &b );

// ./stan/math/prim/fun/size.hpp:16
template <typename T,
    require_stan_scalar_t<T>*>
inline int64_t size (const T &/*x*/);

// ./stan/math/prim/fun/size.hpp:27
template <typename T,
    require_container_t<T>*>
inline int64_t size (const T &m );

// ./stan/math/prim/fun/size.hpp:32
template <typename T,
    require_var_matrix_t<T>*>
inline int64_t size (const T &m );

// ./stan/math/prim/fun/is_inf.hpp:18

inline bool is_inf (double x );

// ./stan/math/prim/fun/isinf.hpp:22
template <typename T,
    typename >
inline bool isinf (const T &v );

// ./stan/math/prim/fun/isnan.hpp:20
template <typename T,
    typename >
inline bool isnan (const T &x );

// ./stan/math/rev/core/operator_multiplication.hpp:76

inline var operator *(const var &a , const var &b );

// ./stan/math/rev/core/operator_multiplication.hpp:92
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline var operator *(const var &a , Arith b );

// ./stan/math/rev/core/operator_multiplication.hpp:112
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline var operator *(Arith a , const var &b );

// ./stan/math/rev/core/operator_multiplication.hpp:127

inline std ::complex <stan ::math ::var >operator *(const std ::complex <stan ::math ::var >&x , const std ::complex <stan ::math ::var >&y );

// ./stan/math/rev/core/operator_multiplication.hpp:141

inline std ::complex <stan ::math ::var >operator *(const std ::complex <double >&x , const std ::complex <stan ::math ::var >&y );

// ./stan/math/rev/core/operator_multiplication.hpp:154

inline std ::complex <stan ::math ::var >operator *(const std ::complex <stan ::math ::var >&x , const std ::complex <double >&y );

// ./stan/math/rev/fun/to_arena.hpp:25
template <typename T,
    require_not_same_t<T,
    arena_t<T>>*>
inline arena_t <T >to_arena (T &&a );

// ./stan/math/rev/fun/to_arena.hpp:46
template <typename T,
    require_same_t<T,
    arena_t<T>>*>
inline std ::remove_reference_t <T >to_arena (T &&a );

// ./stan/math/rev/fun/to_arena.hpp:65
template <typename T,
    require_eigen_t<T>*
;

// ./stan/math/rev/fun/to_arena.hpp:81
template <typename T>
inline std ::vector <T ,arena_allocator <T >>to_arena (const std ::vector <T , arena_allocator <T >>&a );

// ./stan/math/rev/fun/to_arena.hpp:144
template <bool Condition,
    typename T,
    std::enable_if_t<!Condition>*>
inline T to_arena_if (T &&a );

// ./stan/math/rev/fun/to_arena.hpp:149
template <bool Condition,
    typename T,
    std::enable_if_t<Condition>*>
inline arena_t <T >to_arena_if (const T &a );

// ./stan/math/rev/fun/value_of.hpp:23
template <typename T>
inline auto &value_of (const var_value <T >&v );

// ./stan/math/rev/core/operator_division.hpp:61

inline var operator /(const var &dividend , const var &divisor );

// ./stan/math/rev/core/operator_division.hpp:83
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline var operator /(const var &dividend , Arith divisor );

// ./stan/math/rev/core/operator_division.hpp:105
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline var operator /(Arith dividend , const var &divisor );

// ./stan/math/rev/core/operator_division.hpp:122
template <typename Scalar,
    typename Mat,
    require_matrix_t<Mat>*>
inline auto divide (const Mat &m , Scalar c );

// ./stan/math/rev/core/operator_division.hpp:166
template <typename Scalar,
    typename Mat,
    require_matrix_t<Mat>*>
inline auto divide (Scalar c , const Mat &m );

// ./stan/math/rev/core/operator_division.hpp:212
template <typename Mat1,
    typename Mat2,
    require_all_matrix_st<is_var_or_arithmetic,
    Mat1,
    Mat2>*>
inline auto divide (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/rev/core/operator_division.hpp:254
template <typename T1,
    typename T2,
    require_any_var_matrix_t<T1,
    T2>*>
inline auto operator /(const T1 &dividend , const T2 &divisor );

// ./stan/math/rev/core/operator_division.hpp:259

inline std ::complex <var >operator /(const std ::complex <var >&x1 , const std ::complex <var >&x2 );

// ./stan/math/rev/core/operator_equal.hpp:27

inline bool operator ==(const var &a , const var &b );

// ./stan/math/rev/core/operator_equal.hpp:41
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator ==(const var &a , Arith b );

// ./stan/math/rev/core/operator_equal.hpp:55
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator ==(Arith a , const var &b );

// ./stan/math/rev/core/operator_equal.hpp:69

inline bool operator ==(const var &x , const std ::complex <var >&z );

// ./stan/math/rev/core/operator_equal.hpp:82

inline bool operator ==(const std ::complex <var >&z , const var &y );

// ./stan/math/rev/core/operator_greater_than.hpp:26

inline bool operator >(const var &a , const var &b );

// ./stan/math/rev/core/operator_greater_than.hpp:37
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator >(const var &a , Arith b );

// ./stan/math/rev/core/operator_greater_than.hpp:51
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator >(Arith a , const var &b );

// ./stan/math/rev/core/operator_greater_than_or_equal.hpp:28

inline bool operator >=(const var &a , const var &b );

// ./stan/math/rev/core/operator_greater_than_or_equal.hpp:42
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator >=(const var &a , Arith b );

// ./stan/math/rev/core/operator_greater_than_or_equal.hpp:57
template <typename Arith,
    typename Var,
    require_arithmetic_t<Arith>*>
inline bool operator >=(Arith a , const var &b );

// ./stan/math/rev/core/operator_less_than.hpp:25

inline bool operator <(const var &a , const var &b );

// ./stan/math/rev/core/operator_less_than.hpp:36
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator <(const var &a , Arith b );

// ./stan/math/rev/core/operator_less_than.hpp:50
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator <(Arith a , const var &b );

// ./stan/math/rev/core/operator_less_than_or_equal.hpp:27

inline bool operator <=(const var &a , const var &b );

// ./stan/math/rev/core/operator_less_than_or_equal.hpp:41
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator <=(const var &a , Arith b );

// ./stan/math/rev/core/operator_less_than_or_equal.hpp:56
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator <=(Arith a , const var &b );

// ./stan/math/rev/core/operator_logical_and.hpp:18

inline bool operator &&(const var &x , const var &y );

// ./stan/math/rev/core/operator_logical_and.hpp:33
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator &&(const var &x , Arith y );

// ./stan/math/rev/core/operator_logical_and.hpp:49
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator &&(Arith x , const var &y );

// ./stan/math/rev/core/operator_logical_or.hpp:18

inline bool operator ||(const var &x , const var &y );

// ./stan/math/rev/core/operator_logical_or.hpp:32
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator ||(const var &x , Arith y );

// ./stan/math/rev/core/operator_logical_or.hpp:47
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator ||(Arith x , const var &y );

// ./stan/math/rev/core/operator_not_equal.hpp:30

inline bool operator !=(const var &a , const var &b );

// ./stan/math/rev/core/operator_not_equal.hpp:44
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator !=(const var &a , Arith b );

// ./stan/math/rev/core/operator_not_equal.hpp:59
template <typename Arith,
    require_arithmetic_t<Arith>*>
inline bool operator !=(Arith a , const var &b );

// ./stan/math/rev/core/operator_not_equal.hpp:73

inline bool operator !=(const var &x , const std ::complex <var >&z );

// ./stan/math/rev/core/operator_not_equal.hpp:86

inline bool operator !=(const std ::complex <var >&z , const var &y );

// ./stan/math/rev/core/operator_unary_decrement.hpp:26

inline var &operator --(var &a );

// ./stan/math/rev/core/operator_unary_decrement.hpp:42

inline var operator --(var &a , int /*dummy*/);

// ./stan/math/rev/core/operator_unary_increment.hpp:22

inline var &operator ++(var &a );

// ./stan/math/rev/core/operator_unary_increment.hpp:38

inline var operator ++(var &a , int /*dummy*/);

// ./stan/math/rev/core/operator_unary_negative.hpp:37

inline var operator -(const var &a );

// ./stan/math/rev/core/operator_unary_negative.hpp:49
template <typename T,
    require_var_matrix_t<T>*>
inline auto operator -(const T &a );

// ./stan/math/rev/core/operator_unary_not.hpp:17

inline bool operator !(const var &x );

// ./stan/math/rev/core/operator_unary_plus.hpp:44

inline var operator +(const var &a );

// ./stan/math/rev/fun/dims.hpp:20
template <typename T>
inline void dims (const var_value <T >&x , std ::vector <int >&result );

// ./stan/math/rev/fun/dims.hpp:33
template <typename T>
inline void dims (const vari_value <T >&x , std ::vector <int >&result );

// ./stan/math/rev/core/precomputed_gradients.hpp:211
template <typename Arith,
    typename VecVar,
    typename VecArith,
    typename ...ContainerOperands,
    typename ...ContainerGradients>
inline var precomputed_gradients (Arith value , const VecVar &operands , const VecArith &gradients , const std ::tuple <ContainerOperands ...>&container_operands =std ::tuple <>(),const std ::tuple <ContainerGradients ...>&container_gradients =std ::tuple <>());

// ./stan/math/rev/core/print_stack.hpp:20

inline void print_stack (std ::ostream &o );

// ./stan/math/prim/err/throw_domain_error.hpp:25
template <typename T>
inline void throw_domain_error (const char *function , const char *name , const T &y , const char *msg1 , const char *msg2 );

// ./stan/math/prim/err/throw_domain_error.hpp:48
template <typename T>
inline void throw_domain_error (const char *function , const char *name , const T &y , const char *msg1 );

// ./stan/math/prim/fun/value_of_rec.hpp:28
template <typename T,
    typename >
inline double value_of_rec (const T x );

// ./stan/math/prim/fun/value_of_rec.hpp:44

inline double value_of_rec (double x );

// ./stan/math/prim/fun/value_of_rec.hpp:53
template <typename T>
inline std ::complex <double >value_of_rec (const std ::complex <T >&x );

// ./stan/math/prim/fun/value_of_rec.hpp:68
template <typename T,
    require_not_same_t<double ,
    T>*>
inline std ::vector <double >value_of_rec (const std ::vector <T >&x );

// ./stan/math/prim/fun/value_of_rec.hpp:89
template <typename T,
    require_std_vector_t<T>*>
inline T value_of_rec (T &&x );

// ./stan/math/prim/fun/value_of_rec.hpp:105
template <typename T,
    typename >
inline auto value_of_rec (T &&M );

// ./stan/math/prim/fun/value_of_rec.hpp:127
template <typename T,
    typename >
inline T value_of_rec (T &&x );

// ./stan/math/prim/err/elementwise_check.hpp:114
template <typename F,
    typename T,
    typename ...Indexings,
    require_stan_scalar_t<T>*>
inline void elementwise_check (const F &is_good , const char *function , const char *name , const T &x , const char *must_be , const Indexings &...indexings );

// ./stan/math/prim/err/elementwise_check.hpp:145
template <typename F,
    typename T,
    typename ...Indexings,
    require_eigen_t<T>*>
inline void elementwise_check (const F &is_good , const char *function , const char *name , const T &x , const char *must_be , const Indexings &...indexings );

// ./stan/math/prim/err/elementwise_check.hpp:197
template <typename F,
    typename T,
    typename ...Indexings,
    require_eigen_t<T>*>
inline void elementwise_check (const F &is_good , const char *function , const char *name , const T &x , const char *must_be , const Indexings &...indexings );

// ./stan/math/prim/err/elementwise_check.hpp:240
template <typename F,
    typename T,
    typename ...Indexings,
    require_eigen_t<T>*>
inline void elementwise_check (const F &is_good , const char *function , const char *name , const T &x , const char *must_be , const Indexings &...indexings );

// ./stan/math/prim/err/elementwise_check.hpp:283
template <typename F,
    typename T,
    typename ...Indexings,
    require_std_vector_t<T>*>
inline void elementwise_check (const F &is_good , const char *function , const char *name , const T &x , const char *must_be , const Indexings &...indexings );

// ./stan/math/prim/err/elementwise_check.hpp:312
template <typename F,
    typename T,
    typename ...Indexings,
    require_var_matrix_t<T>*>
inline void elementwise_check (const F &is_good , const char *function , const char *name , const T &x , const char *must_be , const Indexings &...indexings );

// ./stan/math/prim/err/elementwise_check.hpp:332
template <typename F,
    typename T>
inline bool elementwise_is (const F &is_good , const T &x );

// ./stan/math/prim/err/check_not_nan.hpp:25
template <typename T_y>
inline void check_not_nan (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/fun/fabs.hpp:13
template <typename T,
    require_arithmetic_t<T>*>
inline auto fabs (T x );

// ./stan/math/prim/fun/fabs.hpp:18
template <typename T,
    require_complex_t<T>*>
inline auto fabs (T x );

// ./stan/math/prim/fun/fabs.hpp:45
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto fabs (const Container &x );

// ./stan/math/prim/fun/fabs.hpp:63
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto fabs (const Container &x );

// ./stan/math/prim/fun/floor.hpp:36
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto floor (const Container &x );

// ./stan/math/prim/fun/floor.hpp:53
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto floor (const Container &x );

// ./stan/math/prim/fun/is_nonpositive_integer.hpp:17
template <typename T>
inline bool is_nonpositive_integer (T x );

// ./stan/math/prim/err/check_2F1_converges.hpp:34
template <typename T_a1,
    typename T_a2,
    typename T_b1,
    typename T_z>
inline void check_2F1_converges (const char *function , const T_a1 &a1 , const T_a2 &a2 , const T_b1 &b1 , const T_z &z );

// ./stan/math/prim/err/check_3F2_converges.hpp:38
template <typename T_a1,
    typename T_a2,
    typename T_a3,
    typename T_b1,
    typename T_b2,
    typename T_z>
inline void check_3F2_converges (const char *function , const T_a1 &a1 , const T_a2 &a2 , const T_a3 &a3 , const T_b1 &b1 , const T_b2 &b2 , const T_z &z );

// ./stan/math/prim/err/throw_domain_error_vec.hpp:30
template <typename T>
inline void throw_domain_error_vec (const char *function , const char *name , const T &y , size_t i , const char *msg1 , const char *msg2 );

// ./stan/math/prim/err/throw_domain_error_vec.hpp:56
template <typename T>
inline void throw_domain_error_vec (const char *function , const char *name , const T &y , size_t i , const char *msg );

// ./stan/math/prim/fun/max_size.hpp:19
template <typename T1,
    typename ...Ts>
inline int64_t max_size (const T1 &x1 , const Ts &...xs );

// ./stan/math/prim/fun/size_zero.hpp:18
template <typename T>
inline bool size_zero (const T &x );

// ./stan/math/prim/fun/size_zero.hpp:31
template <typename T,
    typename ...Ts>
inline bool size_zero (const T &x , const Ts &...xs );

// ./stan/math/prim/err/check_bounded.hpp:74
template <typename T_y,
    typename T_low,
    typename T_high>
inline void check_bounded (const char *function , const char *name , const T_y &y , const T_low &low , const T_high &high );

// ./stan/math/prim/err/check_positive.hpp:26
template <typename T_y>
inline void check_positive (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/err/check_positive.hpp:42

inline void check_positive (const char *function , const char *name , const char *expr , int size );

// ./stan/math/prim/err/check_size_match.hpp:23
template <typename T_size1,
    typename T_size2>
inline void check_size_match (const char *function , const char *name_i , T_size1 i , const char *name_j , T_size2 j );

// ./stan/math/prim/err/check_size_match.hpp:49
template <typename T_size1,
    typename T_size2>
inline void check_size_match (const char *function , const char *expr_i , const char *name_i , T_size1 i , const char *expr_j , const char *name_j , T_size2 j );

// ./stan/math/prim/err/check_matching_sizes.hpp:23
template <typename T_y1,
    typename T_y2>
inline void check_matching_sizes (const char *function , const char *name1 , const T_y1 &y1 , const char *name2 , const T_y2 &y2 );

// ./stan/math/prim/err/throw_domain_error_mat.hpp:31
template <typename T,
    require_eigen_t<T>*>
inline void throw_domain_error_mat (const char *function , const char *name , const T &y , size_t i , size_t j , const char *msg1 , const char *msg2 );

// ./stan/math/prim/err/throw_domain_error_mat.hpp:65
template <typename T>
inline void throw_domain_error_mat (const char *function , const char *name , const T &y , size_t i , size_t j , const char *msg );

// ./stan/math/prim/err/check_less_or_equal.hpp:34
template <typename T_y,
    typename T_high,
    require_all_stan_scalar_t<T_y,
    T_high>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less_or_equal.hpp:69
template <typename T_y,
    typename T_high,
    require_stan_scalar_t<T_y>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less_or_equal.hpp:110
template <typename T_y,
    typename T_high,
    require_stan_scalar_t<T_y>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less_or_equal.hpp:151
template <typename T_y,
    typename T_high,
    require_vector_t<T_y>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less_or_equal.hpp:191
template <typename T_y,
    typename T_high,
    require_dense_dynamic_t<T_y>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less_or_equal.hpp:235
template <typename T_y,
    typename T_high,
    require_all_vector_t<T_y,
    T_high>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less_or_equal.hpp:280
template <typename T_y,
    typename T_high,
    require_all_dense_dynamic_t<T_y,
    T_high>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less_or_equal.hpp:321
template <typename T_y,
    typename T_high,
    require_std_vector_vt<is_container_or_var_matrix,
    T_y>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less_or_equal.hpp:349
template <typename T_y,
    typename T_high,
    require_not_std_vector_t<T_y>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less_or_equal.hpp:379
template <typename T_y,
    typename T_high,
    require_any_std_vector_vt<is_container_or_var_matrix,
    T_y,
    T_high>*>
inline void check_less_or_equal (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_lower_triangular.hpp:26
template <typename T_y,
    require_eigen_t<T_y>*>
inline void check_lower_triangular (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/err/check_cholesky_factor.hpp:32
template <typename Mat,
    require_matrix_t<Mat>*>
inline void check_cholesky_factor (const char *function , const char *name , const Mat &y );

// ./stan/math/prim/err/check_cholesky_factor.hpp:60
template <typename StdVec,
    require_std_vector_t<StdVec>*>
void check_cholesky_factor (const char *function , const char *name , const StdVec &y );

// ./stan/math/prim/err/check_square.hpp:20
template <typename T_y,
    require_any_t<is_matrix<T_y>,
    is_prim_or_rev_kernel_expression<T_y>>*>
inline void check_square (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/err/check_nonzero_size.hpp:21
template <typename T_y>
inline void check_nonzero_size (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/fun/num_elements.hpp:19
template <typename T,
    require_stan_scalar_t<T>*>
inline int64_t num_elements (const T &x );

// ./stan/math/prim/fun/num_elements.hpp:32
template <typename T,
    require_matrix_t<T>*>
inline int64_t num_elements (const T &m );

// ./stan/math/prim/fun/num_elements.hpp:46
template <typename T>
inline int64_t num_elements (const std ::vector <T >&v );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:35
template <typename T1,
    typename T2,
    typename F,
    require_all_stan_scalar_t<T1,
    T2>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:54
template <typename T1,
    typename T2,
    typename F,
    require_all_eigen_t<T1,
    T2>*>
inline auto apply_scalar_binary (T1 &&x , T2 &&y , F &&f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:77
template <typename T1,
    typename T2,
    typename F,
    require_eigen_vector_vt<is_stan_scalar,
    T1>*>
inline auto apply_scalar_binary (T1 &&x , T2 &&y , F &&f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:104
template <typename T1,
    typename T2,
    typename F,
    require_std_vector_vt<std::is_integral,
    T1>*>
inline auto apply_scalar_binary (T1 &&x , T2 &&y , F &&f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:131
template <typename T1,
    typename T2,
    typename F,
    require_eigen_matrix_dynamic_vt<is_stan_scalar,
    T1>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:164
template <typename T1,
    typename T2,
    typename F,
    require_std_vector_vt<is_std_vector,
    T1>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:200
template <typename T1,
    typename T2,
    typename F,
    require_eigen_t<T1>*>
inline auto apply_scalar_binary (T1 &&x , T2 &&y , F &&f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:226
template <typename T1,
    typename T2,
    typename F,
    require_stan_scalar_t<T1>*>
inline auto apply_scalar_binary (T1 &&x , T2 &&y , F &&f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:254
template <typename T1,
    typename T2,
    typename F,
    require_all_std_vector_vt<is_stan_scalar,
    T1,
    T2>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:285
template <typename T1,
    typename T2,
    typename F,
    require_std_vector_vt<is_stan_scalar,
    T1>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:315
template <typename T1,
    typename T2,
    typename F,
    require_stan_scalar_t<T1>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:341
template <typename T1,
    typename T2,
    typename F,
    require_all_std_vector_vt<is_container_or_var_matrix,
    T1,
    T2>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:369
template <typename T1,
    typename T2,
    typename F,
    require_std_vector_vt<is_container_or_var_matrix,
    T1>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/prim/functor/apply_scalar_binary.hpp:396
template <typename T1,
    typename T2,
    typename F,
    require_stan_scalar_t<T1>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/prim/fun/hypot.hpp:23
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline double hypot (T1 x , T2 y );

// ./stan/math/prim/fun/hypot.hpp:39
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto hypot (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/abs.hpp:24
template <typename T,
    require_arithmetic_t<T>*>
inline T abs (T x );

// ./stan/math/prim/fun/abs.hpp:37
template <typename T,
    require_complex_t<T>*>
inline auto abs (T x );

// ./stan/math/prim/fun/abs.hpp:65
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto abs (const Container &x );

// ./stan/math/prim/fun/abs.hpp:81
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto abs (const Container &x );

// ./stan/math/prim/err/check_unit_vector.hpp:35
template <typename Vec,
    require_vector_t<Vec>*>
void check_unit_vector (const char *function , const char *name , const Vec &theta );

// ./stan/math/prim/err/check_unit_vector.hpp:66
template <typename StdVec,
    require_std_vector_t<StdVec>*>
void check_unit_vector (const char *function , const char *name , const StdVec &theta );

// ./stan/math/prim/err/check_cholesky_factor_corr.hpp:36
template <typename Mat,
    require_matrix_t<Mat>*>
void check_cholesky_factor_corr (const char *function , const char *name , const Mat &y );

// ./stan/math/prim/err/check_cholesky_factor_corr.hpp:68
template <typename StdVec,
    require_std_vector_t<StdVec>*>
void check_cholesky_factor_corr (const char *function , const char *name , const StdVec &y );

// ./stan/math/prim/err/out_of_range.hpp:27

inline void out_of_range (const char *function , int max , int index , const char *msg1 ="", const char *msg2 ="");

// ./stan/math/prim/err/check_column_index.hpp:27
template <typename T_y,
    require_any_t<is_matrix<T_y>,
    is_prim_or_rev_kernel_expression<T_y>>*>
inline void check_column_index (const char *function , const char *name , const T_y &y , size_t i );

// ./stan/math/prim/err/check_consistent_size.hpp:24
template <typename T>
inline void check_consistent_size (const char *function , const char *name , const T &x , size_t expected_size );

// ./stan/math/prim/err/check_consistent_sizes.hpp:15

inline void check_consistent_sizes (const char *);

// ./stan/math/prim/err/check_consistent_sizes.hpp:21
template <typename T1>
inline void check_consistent_sizes (const char *, const char *, const T1 &);

// ./stan/math/prim/err/check_consistent_sizes.hpp:45
template <typename T1,
    typename T2,
    typename ...Ts>
inline void check_consistent_sizes (const char *function , const char *name1 , const T1 &x1 , const char *name2 , const T2 &x2 , const Ts &...names_and_xs );

// ./stan/math/prim/err/check_consistent_sizes_mvt.hpp:15

inline void check_consistent_sizes_mvt (const char *);

// ./stan/math/prim/err/check_consistent_sizes_mvt.hpp:21
template <typename T1>
inline void check_consistent_sizes_mvt (const char *, const char *, const T1 &);

// ./stan/math/prim/err/check_consistent_sizes_mvt.hpp:45
template <typename T1,
    typename T2,
    typename ...Ts>
inline void check_consistent_sizes_mvt (const char *function , const char *name1 , const T1 &x1 , const char *name2 , const T2 &x2 , const Ts &...names_and_xs );

// ./stan/math/prim/err/check_symmetric.hpp:31
template <typename EigMat,
    require_matrix_t<EigMat>*>
inline void check_symmetric (const char *function , const char *name , const EigMat &y );

// ./stan/math/prim/err/check_pos_definite.hpp:33
template <typename EigMat,
    require_matrix_t<EigMat>*>
inline void check_pos_definite (const char *function , const char *name , const EigMat &y );

// ./stan/math/prim/err/check_pos_definite.hpp:62
template <typename Derived>
inline void check_pos_definite (const char *function , const char *name , const Eigen ::LDLT <Derived >&cholesky );

// ./stan/math/prim/err/check_pos_definite.hpp:82
template <typename Derived>
inline void check_pos_definite (const char *function , const char *name , const Eigen ::LLT <Derived >&cholesky );

// ./stan/math/prim/err/check_corr_matrix.hpp:33
template <typename Mat,
    require_matrix_t<Mat>*>
inline void check_corr_matrix (const char *function , const char *name , const Mat &y );

// ./stan/math/prim/err/check_corr_matrix.hpp:75
template <typename StdVec,
    require_std_vector_t<StdVec>*>
void check_corr_matrix (const char *function , const char *name , const StdVec &y );

// ./stan/math/prim/err/check_cov_matrix.hpp:28
template <typename Mat,
    require_matrix_t<Mat>*>
inline void check_cov_matrix (const char *function , const char *name , const Mat &y );

// ./stan/math/prim/err/check_cov_matrix.hpp:51
template <typename StdVec,
    require_std_vector_t<StdVec>*>
void check_cov_matrix (const char *function , const char *name , const StdVec &y );

// ./stan/math/prim/err/check_finite.hpp:27
template <typename T_y>
inline void check_finite (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/err/domain_error.hpp:13
template <typename T>
inline void domain_error (const char *function , const char *name , const T &y , const char *msg1 , const char *msg2 );

// ./stan/math/prim/err/domain_error.hpp:22
template <typename T>
inline void domain_error (const char *function , const char *name , const T &y , const char *msg1 );

// ./stan/math/prim/err/check_flag_sundials.hpp:28

inline std ::array <std ::string ,2>cvodes_flag_msg (int flag );

// ./stan/math/prim/err/check_flag_sundials.hpp:186

inline void cvodes_check (int flag , const char *func_name );

// ./stan/math/prim/err/check_flag_sundials.hpp:195

inline std ::array <std ::string ,2>idas_flag_msg (int flag );

// ./stan/math/prim/err/check_flag_sundials.hpp:362

inline void idas_check (int flag , const char *func_name );

// ./stan/math/prim/err/check_flag_sundials.hpp:380

inline void kinsol_check (int flag , const char *func_name );

// ./stan/math/prim/err/check_flag_sundials.hpp:401

inline void kinsol_check (int flag , const char *func_name , long int max_num_steps );

// ./stan/math/prim/err/check_greater.hpp:34
template <typename T_y,
    typename T_low,
    require_all_stan_scalar_t<T_y,
    T_low>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater.hpp:67
template <typename T_y,
    typename T_low,
    require_stan_scalar_t<T_y>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater.hpp:106
template <typename T_y,
    typename T_low,
    require_stan_scalar_t<T_y>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater.hpp:146
template <typename T_y,
    typename T_low,
    require_vector_t<T_y>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater.hpp:184
template <typename T_y,
    typename T_low,
    require_dense_dynamic_t<T_y>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater.hpp:226
template <typename T_y,
    typename T_low,
    require_all_vector_t<T_y,
    T_low>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater.hpp:269
template <typename T_y,
    typename T_low,
    require_all_dense_dynamic_t<T_y,
    T_low>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater.hpp:309
template <typename T_y,
    typename T_low,
    require_std_vector_vt<is_container_or_var_matrix,
    T_y>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater.hpp:336
template <typename T_y,
    typename T_low,
    require_not_std_vector_t<T_y>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater.hpp:365
template <typename T_y,
    typename T_low,
    require_any_std_vector_vt<is_container_or_var_matrix,
    T_y,
    T_low>*>
inline void check_greater (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:34
template <typename T_y,
    typename T_low,
    require_all_stan_scalar_t<T_y,
    T_low>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:69
template <typename T_y,
    typename T_low,
    require_stan_scalar_t<T_y>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:110
template <typename T_y,
    typename T_low,
    require_stan_scalar_t<T_y>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:151
template <typename T_y,
    typename T_low,
    require_vector_t<T_y>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:191
template <typename T_y,
    typename T_low,
    require_dense_dynamic_t<T_y>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:234
template <typename T_y,
    typename T_low,
    require_all_vector_t<T_y,
    T_low>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:279
template <typename T_y,
    typename T_low,
    require_all_dense_dynamic_t<T_y,
    T_low>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:320
template <typename T_y,
    typename T_low,
    require_std_vector_vt<is_container_or_var_matrix,
    T_y>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:348
template <typename T_y,
    typename T_low,
    require_not_std_vector_t<T_y>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/err/check_greater_or_equal.hpp:376
template <typename T_y,
    typename T_low,
    require_any_std_vector_vt<is_container_or_var_matrix,
    T_y,
    T_low>*>
inline void check_greater_or_equal (const char *function , const char *name , const T_y &y , const T_low &low , Idxs ...idxs );

// ./stan/math/prim/fun/log.hpp:46
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto log (const Container &x );

// ./stan/math/prim/fun/log.hpp:63
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto log (const Container &x );

// ./stan/math/prim/fun/LDLT_factor.hpp:61
template <typename T,
    require_matrix_t<T>*>
inline auto make_ldlt_factor (const T &A );

// ./stan/math/prim/err/check_ldlt_factor.hpp:24
template <typename T>
inline void check_ldlt_factor (const char *function , const char *name , LDLT_factor <T >&A );

// ./stan/math/prim/err/check_less.hpp:33
template <typename T_y,
    typename T_high,
    require_all_stan_scalar_t<T_y,
    T_high>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less.hpp:66
template <typename T_y,
    typename T_high,
    require_stan_scalar_t<T_y>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less.hpp:105
template <typename T_y,
    typename T_high,
    require_stan_scalar_t<T_y>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less.hpp:145
template <typename T_y,
    typename T_high,
    require_vector_t<T_y>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less.hpp:183
template <typename T_y,
    typename T_high,
    require_dense_dynamic_t<T_y>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less.hpp:225
template <typename T_y,
    typename T_high,
    require_all_vector_t<T_y,
    T_high>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less.hpp:268
template <typename T_y,
    typename T_high,
    require_all_dense_dynamic_t<T_y,
    T_high>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less.hpp:308
template <typename T_y,
    typename T_high,
    require_std_vector_vt<is_container_or_var_matrix,
    T_y>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less.hpp:335
template <typename T_y,
    typename T_high,
    require_not_std_vector_t<T_y>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_less.hpp:364
template <typename T_y,
    typename T_high,
    require_any_std_vector_vt<is_container_or_var_matrix,
    T_y,
    T_high>*>
inline void check_less (const char *function , const char *name , const T_y &y , const T_high &high , Idxs ...idxs );

// ./stan/math/prim/err/check_multiplicable.hpp:29
template <typename T1,
    typename T2>
inline void check_multiplicable (const char *function , const char *name1 , const T1 &y1 , const char *name2 , const T2 &y2 );

// ./stan/math/prim/err/check_nonnegative.hpp:23
template <typename T_y>
inline void check_nonnegative (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/err/check_ordered.hpp:28
template <typename T_y,
    require_vector_t<T_y>*>
void check_ordered (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/err/check_ordered.hpp:60
template <typename T_y,
    require_std_vector_vt<is_stan_scalar,
    T_y>*>
void check_ordered (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/err/check_ordered.hpp:90
template <typename T_y,
    require_std_vector_t<T_y>*>
void check_ordered (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/err/check_sorted.hpp:25
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
void check_sorted (const char *function , const char *name , const EigVec &y );

// ./stan/math/prim/err/check_sorted.hpp:55
template <typename T_y>
void check_sorted (const char *function , const char *name , const std ::vector <T_y >&y );

// ./stan/math/prim/err/check_pos_semidefinite.hpp:29
template <typename EigMat,
    require_matrix_t<EigMat>*>
inline void check_pos_semidefinite (const char *function , const char *name , const EigMat &y );

// ./stan/math/prim/err/check_pos_semidefinite.hpp:61
template <typename Derived>
inline void check_pos_semidefinite (const char *function , const char *name , const Eigen ::LDLT <Derived >&cholesky );

// ./stan/math/prim/err/check_positive_finite.hpp:21
template <typename T_y>
inline void check_positive_finite (const char *function , const char *name , const T_y &y );

// ./stan/math/prim/err/check_positive_ordered.hpp:29
template <typename Vec,
    require_vector_t<Vec>*>
void check_positive_ordered (const char *function , const char *name , const Vec &y );

// ./stan/math/prim/err/check_positive_ordered.hpp:62
template <typename StdVec,
    require_std_vector_t<StdVec>*>
void check_positive_ordered (const char *function , const char *name , const StdVec &y );

// ./stan/math/prim/err/check_range.hpp:24

inline void check_range (const char *function , const char *name , int max , int index , int nested_level , const char *error_msg );

// ./stan/math/prim/err/check_range.hpp:49

inline void check_range (const char *function , const char *name , int max , int index , const char *error_msg );

// ./stan/math/prim/err/check_range.hpp:68

inline void check_range (const char *function , const char *name , int max , int index );

// ./stan/math/prim/err/check_row_index.hpp:24
template <typename T_y,
    typename >
inline void check_row_index (const char *function , const char *name , const T_y &y , size_t i );

// ./stan/math/prim/err/check_simplex.hpp:33
template <typename T,
    require_matrix_t<T>*>
void check_simplex (const char *function , const char *name , const T &theta );

// ./stan/math/prim/err/check_simplex.hpp:80
template <typename T,
    require_std_vector_t<T>*>
void check_simplex (const char *function , const char *name , const T &theta );

// ./stan/math/prim/err/check_std_vector_index.hpp:25
template <typename T>
inline void check_std_vector_index (const char *function , const char *name , const std ::vector <T >&y , int i );

// ./stan/math/prim/err/check_stochastic_column.hpp:34
template <typename T,
    require_matrix_t<T>*>
void check_stochastic_column (const char *function , const char *name , const T &theta );

// ./stan/math/prim/err/check_stochastic_column.hpp:88
template <typename T,
    require_std_vector_t<T>*>
void check_stochastic_column (const char *function , const char *name , const T &theta );

// ./stan/math/prim/err/check_stochastic_row.hpp:34
template <typename T,
    require_matrix_t<T>*>
void check_stochastic_row (const char *function , const char *name , const T &theta );

// ./stan/math/prim/err/check_stochastic_row.hpp:88
template <typename T,
    require_std_vector_t<T>*>
void check_stochastic_row (const char *function , const char *name , const T &theta );

// ./stan/math/prim/err/check_sum_to_zero.hpp:30
template <typename T,
    require_matrix_t<T>*>
void check_sum_to_zero (const char *function , const char *name , const T &theta );

// ./stan/math/prim/err/check_sum_to_zero.hpp:62
template <typename T,
    require_std_vector_t<T>*>
void check_sum_to_zero (const char *function , const char *name , const T &theta );

// ./stan/math/prim/err/check_vector.hpp:29
template <typename Mat,
    require_any_t<is_matrix<Mat>,
    is_prim_or_rev_kernel_expression<Mat>>*>
inline void check_vector (const char *function , const char *name , const Mat &x );

// ./stan/math/prim/err/check_vector_index.hpp:24
template <typename T,
    require_any_t<is_vector<T>,
    is_prim_or_rev_kernel_expression<T>>*>
inline void check_vector_index (const char *function , const char *name , const T &y , size_t i );

// ./stan/math/prim/err/domain_error_vec.hpp:12
template <typename T>
inline void domain_error_vec (const char *function , const char *name , const T &y , size_t i , const char *msg1 , const char *msg2 );

// ./stan/math/prim/err/domain_error_vec.hpp:21
template <typename T>
inline void domain_error_vec (const char *function , const char *name , const T &y , size_t i , const char *msg1 );

// ./stan/math/prim/err/hmm_check.hpp:27
template <typename T_omega,
    typename T_Gamma,
    typename T_rho,
    require_all_eigen_t<T_omega,
    T_Gamma>*>
inline void hmm_check (const T_omega &log_omegas , const T_Gamma &Gamma , const T_rho &rho , const char *function );

// ./stan/math/prim/err/invalid_argument_vec.hpp:32
template <typename T>
inline void invalid_argument_vec (const char *function , const char *name , const T &y , size_t i , const char *msg1 , const char *msg2 );

// ./stan/math/prim/err/invalid_argument_vec.hpp:60
template <typename T>
inline void invalid_argument_vec (const char *function , const char *name , const T &y , size_t i , const char *msg );

// ./stan/math/prim/err/is_positive.hpp:19
template <typename T_y>
inline bool is_positive (const T_y &y );

// ./stan/math/prim/err/is_positive.hpp:34

inline bool is_positive (int size );

// ./stan/math/prim/err/is_not_nan.hpp:22
template <typename T_y>
inline bool is_not_nan (const T_y &y );

// ./stan/math/prim/err/is_lower_triangular.hpp:20
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
inline bool is_lower_triangular (const EigMat &y );

// ./stan/math/prim/err/is_less_or_equal.hpp:25
template <typename T_y,
    typename T_high>
inline bool is_less_or_equal (const T_y &y , const T_high &high );

// ./stan/math/prim/err/is_cholesky_factor.hpp:28
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
inline bool is_cholesky_factor (const EigMat &y );

// ./stan/math/prim/err/is_nonzero_size.hpp:14
template <typename T_y>
inline bool is_nonzero_size (const T_y &y );

// ./stan/math/prim/err/is_unit_vector.hpp:27
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
inline bool is_unit_vector (const EigVec &theta );

// ./stan/math/prim/err/is_cholesky_factor_corr.hpp:25
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
inline bool is_cholesky_factor_corr (const EigMat &y );

// ./stan/math/prim/err/is_column_index.hpp:19
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline bool is_column_index (const EigMat &y , size_t i );

// ./stan/math/prim/err/is_size_match.hpp:16
template <typename T_size1,
    typename T_size2>
inline bool is_size_match (T_size1 i , T_size2 j );

// ./stan/math/prim/err/is_square.hpp:18
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline bool is_square (const EigMat &y );

// ./stan/math/prim/err/is_symmetric.hpp:23
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
inline bool is_symmetric (const EigMat &y );

// ./stan/math/prim/err/is_pos_definite.hpp:26
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
inline bool is_pos_definite (const EigMat &y );

// ./stan/math/prim/err/is_pos_definite.hpp:55
template <typename Derived>
inline bool is_pos_definite (const Eigen ::LDLT <Derived >&cholesky );

// ./stan/math/prim/err/is_pos_definite.hpp:70
template <typename Derived>
inline bool is_pos_definite (const Eigen ::LLT <Derived >&cholesky );

// ./stan/math/prim/err/is_corr_matrix.hpp:29
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
inline bool is_corr_matrix (const EigMat &y );

// ./stan/math/prim/err/is_ldlt_factor.hpp:20
template <typename T>
inline bool is_ldlt_factor (LDLT_factor <T >&A );

// ./stan/math/prim/err/is_mat_finite.hpp:19
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline bool is_mat_finite (const EigMat &y );

// ./stan/math/prim/err/is_matching_dims.hpp:20
template <typename EigMat1,
    typename EigMat2,
    require_all_matrix_t<EigMat1,
    EigMat2>*>
inline bool is_matching_dims (const EigMat1 &y1 , const EigMat2 &y2 );

// ./stan/math/prim/err/is_matching_dims.hpp:39
template <bool check_compile,
    typename EigMat1,
    typename EigMat2,
    require_all_matrix_t<EigMat1,
    EigMat2>*>
inline bool is_matching_dims (const EigMat1 &y1 , const EigMat2 &y2 );

// ./stan/math/prim/err/is_matching_size.hpp:23
template <typename T_y1,
    typename T_y2>
inline bool is_matching_size (const T_y1 &y1 , const T_y2 &y2 );

// ./stan/math/prim/err/is_ordered.hpp:18
template <typename T_y>
inline bool is_ordered (const std ::vector <T_y >&y );

// ./stan/math/prim/err/is_scal_finite.hpp:21
template <typename T_y>
inline bool is_scal_finite (const T_y &y );

// ./stan/math/prim/err/system_error.hpp:25

inline void system_error (const char *function , const char *name , const int &y , const char *msg1 , const char *msg2 );

// ./stan/math/prim/err/system_error.hpp:45

inline void system_error (const char *function , const char *name , const int &y , const char *msg1 );

// ./stan/math/prim/err/validate_non_negative_index.hpp:12

inline void validate_non_negative_index (const char *var_name , const char *expr , int val );

// ./stan/math/prim/err/validate_positive_index.hpp:20

inline void validate_positive_index (const char *var_name , const char *expr , int val );

// ./stan/math/prim/err/validate_unit_vector_index.hpp:20

inline void validate_unit_vector_index (const char *var_name , const char *expr , int val );

// ./stan/math/rev/core/recover_memory.hpp:18

static inline void recover_memory ();

// ./stan/math/rev/core/set_zero_all_adjoints.hpp:14

static inline void set_zero_all_adjoints ();

// ./stan/math/rev/core/save_varis.hpp:15
template <typename ...Pargs>
inline vari **save_varis (vari **dest , const var &x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:18
template <typename VarVec,
    require_std_vector_vt<is_var,
    VarVec>*>
inline vari **save_varis (vari **dest , VarVec &&x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:22
template <typename VecContainer,
    require_std_vector_st<is_var,
    VecContainer>*>
inline vari **save_varis (vari **dest , VecContainer &&x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:28
template <typename EigT,
    require_eigen_vt<is_var,
    EigT>*>
inline vari **save_varis (vari **dest , EigT &&x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:32
template <typename Arith,
    require_st_arithmetic<Arith>*>
inline vari **save_varis (vari **dest , Arith &&x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:36

inline vari **save_varis (vari **dest );

// ./stan/math/rev/core/save_varis.hpp:50
template <typename ...Pargs>
inline vari **save_varis (vari **dest , const var &x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:69
template <typename VarVec,
    require_std_vector_vt<is_var,
    VarVec>*,
    typename ...Pargs>
inline vari **save_varis (vari **dest , VarVec &&x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:91
template <typename VecContainer,
    require_std_vector_st<is_var,
    VecContainer>*,
    require_std_vector_vt<is_container,
    VecContainer>*,
    typename ...Pargs>
inline vari **save_varis (vari **dest , VecContainer &&x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:113
template <typename EigT,
    require_eigen_vt<is_var,
    EigT>*,
    typename ...Pargs>
inline vari **save_varis (vari **dest , EigT &&x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:134
template <typename Arith,
    require_st_arithmetic<Arith>*,
    typename ...Pargs>
inline vari **save_varis (vari **dest , Arith &&x , Pargs &&...args );

// ./stan/math/rev/core/save_varis.hpp:144

inline vari **save_varis (vari **dest );

// ./stan/math/rev/core/zero_adjoints.hpp:14

inline void zero_adjoints ()noexcept ;

// ./stan/math/rev/core/zero_adjoints.hpp:26
template <typename T,
    require_st_arithmetic<T>*>
inline void zero_adjoints (T &x )noexcept ;

// ./stan/math/rev/core/zero_adjoints.hpp:39

inline void zero_adjoints (var &x );

// ./stan/math/rev/core/zero_adjoints.hpp:51
template <typename EigMat,
    require_eigen_vt<is_autodiff,
    EigMat>*>
inline void zero_adjoints (EigMat &x );

// ./stan/math/rev/core/zero_adjoints.hpp:67
template <typename StdVec,
    require_std_vector_st<is_autodiff,
    StdVec>*>
inline void zero_adjoints (StdVec &x );

// ./stan/math/fwd/fun/read_fvar.hpp:43
template <typename EigFvar,
    typename EigOut>
inline void read_fvar (const EigFvar &FvarMat , EigOut &ValMat , EigOut &DMat );

// ./stan/math/mix/functor/derivative.hpp:23
template <typename T,
    typename F>
void derivative (const F &f , const T &x , T &fx , T &dfx_dx );

// ./stan/math/mix/functor/hessian.hpp:41
template <typename F>
void hessian (const F &f , const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&x , double &fx , Eigen ::Matrix <double , Eigen ::Dynamic , 1>&grad , Eigen ::Matrix <double , Eigen ::Dynamic , Eigen ::Dynamic >&H );

// ./stan/math/mix/functor/finite_diff_grad_hessian.hpp:40
template <typename F>
void finite_diff_grad_hessian (const F &f , const Eigen ::VectorXd &x , double &fx , Eigen ::MatrixXd &hess , std ::vector <Eigen ::MatrixXd >&grad_hess_fx , double epsilon =1e-04);

// ./stan/math/prim/fun/finite_diff_stepsize.hpp:21

inline double finite_diff_stepsize (double u );

// ./stan/math/prim/functor/finite_diff_gradient_auto.hpp:49
template <typename F,
    typename VectorT,
    typename GradVectorT,
    typename ScalarT>
void finite_diff_gradient_auto (const F &f , VectorT &&x , ScalarT &fx , GradVectorT &grad_fx );

// ./stan/math/mix/functor/finite_diff_grad_hessian_auto.hpp:43
template <typename F>
void finite_diff_grad_hessian_auto (const F &f , const Eigen ::VectorXd &x , double &fx , Eigen ::MatrixXd &hess , std ::vector <Eigen ::MatrixXd >&grad_hess_fx );

// ./stan/math/mix/functor/grad_hessian.hpp:41
template <typename F>
void grad_hessian (const F &f , const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&x , double &fx , Eigen ::Matrix <double , Eigen ::Dynamic , Eigen ::Dynamic >&H , std ::vector <Eigen ::Matrix <double , Eigen ::Dynamic , Eigen ::Dynamic >>&grad_H );

// ./stan/math/mix/functor/gradient_dot_vector.hpp:12
template <typename T1,
    typename T2,
    typename F>
void gradient_dot_vector (const F &f , const Eigen ::Matrix <T1 , Eigen ::Dynamic , 1>&x , const Eigen ::Matrix <T2 , Eigen ::Dynamic , 1>&v , T1 &fx , T1 &grad_fx_dot_v );

// ./stan/math/mix/functor/grad_tr_mat_times_hessian.hpp:14
template <typename F>
void grad_tr_mat_times_hessian (const F &f , const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&x , const Eigen ::Matrix <double , Eigen ::Dynamic , Eigen ::Dynamic >&M , Eigen ::Matrix <double , Eigen ::Dynamic , 1>&grad_tr_MH );

// ./stan/math/mix/functor/hessian_times_vector.hpp:13
template <typename F>
void hessian_times_vector (const F &f , const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&x , const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&v , double &fx , Eigen ::Matrix <double , Eigen ::Dynamic , 1>&Hv );

// ./stan/math/mix/functor/hessian_times_vector.hpp:38
template <typename T,
    typename F>
void hessian_times_vector (const F &f , const Eigen ::Matrix <T , Eigen ::Dynamic , 1>&x , const Eigen ::Matrix <T , Eigen ::Dynamic , 1>&v , T &fx , Eigen ::Matrix <T , Eigen ::Dynamic , 1>&Hv );

// ./stan/math/mix/functor/partial_derivative.hpp:24
template <typename T,
    typename F>
void partial_derivative (const F &f , const Eigen ::Matrix <T , Eigen ::Dynamic , 1>&x , int n , T &fx , T &dfx_dxn );

// ./stan/math/prim/fun/transpose.hpp:16
template <typename T,
    require_matrix_t<T>*>
auto inline transpose (const T &m );

// ./stan/math/prim/fun/dot_product.hpp:21
template <typename Vec1,
    typename Vec2,
    require_all_eigen_vector_t<Vec1,
    Vec2>*
;

// ./stan/math/prim/fun/dot_product.hpp:36
template <typename Scalar1,
    typename Scalar2,
    require_all_stan_scalar_t<Scalar1,
    Scalar2>*>
inline auto dot_product (const Scalar1 *v1 , const Scalar2 *v2 , size_t length );

// ./stan/math/prim/fun/dot_product.hpp:54
template <typename Scalar1,
    typename Scalar2,
    typename Alloc1,
    typename Alloc2,
    require_all_stan_scalar_t<Scalar1,
    Scalar2>*>
inline return_type_t <Scalar1 ,Scalar2 >dot_product (const std ::vector <Scalar1 , Alloc1 >&v1 , const std ::vector <Scalar2 , Alloc2 >&v2 );

// ./stan/math/fwd/fun/multiply.hpp:14
template <typename Mat1,
    typename Mat2,
    require_all_eigen_vt<is_fvar,
    Mat1,
    Mat2>*>
inline auto multiply (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/fwd/fun/multiply.hpp:23
template <typename Mat1,
    typename Mat2,
    require_eigen_vt<is_fvar,
    Mat1>*>
inline auto multiply (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/fwd/fun/multiply.hpp:42
template <typename Mat1,
    typename Mat2,
    require_eigen_vt<std::is_floating_point,
    Mat1>*>
inline auto multiply (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/fwd/fun/tcrossprod.hpp:13
template <typename EigMat,
    require_eigen_vt<is_fvar,
    EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,EigMat ::RowsAtCompileTime ,EigMat ::RowsAtCompileTime >tcrossprod (const EigMat &m );

// ./stan/math/prim/fun/sqrt.hpp:36
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto sqrt (const Container &x );

// ./stan/math/prim/fun/sqrt.hpp:53
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto sqrt (const Container &x );

// ./stan/math/prim/fun/inv_sqrt.hpp:15
template <typename T,
    require_stan_scalar_t<T>*>
inline auto inv_sqrt (T x );

// ./stan/math/prim/fun/inv_sqrt.hpp:42
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto inv_sqrt (const Container &x );

// ./stan/math/prim/fun/inv_sqrt.hpp:60
template <typename Container,
    require_not_var_matrix_t<Container>*>
inline auto inv_sqrt (const Container &x );

// ./stan/math/fwd/fun/inv_sqrt.hpp:11
template <typename T>
inline fvar <T >inv_sqrt (const fvar <T >&x );

// ./stan/math/fwd/fun/sqrt.hpp:16
template <typename T>
inline fvar <T >sqrt (const fvar <T >&x );

// ./stan/math/fwd/fun/sqrt.hpp:32
template <typename T>
inline std ::complex <fvar <T >>sqrt (const std ::complex <fvar <T >>&z );

// ./stan/math/prim/fun/divide.hpp:22
template <typename Scal1,
    typename Scal2,
    require_all_stan_scalar_t<Scal1,
    Scal2>*>
inline return_type_t <Scal1 ,Scal2 >divide (const Scal1 &x , const Scal2 &y );

// ./stan/math/prim/fun/divide.hpp:28

inline int divide (int x , int y );

// ./stan/math/prim/fun/divide.hpp:44
template <typename T1,
    typename T2,
    require_any_eigen_t<T1,
    T2>*>
inline auto divide (const T1 &m , const T2 &c );

// ./stan/math/prim/fun/dot_self.hpp:13

inline double dot_self (const std ::vector <double >&x );

// ./stan/math/prim/fun/dot_self.hpp:28
template <typename T,
    require_eigen_t<T>*>
inline value_type_t <T >dot_self (const T &v );

// ./stan/math/prim/constraint/unit_vector_constrain.hpp:25
template <typename T,
    require_eigen_col_vector_t<T>*>
inline plain_type_t <T >unit_vector_constrain (const T &y );

// ./stan/math/prim/constraint/unit_vector_constrain.hpp:47
template <typename T1,
    typename T2,
    require_eigen_col_vector_t<T1>*>
inline plain_type_t <T1 >unit_vector_constrain (const T1 &y , T2 &lp );

// ./stan/math/prim/constraint/unit_vector_constrain.hpp:75
template <bool Jacobian,
    typename T,
    require_not_std_vector_t<T>*>
inline auto unit_vector_constrain (const T &y , return_type_t <T >&lp );

// ./stan/math/prim/constraint/unit_vector_constrain.hpp:100
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto unit_vector_constrain (const T &y , return_type_t <T >&lp );

// ./stan/math/prim/fun/tcrossprod.hpp:20
template <typename T,
    require_eigen_vt<std::is_arithmetic,
    T>*>
inline Eigen ::Matrix <value_type_t <T >,T ::RowsAtCompileTime ,T ::RowsAtCompileTime >tcrossprod (const T &M );

// ./stan/math/fwd/constraint/unit_vector_constrain.hpp:18
template <typename EigMat,
    require_eigen_col_vector_vt<is_fvar,
    EigMat>*>
inline auto unit_vector_constrain (const EigMat &y );

// ./stan/math/fwd/constraint/unit_vector_constrain.hpp:41
template <typename EigMat,
    typename T,
    require_eigen_vt<is_fvar,
    EigMat>*>
inline auto unit_vector_constrain (const EigMat &y , T &lp );

// ./stan/math/fwd/fun/abs.hpp:14
template <typename T>
inline fvar <T >abs (const fvar <T >&x );

// ./stan/math/fwd/fun/abs.hpp:34
template <typename T>
inline fvar <T >abs (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/sum.hpp:21
template <typename T>
inline fvar <T >sum (const std ::vector <fvar <T >>&m );

// ./stan/math/fwd/fun/sum.hpp:43
template <typename T,
    require_eigen_vt<is_fvar,
    T>*>
inline value_type_t <T >sum (const T &m );

// ./stan/math/prim/fun/arg.hpp:16
template <typename V>
inline V arg (const std ::complex <V >&z );

// ./stan/math/prim/fun/copysign.hpp:29
template <typename T,
    typename U>
inline T copysign (const T &x , const U &y );

// ./stan/math/prim/fun/copysign.hpp:48
template <typename T,
    typename U>
inline T copysign_non_zero (const T &x , const U &y );

// ./stan/math/prim/fun/copysign.hpp:71
template <typename T,
    typename U>
inline std ::complex <T >copysign (const std ::complex <T >&x , const std ::complex <U >&y );

// ./stan/math/prim/fun/isfinite.hpp:20
template <typename ADType,
    require_autodiff_t<ADType>*>
inline bool isfinite (ADType &&v );

// ./stan/math/prim/fun/exp.hpp:42
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto exp (const Container &x );

// ./stan/math/prim/fun/exp.hpp:59
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto exp (const Container &x );

// ./stan/math/prim/fun/cosh.hpp:38
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto cosh (const Container &x );

// ./stan/math/prim/fun/cosh.hpp:55
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto cosh (const Container &x );

// ./stan/math/prim/fun/i_times.hpp:19
template <typename T>
inline std ::complex <T >i_times (const std ::complex <T >&z );

// ./stan/math/prim/fun/i_times.hpp:35
template <typename T>
inline std ::complex <T >neg_i_times (const std ::complex <T >&z );

// ./stan/math/prim/fun/i_times.hpp:47
template <typename V>
inline std ::complex <V >complex_negate (const std ::complex <V >&z );

// ./stan/math/prim/fun/cos.hpp:39
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto cos (const Container &x );

// ./stan/math/prim/fun/cos.hpp:56
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto cos (const Container &x );

// ./stan/math/prim/fun/polar.hpp:39

inline std ::complex <double >polar (double r , double theta );

// ./stan/math/prim/fun/polar.hpp:42

inline std ::complex <double >polar (double r , int theta );

// ./stan/math/prim/fun/polar.hpp:45

inline std ::complex <double >polar (int r , double theta );

// ./stan/math/prim/fun/polar.hpp:48

inline std ::complex <double >polar (int r , int theta );

// ./stan/math/prim/fun/asinh.hpp:46
template <typename T,
    require_not_var_matrix_t<T>*>
inline auto asinh (const T &x );

// ./stan/math/prim/fun/asin.hpp:42
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto asin (const Container &x );

// ./stan/math/prim/fun/asin.hpp:59
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto asin (const Container &x );

// ./stan/math/prim/fun/acos.hpp:44
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto acos (const Container &x );

// ./stan/math/prim/fun/acos.hpp:61
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto acos (const Container &x );

// ./stan/math/fwd/fun/acos.hpp:14
template <typename T>
inline fvar <T >acos (const fvar <T >&x );

// ./stan/math/fwd/fun/acos.hpp:28
template <typename T>
inline std ::complex <fvar <T >>acos (const std ::complex <fvar <T >>&x );

// ./stan/math/prim/fun/acosh.hpp:26

inline double acosh (double x );

// ./stan/math/prim/fun/acosh.hpp:46

inline double acosh (int x );

// ./stan/math/prim/fun/acosh.hpp:86
template <typename T,
    require_not_var_matrix_t<T>*>
inline auto acosh (const T &x );

// ./stan/math/fwd/fun/acosh.hpp:15
template <typename T>
inline fvar <T >acosh (const fvar <T >&x );

// ./stan/math/fwd/fun/acosh.hpp:28
template <typename T>
inline std ::complex <fvar <T >>acosh (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/asin.hpp:14
template <typename T>
inline fvar <T >asin (const fvar <T >&x );

// ./stan/math/fwd/fun/asin.hpp:28
template <typename T>
inline std ::complex <fvar <T >>asin (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/arg.hpp:18
template <typename T>
inline fvar <T >arg (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/asinh.hpp:14
template <typename T>
inline fvar <T >asinh (const fvar <T >&x );

// ./stan/math/fwd/fun/asinh.hpp:27
template <typename T>
inline std ::complex <fvar <T >>asinh (const std ::complex <fvar <T >>&z );

// ./stan/math/prim/fun/atanh.hpp:27

inline double atanh (double x );

// ./stan/math/prim/fun/atanh.hpp:43

inline double atanh (int x );

// ./stan/math/prim/fun/atanh.hpp:75
template <typename T,
    require_not_var_matrix_t<T>*>
inline auto atanh (const T &x );

// ./stan/math/prim/fun/atan.hpp:40
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto atan (const Container &x );

// ./stan/math/prim/fun/atan.hpp:57
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto atan (const Container &x );

// ./stan/math/fwd/fun/atan.hpp:14
template <typename T>
inline fvar <T >atan (const fvar <T >&x );

// ./stan/math/fwd/fun/atan.hpp:27
template <typename T>
inline std ::complex <fvar <T >>atan (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/atan2.hpp:12
template <typename T>
inline fvar <T >atan2 (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/atan2.hpp:20
template <typename T>
inline fvar <T >atan2 (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/atan2.hpp:27
template <typename T>
inline fvar <T >atan2 (const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/fun/atanh.hpp:23
template <typename T>
inline fvar <T >atanh (const fvar <T >&x );

// ./stan/math/fwd/fun/atanh.hpp:35
template <typename T>
inline std ::complex <fvar <T >>atanh (const std ::complex <fvar <T >>&z );

// ./stan/math/prim/fun/bessel_first_kind.hpp:38
template <typename T2,
    require_arithmetic_t<T2>*>
inline T2 bessel_first_kind (int v , const T2 z );

// ./stan/math/prim/fun/bessel_first_kind.hpp:54
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto bessel_first_kind (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/bessel_first_kind.hpp:11
template <typename T>
inline fvar <T >bessel_first_kind (int v , const fvar <T >&z );

// ./stan/math/prim/fun/bessel_second_kind.hpp:40
template <typename T2,
    require_arithmetic_t<T2>*>
inline T2 bessel_second_kind (int v , const T2 z );

// ./stan/math/prim/fun/bessel_second_kind.hpp:55
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto bessel_second_kind (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/bessel_second_kind.hpp:11
template <typename T>
inline fvar <T >bessel_second_kind (int v , const fvar <T >&z );

// ./stan/math/prim/fun/lgamma.hpp:63

inline double lgamma (double x );

// ./stan/math/prim/fun/lgamma.hpp:82

inline double lgamma (int x );

// ./stan/math/prim/fun/lgamma.hpp:117
template <typename T,
    require_not_var_matrix_t<T>*>
inline auto lgamma (const T &x );

// ./stan/math/prim/fun/beta.hpp:52
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline return_type_t <T1 ,T2 >beta (const T1 a , const T2 b );

// ./stan/math/prim/fun/beta.hpp:68
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto beta (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/digamma.hpp:47

inline double digamma (double x );

// ./stan/math/prim/fun/digamma.hpp:74
template <typename T,
    require_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto digamma (const T &x );

// ./stan/math/fwd/fun/beta.hpp:50
template <typename T>
inline fvar <T >beta (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/beta.hpp:59
template <typename T>
inline fvar <T >beta (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/beta.hpp:66
template <typename T>
inline fvar <T >beta (const fvar <T >&x1 , double x2 );

// ./stan/math/prim/fun/log1p.hpp:29

inline double log1p (double x );

// ./stan/math/prim/fun/log1p.hpp:47

inline double log1p (int x );

// ./stan/math/prim/fun/log1p.hpp:79
template <typename T,
    require_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto log1p (const T &x );

// ./stan/math/prim/fun/log1m.hpp:42

inline double log1m (double x );

// ./stan/math/prim/fun/log1m.hpp:70
template <typename T,
    require_not_var_matrix_t<T>*>
inline auto log1m (const T &x );

// ./stan/math/prim/fun/binary_log_loss.hpp:29
template <typename T,
    require_arithmetic_t<T>*>
inline T binary_log_loss (int y , const T &y_hat );

// ./stan/math/prim/fun/binary_log_loss.hpp:45
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto binary_log_loss (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/binary_log_loss.hpp:11
template <typename T>
inline fvar <T >binary_log_loss (int y , const fvar <T >&y_hat );

// ./stan/math/prim/fun/cbrt.hpp:34
template <typename T,
    require_not_var_matrix_t<T>*>
inline auto cbrt (const T &x );

// ./stan/math/fwd/fun/cbrt.hpp:19
template <typename T>
inline fvar <T >cbrt (const fvar <T >&x );

// ./stan/math/fwd/fun/ceil.hpp:11
template <typename T>
inline fvar <T >ceil (const fvar <T >&x );

// ./stan/math/prim/fun/conj.hpp:17
template <typename V>
inline std ::complex <V >conj (const std ::complex <V >&z );

// ./stan/math/prim/fun/conj.hpp:29
template <typename Eig,
    require_eigen_vt<is_complex,
    Eig>*>
inline auto conj (const Eig &z );

// ./stan/math/prim/fun/conj.hpp:41
template <typename StdVec,
    require_std_vector_st<is_complex,
    StdVec>*>
inline auto conj (const StdVec &z );

// ./stan/math/fwd/fun/conj.hpp:18
template <typename T>
inline std ::complex <fvar <T >>conj (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/cos.hpp:13
template <typename T>
inline fvar <T >cos (const fvar <T >&x );

// ./stan/math/fwd/fun/cos.hpp:27
template <typename T>
inline std ::complex <fvar <T >>cos (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/exp.hpp:12
template <typename T>
inline fvar <T >exp (const fvar <T >&x );

// ./stan/math/fwd/fun/exp.hpp:25
template <typename T>
inline std ::complex <fvar <T >>exp (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/cosh.hpp:14
template <typename T>
inline fvar <T >cosh (const fvar <T >&x );

// ./stan/math/fwd/fun/cosh.hpp:28
template <typename T>
inline std ::complex <fvar <T >>cosh (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/determinant.hpp:13
template <typename EigMat,
    require_eigen_vt<is_fvar,
    EigMat>*>
inline value_type_t <EigMat >determinant (const EigMat &m );

// ./stan/math/prim/fun/inv_square.hpp:19
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto inv_square (const Container &x );

// ./stan/math/prim/fun/inv_square.hpp:35
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto inv_square (const Container &x );

// ./stan/math/prim/fun/sinh.hpp:36
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto sinh (const Container &x );

// ./stan/math/prim/fun/sinh.hpp:53
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto sinh (const Container &x );

// ./stan/math/prim/fun/sin.hpp:38
template <typename T,
    require_not_container_st<std::is_arithmetic,
    T>*>
inline auto sin (const T &x );

// ./stan/math/prim/fun/sin.hpp:54
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto sin (const Container &x );

// ./stan/math/prim/fun/trigamma.hpp:36
template <typename T>
inline T trigamma_impl (const T &x );

// ./stan/math/prim/fun/trigamma.hpp:120

inline double trigamma (double u );

// ./stan/math/prim/fun/trigamma.hpp:129

inline double trigamma (int u );

// ./stan/math/prim/fun/trigamma.hpp:159
template <typename T,
    require_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto trigamma (const T &x );

// ./stan/math/fwd/fun/digamma.hpp:22
template <typename T>
inline fvar <T >digamma (const fvar <T >&x );

// ./stan/math/prim/fun/erf.hpp:34
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto erf (const T &x );

// ./stan/math/fwd/fun/erf.hpp:14
template <typename T>
inline fvar <T >erf (const fvar <T >&x );

// ./stan/math/prim/fun/erfc.hpp:35
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto erfc (const T &x );

// ./stan/math/fwd/fun/erfc.hpp:14
template <typename T>
inline fvar <T >erfc (const fvar <T >&x );

// ./stan/math/prim/fun/exp2.hpp:38
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto exp2 (const T &x );

// ./stan/math/fwd/fun/exp2.hpp:13
template <typename T>
inline fvar <T >exp2 (const fvar <T >&x );

// ./stan/math/prim/fun/expm1.hpp:35
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto expm1 (const T &x );

// ./stan/math/fwd/fun/expm1.hpp:12
template <typename T>
inline fvar <T >expm1 (const fvar <T >&x );

// ./stan/math/fwd/fun/fabs.hpp:14
template <typename T>
inline fvar <T >fabs (const fvar <T >&x );

// ./stan/math/prim/fun/falling_factorial.hpp:63
template <typename T,
    require_arithmetic_t<T>*>
inline return_type_t <T >falling_factorial (const T &x , int n );

// ./stan/math/prim/fun/falling_factorial.hpp:81
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto falling_factorial (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/falling_factorial.hpp:25
template <typename T>
inline fvar <T >falling_factorial (const fvar <T >&x , int n );

// ./stan/math/prim/fun/fdim.hpp:21
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline double fdim (T1 x , T2 y );

// ./stan/math/prim/fun/fdim.hpp:37
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto fdim (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/fdim.hpp:20
template <typename T>
inline fvar <T >fdim (const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/fun/fdim.hpp:38
template <typename T>
inline fvar <T >fdim (const fvar <T >&x , double y );

// ./stan/math/fwd/fun/fdim.hpp:56
template <typename T>
inline fvar <T >fdim (double x , const fvar <T >&y );

// ./stan/math/fwd/fun/floor.hpp:11
template <typename T>
inline fvar <T >floor (const fvar <T >&x );

// ./stan/math/prim/fun/fma.hpp:24
template <typename T1,
    typename T2,
    typename T3,
    require_all_arithmetic_t<T1,
    T2,
    T3>*>
inline double fma (T1 x , T2 y , T3 z );

// ./stan/math/prim/fun/fma.hpp:31
template <typename T1,
    typename T2,
    typename T3,
    require_any_matrix_t<T1,
    T2,
    T3>*
;

// ./stan/math/fwd/fun/fma.hpp:58
template <typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline fvar <return_type_t <T1 ,T2 ,T3 >>fma (const fvar <T1 >&x1 , const fvar <T2 >&x2 , const fvar <T3 >&x3 );

// ./stan/math/fwd/fun/fma.hpp:71
template <typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline fvar <return_type_t <T1 ,T2 ,T3 >>fma (const T1 &x1 , const fvar <T2 >&x2 , const fvar <T3 >&x3 );

// ./stan/math/fwd/fun/fma.hpp:82
template <typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline fvar <return_type_t <T1 ,T2 ,T3 >>fma (const fvar <T1 >&x1 , const T2 &x2 , const fvar <T3 >&x3 );

// ./stan/math/fwd/fun/fma.hpp:93
template <typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline fvar <return_type_t <T1 ,T2 ,T3 >>fma (const fvar <T1 >&x1 , const fvar <T2 >&x2 , const T3 &x3 );

// ./stan/math/fwd/fun/fma.hpp:104
template <typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline fvar <return_type_t <T1 ,T2 ,T3 >>fma (const T1 &x1 , const T2 &x2 , const fvar <T3 >&x3 );

// ./stan/math/fwd/fun/fma.hpp:114
template <typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline fvar <return_type_t <T1 ,T2 ,T3 >>fma (const fvar <T1 >&x1 , const T2 &x2 , const T3 &x3 );

// ./stan/math/fwd/fun/fma.hpp:124
template <typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline fvar <return_type_t <T1 ,T2 ,T3 >>fma (const T1 &x1 , const fvar <T2 >&x2 , const T3 &x3 );

// ./stan/math/prim/fun/fmax.hpp:19
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline double fmax (T1 x , T2 y );

// ./stan/math/prim/fun/fmax.hpp:35
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto fmax (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/fmax.hpp:22
template <typename T>
inline fvar <T >fmax (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/fmax.hpp:45
template <typename T>
inline fvar <T >fmax (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/fmax.hpp:68
template <typename T>
inline fvar <T >fmax (const fvar <T >&x1 , double x2 );

// ./stan/math/prim/fun/fmin.hpp:19
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline double fmin (T1 x , T2 y );

// ./stan/math/prim/fun/fmin.hpp:35
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto fmin (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/fmin.hpp:13
template <typename T>
inline fvar <T >fmin (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/fmin.hpp:27
template <typename T>
inline fvar <T >fmin (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/fmin.hpp:41
template <typename T>
inline fvar <T >fmin (const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/fun/fmod.hpp:14
template <typename T>
inline fvar <T >fmod (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/fmod.hpp:22
template <typename T>
inline fvar <T >fmod (const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/fun/fmod.hpp:32
template <typename T>
inline fvar <T >fmod (double x1 , const fvar <T >&x2 );

// ./stan/math/prim/fun/gamma_p.hpp:67

inline double gamma_p (double z , double a );

// ./stan/math/prim/fun/gamma_p.hpp:89
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto gamma_p (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/gamma_q.hpp:55

inline double gamma_q (double x , double a );

// ./stan/math/prim/fun/gamma_q.hpp:67
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto gamma_q (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/multiply_log.hpp:47
template <typename T_a,
    typename T_b,
    require_all_arithmetic_t<T_a,
    T_b>*>
inline return_type_t <T_a ,T_b >multiply_log (const T_a a , const T_b b );

// ./stan/math/prim/fun/multiply_log.hpp:68
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto multiply_log (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/grad_reg_inc_gamma.hpp:51
template <typename T1,
    typename T2>
return_type_t <T1 ,T2 >grad_reg_inc_gamma (T1 a , T2 z , T1 g , T1 dig , double precision =1e-6, int max_steps =1e5);

// ./stan/math/prim/fun/tgamma.hpp:19

inline double tgamma (double x );

// ./stan/math/prim/fun/tgamma.hpp:49
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto tgamma (const T &x );

// ./stan/math/prim/fun/grad_reg_lower_inc_gamma.hpp:109
template <typename T1,
    typename T2>
return_type_t <T1 ,T2 >grad_reg_lower_inc_gamma (const T1 &a , const T2 &z , double precision =1e-10, int max_steps =1e5);

// ./stan/math/fwd/fun/gamma_p.hpp:15
template <typename T>
inline fvar <T >gamma_p (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/gamma_p.hpp:34
template <typename T>
inline fvar <T >gamma_p (const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/fun/gamma_p.hpp:49
template <typename T>
inline fvar <T >gamma_p (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/gamma_q.hpp:13
template <typename T>
inline fvar <T >gamma_q (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/gamma_q.hpp:45
template <typename T>
inline fvar <T >gamma_q (const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/fun/gamma_q.hpp:76
template <typename T>
inline fvar <T >gamma_q (double x1 , const fvar <T >&x2 );

// ./stan/math/prim/fun/sign.hpp:11
template <typename T,
    require_stan_scalar_t<T>*>
inline int sign (const T &z );

// ./stan/math/prim/fun/sign.hpp:41
template <typename T,
    require_container_t<T>*>
inline auto sign (const T &x );

// ./stan/math/prim/fun/pow.hpp:61
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto pow (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/hypergeometric_pFq.hpp:28
template <typename Ta,
    typename Tb,
    typename Tz,
    require_all_eigen_st<std::is_arithmetic,
    Ta,
    Tb>*>
return_type_t <Ta ,Tb ,Tz >hypergeometric_pFq (const Ta &a , const Tb &b , const Tz &z );

// ./stan/math/prim/fun/hypergeometric_2F1.hpp:150
template <typename Ta1,
    typename Ta2,
    typename Tb,
    typename Tz,
    typename ScalarT>
inline return_type_t <Ta1 ,Ta1 ,Tb ,Tz >hypergeometric_2F1 (const Ta1 &a1 , const Ta2 &a2 , const Tb &b , const Tz &z );

// ./stan/math/prim/fun/grad_2F1.hpp:328
template <typename T1,
    typename T2,
    typename T3,
    typename T_z>
auto grad_2F1 (const T1 &a1 , const T2 &a2 , const T3 &b1 , const T_z &z , double precision =1e-14, int max_steps =1e6);

// ./stan/math/prim/functor/apply_scalar_ternary.hpp:37
template <typename F,
    typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline auto apply_scalar_ternary (const F &f , const T1 &x , const T2 &y , const T3 &z );

// ./stan/math/prim/functor/apply_scalar_ternary.hpp:59
template <typename F,
    typename T1,
    typename T2,
    typename T3,
    require_all_eigen_t<T1,
    T2,
    T3>*>
inline auto apply_scalar_ternary (F &&f , T1 &&x , T2 &&y , T3 &&z );

// ./stan/math/prim/functor/apply_scalar_ternary.hpp:94
template <typename F,
    typename T1,
    typename T2,
    typename T3,
    require_all_std_vector_vt<is_stan_scalar,
    T1,
    T2,
    T3>*>
inline auto apply_scalar_ternary (const F &f , const T1 &x , const T2 &y , const T3 &z );

// ./stan/math/prim/functor/apply_scalar_ternary.hpp:126
template <typename F,
    typename T1,
    typename T2,
    typename T3,
    require_all_std_vector_vt<is_container_or_var_matrix,
    T1,
    T2,
    T3>*>
inline auto apply_scalar_ternary (const F &f , const T1 &x , const T2 &y , const T3 &z );

// ./stan/math/prim/functor/apply_scalar_ternary.hpp:159
template <typename F,
    typename T1,
    typename T2,
    typename T3,
    require_any_container_t<T1,
    T2>*>
inline auto apply_scalar_ternary (const F &f , const T1 &x , const T2 &y , const T3 &z );

// ./stan/math/prim/functor/apply_scalar_ternary.hpp:184
template <typename F,
    typename T1,
    typename T2,
    typename T3,
    require_all_container_t<T1,
    T3>*>
inline auto apply_scalar_ternary (const F &f , const T1 &x , const T2 &y , const T3 &z );

// ./stan/math/prim/functor/apply_scalar_ternary.hpp:209
template <typename F,
    typename T1,
    typename T2,
    typename T3,
    require_container_t<T3>*>
inline auto apply_scalar_ternary (const F &f , const T1 &x , const T2 &y , const T3 &z );

// ./stan/math/prim/fun/inc_beta.hpp:25

inline double inc_beta (double a , double b , double x );

// ./stan/math/prim/fun/inc_beta.hpp:44
template <typename T1,
    typename T2,
    typename T3,
    require_any_container_t<T1,
    T2,
    T3>*>
inline auto inc_beta (const T1 &a , const T2 &b , const T3 &c );

// ./stan/math/prim/fun/grad_inc_beta.hpp:22

inline void grad_inc_beta (double &g1 , double &g2 , double a , double b , double z );

// ./stan/math/prim/fun/grad_reg_inc_beta.hpp:32
template <typename T>
void grad_reg_inc_beta (T &g1 , T &g2 , const T &a , const T &b , const T &z , const T &digammaA , const T &digammaB , const T &digammaSum , const T &betaAB );

// ./stan/math/prim/fun/lgamma_stirling.hpp:25
template <typename T>
return_type_t <T >lgamma_stirling (const T x );

// ./stan/math/prim/fun/lgamma_stirling_diff.hpp:42
template <typename T>
return_type_t <T >lgamma_stirling_diff (const T x );

// ./stan/math/prim/fun/lbeta.hpp:64
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
return_type_t <T1 ,T2 >lbeta (const T1 a , const T2 b );

// ./stan/math/prim/fun/lbeta.hpp:129
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto lbeta (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/inv.hpp:11
template <typename T>
inline fvar <T >inv (const fvar <T >&x );

// ./stan/math/fwd/fun/inv_square.hpp:11
template <typename T>
inline fvar <T >inv_square (const fvar <T >&x );

// ./stan/math/fwd/fun/pow.hpp:18
template <typename T>
inline fvar <T >pow (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/pow.hpp:27
template <typename T,
    typename U,
    typename >
inline fvar <T >pow (U x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/pow.hpp:35
template <typename T,
    typename U,
    typename >
inline fvar <T >pow (const fvar <T >&x1 , U x2 );

// ./stan/math/fwd/fun/pow.hpp:80
template <typename V>
inline std ::complex <fvar <V >>pow (const std ::complex <fvar <V >>&x , const std ::complex <fvar <V >>&y );

// ./stan/math/fwd/fun/pow.hpp:95
template <typename V,
    typename T,
    typename >
inline std ::complex <fvar <V >>pow (const std ::complex <fvar <V >>&x , const std ::complex <T >&y );

// ./stan/math/fwd/fun/pow.hpp:109
template <typename V>
inline std ::complex <fvar <V >>pow (const std ::complex <fvar <V >>&x , const fvar <V >&y );

// ./stan/math/fwd/fun/pow.hpp:124
template <typename V,
    typename T,
    typename >
inline std ::complex <fvar <V >>pow (const std ::complex <fvar <V >>&x , const T &y );

// ./stan/math/fwd/fun/pow.hpp:138
template <typename V,
    typename T,
    typename >
inline std ::complex <fvar <V >>pow (const std ::complex <T >&x , const std ::complex <fvar <V >>&y );

// ./stan/math/fwd/fun/pow.hpp:153
template <typename V,
    typename T,
    typename >
inline std ::complex <fvar <V >>pow (const std ::complex <T >&x , const fvar <V >&y );

// ./stan/math/fwd/fun/pow.hpp:166
template <typename V>
inline std ::complex <fvar <V >>pow (const fvar <V >&x , const std ::complex <fvar <V >>&y );

// ./stan/math/fwd/fun/pow.hpp:181
template <typename V,
    typename T,
    typename >
inline std ::complex <fvar <V >>pow (const fvar <V >&x , const std ::complex <T >&y );

// ./stan/math/fwd/fun/pow.hpp:195
template <typename T,
    typename V,
    typename >
inline std ::complex <fvar <V >>pow (T x , const std ::complex <fvar <V >>&y );

// ./stan/math/fwd/fun/pow.hpp:212
template <typename T>
inline std ::complex <fvar <T >>pow (const std ::complex <fvar <T >>&x , int y );

// ./stan/math/fwd/fun/inc_beta.hpp:18
template <typename T>
inline fvar <T >inc_beta (const fvar <T >&a , const fvar <T >&b , const fvar <T >&x );

// ./stan/math/fwd/fun/inc_beta.hpp:31
template <typename T>
inline fvar <T >inc_beta (double a , const fvar <T >&b , const fvar <T >&x );

// ./stan/math/fwd/fun/inc_beta.hpp:35
template <typename T>
inline fvar <T >inc_beta (const fvar <T >&a , double b , const fvar <T >&x );

// ./stan/math/fwd/fun/inc_beta.hpp:39
template <typename T>
inline fvar <T >inc_beta (const fvar <T >&a , const fvar <T >&b , double x );

// ./stan/math/fwd/fun/inc_beta.hpp:44
template <typename T>
inline fvar <T >inc_beta (double a , double b , const fvar <T >&x );

// ./stan/math/fwd/fun/inc_beta.hpp:48
template <typename T>
inline fvar <T >inc_beta (const fvar <T >&a , double b , double x );

// ./stan/math/fwd/fun/inc_beta.hpp:52
template <typename T>
inline fvar <T >inc_beta (double a , const fvar <T >&b , double x );

// ./stan/math/fwd/fun/log.hpp:14
template <typename T>
inline fvar <T >log (const fvar <T >&x );

// ./stan/math/fwd/fun/log.hpp:31
template <typename T>
inline std ::complex <fvar <T >>log (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/log1m.hpp:11
template <typename T>
inline fvar <T >log1m (const fvar <T >&x );

// ./stan/math/fwd/fun/value_of.hpp:17
template <typename T>
inline T value_of (const fvar <T >&v );

// ./stan/math/fwd/fun/grad_inc_beta.hpp:33
template <typename T>
void grad_inc_beta (fvar <T >&g1 , fvar <T >&g2 , fvar <T >a , fvar <T >b , fvar <T >z );

// ./stan/math/prim/fun/hypergeometric_1F0.hpp:30
template <typename Ta,
    typename Tz,
    require_all_arithmetic_t<Ta,
    Tz>*>
return_type_t <Ta ,Tz >hypergeometric_1f0 (const Ta &a , const Tz &z );

// ./stan/math/fwd/fun/hypergeometric_1F0.hpp:28
template <typename Ta,
    typename Tz,
    typename FvarT>
FvarT hypergeometric_1f0 (const Ta &a , const Tz &z );

// ./stan/math/fwd/fun/hypergeometric_2F1.hpp:29
template <typename Ta1,
    typename Ta2,
    typename Tb,
    typename Tz,
    require_all_stan_scalar_t<Ta1,
    Ta2,
    Tb,
    Tz>*>
inline return_type_t <Ta1 ,Ta1 ,Tb ,Tz >hypergeometric_2F1 (const Ta1 &a1 , const Ta2 &a2 , const Tb &b , const Tz &z );

// ./stan/math/prim/fun/add.hpp:20
template <typename ScalarA,
    typename ScalarB,
    require_all_stan_scalar_t<ScalarA,
    ScalarB>*>
inline return_type_t <ScalarA ,ScalarB >add (const ScalarA &a , const ScalarB &b );

// ./stan/math/prim/fun/add.hpp:40
template <typename Mat1,
    typename Mat2,
    require_all_eigen_t<Mat1,
    Mat2>*>
inline auto add (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/prim/fun/add.hpp:57
template <typename Mat,
    typename Scal,
    require_eigen_t<Mat>*>
inline auto add (const Mat &m , const Scal c );

// ./stan/math/prim/fun/add.hpp:73
template <typename Scal,
    typename Mat,
    require_stan_scalar_t<Scal>*>
inline auto add (const Scal c , const Mat &m );

// ./stan/math/prim/fun/prod.hpp:18
template <typename T,
    require_stan_scalar_t<T>*>
T prod (const T &v );

// ./stan/math/prim/fun/prod.hpp:31
template <typename T>
inline T prod (const std ::vector <T >&v );

// ./stan/math/prim/fun/prod.hpp:48
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline value_type_t <EigMat >prod (const EigMat &v );

// ./stan/math/prim/fun/max.hpp:24
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
auto max (T1 x , T2 y );

// ./stan/math/prim/fun/max.hpp:41
template <typename T,
    require_container_t<T>*>
inline value_type_t <T >max (const T &m );

// ./stan/math/prim/fun/select.hpp:23
template <typename T_true,
    typename T_false,
    typename ReturnT>
inline ReturnT select (const bool c , const T_true y_true , const T_false y_false );

// ./stan/math/prim/fun/select.hpp:47
template <typename T_true,
    typename T_false,
    typename T_return
;

// ./stan/math/prim/fun/select.hpp:146
template <typename T_bool,
    typename T_true,
    typename T_false,
    require_eigen_array_vt<std::is_integral,
    T_bool>*>
inline auto select (const T_bool c , const T_true y_true , const T_false y_false );

// ./stan/math/prim/fun/select.hpp:170
template <typename T_bool,
    typename T_true,
    typename T_false,
    require_eigen_array_t<T_bool>*>
inline auto select (const T_bool c , const T_true y_true , const T_false y_false );

// ./stan/math/prim/fun/grad_pFq.hpp:87
template <bool calc_a
;

// ./stan/math/fwd/fun/hypergeometric_pFq.hpp:24
template <typename Ta,
    typename Tb,
    typename Tz,
    typename FvarT>
inline FvarT hypergeometric_pFq (const Ta &a , const Tb &b , const Tz &z );

// ./stan/math/fwd/fun/hypot.hpp:25
template <typename T>
inline fvar <T >hypot (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/hypot.hpp:45
template <typename T>
inline fvar <T >hypot (const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/fun/hypot.hpp:65
template <typename T>
inline fvar <T >hypot (double x1 , const fvar <T >&x2 );

// ./stan/math/prim/fun/inv_erfc.hpp:18
template <typename T,
    require_arithmetic_t<T>*>
inline auto inv_erfc (const T &x );

// ./stan/math/prim/fun/inv_erfc.hpp:46
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto inv_erfc (const T &x );

// ./stan/math/fwd/fun/square.hpp:11
template <typename T>
inline fvar <T >square (const fvar <T >&x );

// ./stan/math/fwd/fun/inv_erfc.hpp:15
template <typename T>
inline fvar <T >inv_erfc (const fvar <T >&x );

// ./stan/math/prim/fun/Phi.hpp:32

inline double Phi (double x );

// ./stan/math/prim/fun/Phi.hpp:66
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto Phi (const T &x );

// ./stan/math/prim/fun/inv_Phi.hpp:147

inline double inv_Phi (double p );

// ./stan/math/prim/fun/inv_Phi.hpp:176
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto inv_Phi (const T &x );

// ./stan/math/fwd/fun/inv_Phi.hpp:14
template <typename T>
inline fvar <T >inv_Phi (const fvar <T >&p );

// ./stan/math/prim/fun/inv_cloglog.hpp:48

inline double inv_cloglog (double x );

// ./stan/math/prim/fun/inv_cloglog.hpp:74
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto inv_cloglog (const Container &x );

// ./stan/math/prim/fun/inv_cloglog.hpp:91
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto inv_cloglog (const Container &x );

// ./stan/math/fwd/fun/inv_cloglog.hpp:12
template <typename T>
inline fvar <T >inv_cloglog (const fvar <T >&x );

// ./stan/math/prim/fun/inv_inc_beta.hpp:25

inline double inv_inc_beta (double a , double b , double p );

// ./stan/math/prim/fun/inv_inc_beta.hpp:47
template <typename T1,
    typename T2,
    typename T3,
    require_any_container_t<T1,
    T2,
    T3>*>
inline auto inv_inc_beta (const T1 &a , const T2 &b , const T3 &c );

// ./stan/math/prim/fun/log1m_exp.hpp:47

inline double log1m_exp (double a );

// ./stan/math/prim/fun/log1m_exp.hpp:80
template <typename T,
    require_not_var_matrix_t<T>*>
inline auto log1m_exp (const T &x );

// ./stan/math/prim/fun/log_diff_exp.hpp:50
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline return_type_t <T1 ,T2 >log_diff_exp (const T1 x , const T2 y );

// ./stan/math/prim/fun/log_diff_exp.hpp:68
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto log_diff_exp (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/append_row.hpp:30
template <typename T1,
    typename T2,
    require_all_eigen_t<T1,
    T2>*>
inline auto append_row (const T1 &A , const T2 &B );

// ./stan/math/prim/fun/append_row.hpp:65
template <typename Scal,
    typename ColVec,
    require_stan_scalar_t<Scal>*
;

// ./stan/math/prim/fun/append_row.hpp:92
template <typename ColVec,
    typename Scal,
    require_t<is_eigen_col_vector<ColVec>>*>
inline Eigen ::Matrix <return_type_t <ColVec ,Scal >,Eigen ::Dynamic ,1>append_row (const ColVec &A , const Scal &B );

// ./stan/math/prim/fun/to_vector.hpp:14
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,1>to_vector (const EigMat &matrix );

// ./stan/math/prim/fun/to_vector.hpp:27
template <typename T>
inline Eigen ::Matrix <T ,Eigen ::Dynamic ,1>to_vector (const std ::vector <T >&vec );

// ./stan/math/prim/fun/to_vector.hpp:34

inline Eigen ::Matrix <double ,Eigen ::Dynamic ,1>to_vector (const std ::vector <int >&vec );

// ./stan/math/prim/fun/hypergeometric_3F2.hpp:115
template <typename Ta,
    typename Tb,
    typename Tz,
    require_all_vector_t<Ta,
    Tb>*>
auto hypergeometric_3F2 (const Ta &a , const Tb &b , const Tz &z );

// ./stan/math/prim/fun/hypergeometric_3F2.hpp:144
template <typename Ta,
    typename Tb,
    typename Tz,
    require_all_stan_scalar_t<Ta,
    Tb,
    Tz>*>
auto hypergeometric_3F2 (const std ::initializer_list <Ta >&a , const std ::initializer_list <Tb >&b , const Tz &z );

// ./stan/math/fwd/fun/inv_inc_beta.hpp:32
template <typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline fvar <partials_return_t <T1 ,T2 ,T3 >>inv_inc_beta (const T1 &a , const T2 &b , const T3 &p );

// ./stan/math/prim/fun/inv_logit.hpp:51

inline double inv_logit (double a );

// ./stan/math/prim/fun/inv_logit.hpp:84
template <typename T,
    require_not_var_matrix_t<T>*>
inline auto inv_logit (const T &x );

// ./stan/math/fwd/fun/inv_logit.hpp:19
template <typename T>
inline fvar <T >inv_logit (const fvar <T >&x );

// ./stan/math/prim/fun/multiply.hpp:23
template <typename Mat,
    typename Scal,
    require_stan_scalar_t<Scal>*
;

// ./stan/math/prim/fun/multiply.hpp:43
template <typename Mat,
    typename Scal,
    require_any_complex_t<value_type_t<Mat>,
    Scal>*>
inline auto multiply (const Mat &m , Scal c );

// ./stan/math/prim/fun/multiply.hpp:62
template <typename Mat,
    typename Scal,
    require_any_complex_t<value_type_t<Mat>,
    Scal>*>
inline auto multiply (const Scal &m , const Mat &c );

// ./stan/math/prim/fun/multiply.hpp:79
template <typename Scal,
    typename Mat,
    require_stan_scalar_t<Scal>*
;

// ./stan/math/prim/fun/multiply.hpp:101
template <typename Mat1,
    typename Mat2,
    require_all_eigen_vt<std::is_arithmetic,
    Mat1,
    Mat2>*>
inline auto multiply (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/prim/fun/multiply.hpp:124
template <typename Mat1,
    typename Mat2,
    require_all_eigen_t<Mat1,
    Mat2>*
;

// ./stan/math/prim/fun/multiply.hpp:147
template <typename RowVec,
    typename ColVec,
    require_not_var_t<return_type_t<RowVec,
    ColVec>>*>
inline auto multiply (const RowVec &rv , const ColVec &v );

// ./stan/math/prim/fun/multiply.hpp:164
template <typename Scalar1,
    typename Scalar2,
    require_all_stan_scalar_t<Scalar1,
    Scalar2>*>
inline return_type_t <Scalar1 ,Scalar2 >multiply (Scalar1 m , Scalar2 c );

// ./stan/math/prim/fun/inverse.hpp:23
template <typename EigMat,
    require_eigen_vt<std::is_arithmetic,
    EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,EigMat ::RowsAtCompileTime ,EigMat ::ColsAtCompileTime >inverse (const EigMat &m );

// ./stan/math/fwd/fun/to_fvar.hpp:13
template <typename T,
    require_stan_scalar_t<T>*>
inline fvar <T >to_fvar (const T &x );

// ./stan/math/fwd/fun/to_fvar.hpp:25
template <typename T,
    require_fvar_t<scalar_type_t<T>>*>
inline T &&to_fvar (T &&x );

// ./stan/math/fwd/fun/to_fvar.hpp:30
template <typename T>
inline std ::vector <fvar <T >>to_fvar (const std ::vector <T >&v );

// ./stan/math/fwd/fun/to_fvar.hpp:39
template <typename T>
inline std ::vector <fvar <T >>to_fvar (const std ::vector <T >&v , const std ::vector <T >&d );

// ./stan/math/fwd/fun/to_fvar.hpp:49
template <typename T,
    require_eigen_t<T>*>
inline promote_scalar_t <fvar <value_type_t <T >>,T >to_fvar (const T &m );

// ./stan/math/fwd/fun/to_fvar.hpp:58
template <typename T1,
    typename T2,
    require_all_eigen_t<T1,
    T2>*>
inline promote_scalar_t <fvar <value_type_t <T1 >>,T1 >to_fvar (const T1 &val , const T2 &deriv );

// ./stan/math/fwd/fun/inverse.hpp:26
template <typename EigMat,
    require_eigen_vt<is_fvar,
    EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,EigMat ::RowsAtCompileTime ,EigMat ::ColsAtCompileTime >inverse (const EigMat &m );

// ./stan/math/fwd/fun/is_inf.hpp:20
template <typename T>
inline int is_inf (const fvar <T >&x );

// ./stan/math/fwd/fun/is_nan.hpp:21
template <typename T,
    require_fvar_t<T>*>
inline bool is_nan (T &&x );

// ./stan/math/prim/fun/lambert_w.hpp:20
template <typename T,
    require_arithmetic_t<T>*>
inline double lambert_w0 (const T &x );

// ./stan/math/prim/fun/lambert_w.hpp:34
template <typename T,
    require_arithmetic_t<T>*>
inline double lambert_wm1 (const T &x );

// ./stan/math/prim/fun/lambert_w.hpp:81
template <typename T,
    require_not_stan_scalar_t<T>*>
inline auto lambert_w0 (const T &x );

// ./stan/math/prim/fun/lambert_w.hpp:96
template <typename T,
    require_not_stan_scalar_t<T>*>
inline auto lambert_wm1 (const T &x );

// ./stan/math/fwd/fun/lambert_w.hpp:12
template <typename T>
inline fvar <T >lambert_w0 (const fvar <T >&x );

// ./stan/math/fwd/fun/lambert_w.hpp:18
template <typename T>
inline fvar <T >lambert_wm1 (const fvar <T >&x );

// ./stan/math/fwd/fun/lbeta.hpp:13
template <typename T>
inline fvar <T >lbeta (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/lbeta.hpp:20
template <typename T>
inline fvar <T >lbeta (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/lbeta.hpp:26
template <typename T>
inline fvar <T >lbeta (const fvar <T >&x1 , double x2 );

// ./stan/math/prim/fun/ldexp.hpp:20
template <typename T1,
    require_arithmetic_t<T1>*>
inline double ldexp (T1 a , int b );

// ./stan/math/prim/fun/ldexp.hpp:36
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto ldexp (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/ldexp.hpp:20
template <typename T>
inline fvar <T >ldexp (const fvar <T >&a , int b );

// ./stan/math/fwd/fun/lgamma.hpp:20
template <typename T>
inline fvar <T >lgamma (const fvar <T >&x );

// ./stan/math/prim/fun/lmgamma.hpp:54
template <typename T,
    require_arithmetic_t<T>*>
inline return_type_t <T >lmgamma (int k , T x );

// ./stan/math/prim/fun/lmgamma.hpp:72
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto lmgamma (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/lmgamma.hpp:13
template <typename T>
inline fvar <return_type_t <T ,int >>lmgamma (int x1 , const fvar <T >&x2 );

// ./stan/math/prim/fun/lmultiply.hpp:24
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline return_type_t <T1 ,T2 >lmultiply (const T1 a , const T2 b );

// ./stan/math/prim/fun/lmultiply.hpp:44
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto lmultiply (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/lmultiply.hpp:12
template <typename T>
inline fvar <T >lmultiply (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/lmultiply.hpp:19
template <typename T>
inline fvar <T >lmultiply (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/lmultiply.hpp:25
template <typename T>
inline fvar <T >lmultiply (const fvar <T >&x1 , double x2 );

// ./stan/math/prim/fun/log10.hpp:37
template <typename Container,
    require_not_var_matrix_t<Container>*>
inline auto log10 (const Container &x );

// ./stan/math/prim/fun/log10.hpp:53
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto log10 (const Container &x );

// ./stan/math/fwd/fun/log10.hpp:14
template <typename T>
inline fvar <T >log10 (const fvar <T >&x );

// ./stan/math/fwd/fun/log10.hpp:32
template <typename T>
inline std ::complex <fvar <T >>log10 (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/log1m_exp.hpp:22
template <typename T>
inline fvar <T >log1m_exp (const fvar <T >&x );

// ./stan/math/prim/fun/log1p_exp.hpp:45

inline double log1p_exp (double a );

// ./stan/math/prim/fun/log1p_exp.hpp:75
template <typename T,
    require_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto log1p_exp (const T &x );

// ./stan/math/prim/fun/log1m_inv_logit.hpp:36

inline double log1m_inv_logit (double u );

// ./stan/math/prim/fun/log1m_inv_logit.hpp:51

inline double log1m_inv_logit (int u );

// ./stan/math/prim/fun/log1m_inv_logit.hpp:83
template <typename T,
    require_not_var_matrix_t<T>*>
inline typename apply_scalar_unary <log1m_inv_logit_fun ,T >::return_t log1m_inv_logit (const T &x );

// ./stan/math/fwd/fun/log1m_inv_logit.hpp:20
template <typename T>
inline fvar <T >log1m_inv_logit (const fvar <T >&x );

// ./stan/math/fwd/fun/log1p.hpp:11
template <typename T>
inline fvar <T >log1p (const fvar <T >&x );

// ./stan/math/fwd/fun/log1p_exp.hpp:12
template <typename T>
inline fvar <T >log1p_exp (const fvar <T >&x );

// ./stan/math/prim/fun/log2.hpp:17

inline constexpr double log2 ();

// ./stan/math/prim/fun/log2.hpp:47
template <typename T,
    require_not_var_matrix_t<T>*>
inline auto log2 (const T &x );

// ./stan/math/fwd/fun/log2.hpp:19
template <typename T>
inline fvar <T >log2 (const fvar <T >&x );

// ./stan/math/fwd/fun/log_determinant.hpp:15
template <typename EigMat,
    require_eigen_vt<is_fvar,
    EigMat>*>
inline value_type_t <EigMat >log_determinant (const EigMat &m );

// ./stan/math/fwd/fun/log_diff_exp.hpp:13
template <typename T>
inline fvar <T >log_diff_exp (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/log_diff_exp.hpp:26
template <typename T1,
    typename T2,
    require_arithmetic_t<T1>*>
inline fvar <T2 >log_diff_exp (const T1 &x1 , const fvar <T2 >&x2 );

// ./stan/math/fwd/fun/log_diff_exp.hpp:37
template <typename T1,
    typename T2,
    require_arithmetic_t<T2>*>
inline fvar <T1 >log_diff_exp (const fvar <T1 >&x1 , const T2 &x2 );

// ./stan/math/prim/fun/log_falling_factorial.hpp:54
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline return_type_t <T1 ,T2 >log_falling_factorial (const T1 x , const T2 n );

// ./stan/math/prim/fun/log_falling_factorial.hpp:74
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto log_falling_factorial (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/log_falling_factorial.hpp:13
template <typename T>
inline fvar <T >log_falling_factorial (const fvar <T >&x , const fvar <T >&n );

// ./stan/math/fwd/fun/log_falling_factorial.hpp:20
template <typename T>
inline fvar <T >log_falling_factorial (double x , const fvar <T >&n );

// ./stan/math/fwd/fun/log_falling_factorial.hpp:26
template <typename T>
inline fvar <T >log_falling_factorial (const fvar <T >&x , double n );

// ./stan/math/prim/fun/log_inv_logit.hpp:34

inline double log_inv_logit (double u );

// ./stan/math/prim/fun/log_inv_logit.hpp:49

inline double log_inv_logit (int u );

// ./stan/math/prim/fun/log_inv_logit.hpp:81
template <typename T,
    require_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto log_inv_logit (const T &x );

// ./stan/math/fwd/fun/log_inv_logit.hpp:13
template <typename T>
inline fvar <T >log_inv_logit (const fvar <T >&x );

// ./stan/math/prim/fun/log_inv_logit_diff.hpp:36
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline return_type_t <T1 ,T2 >log_inv_logit_diff (const T1 &x , const T2 &y );

// ./stan/math/prim/fun/log_inv_logit_diff.hpp:54
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto log_inv_logit_diff (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/log_inv_logit_diff.hpp:37
template <typename T>
inline fvar <T >log_inv_logit_diff (const fvar <T >&x , const fvar <T >&y );

// ./stan/math/fwd/fun/log_inv_logit_diff.hpp:45
template <typename T>
inline fvar <T >log_inv_logit_diff (const fvar <T >&x , double y );

// ./stan/math/fwd/fun/log_inv_logit_diff.hpp:51
template <typename T>
inline fvar <T >log_inv_logit_diff (double x , const fvar <T >&y );

// ./stan/math/prim/fun/log_sum_exp.hpp:51
template <typename T1,
    typename T2,
    require_all_not_st_var<T1,
    T2>*>
inline return_type_t <T1 ,T2 >log_sum_exp (const T2 &a , const T1 &b );

// ./stan/math/prim/fun/log_sum_exp.hpp:81
template <typename T,
    require_container_st<std::is_arithmetic,
    T>*>
inline auto log_sum_exp (const T &x );

// ./stan/math/prim/fun/log_sum_exp.hpp:106
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto log_sum_exp (const T1 &a , const T2 &b );

// ./stan/math/prim/functor/partials_propagator.hpp:82
template <std::size_tI,
    class ...Types>
inline constexpr auto &edge (internal ::partials_propagator <Types ...>&x )noexcept ;

// ./stan/math/prim/functor/partials_propagator.hpp:94
template <std::size_tI,
    class ...Types>
inline constexpr auto &partials (internal ::partials_propagator <Types ...>&x )noexcept ;

// ./stan/math/prim/functor/partials_propagator.hpp:106
template <std::size_tI,
    class ...Types>
inline constexpr auto &partials_vec (internal ::partials_propagator <Types ...>&x )noexcept ;

// ./stan/math/prim/functor/partials_propagator.hpp:118
template <typename ...Ops>
inline auto make_partials_propagator (Ops &&...ops );

// ./stan/math/prim/fun/log_mix.hpp:39
template <typename T_theta,
    typename T_lambda1,
    typename T_lambda2,
    require_all_arithmetic_t<T_theta,
    T_lambda1,
    T_lambda2>*>
inline double log_mix (T_theta theta , T_lambda1 lambda1 , T_lambda2 lambda2 );

// ./stan/math/prim/fun/log_mix.hpp:75
template <typename T_theta,
    typename T_lam,
    require_any_vector_t<T_theta,
    T_lam>*>
return_type_t <T_theta ,T_lam >log_mix (const T_theta &theta , const T_lam &lambda );

// ./stan/math/fwd/fun/log_mix.hpp:27
template <typename T_theta,
    typename T_lambda1,
    typename T_lambda2,
    int N>
inline void log_mix_partial_helper (const T_theta &theta , const T_lambda1 &lambda1 , const T_lambda2 &lambda2 , promote_args_t <T_theta , T_lambda1 , T_lambda2 >(&partials_array )[N ]);

// ./stan/math/fwd/fun/log_mix.hpp:97
template <typename T>
inline fvar <T >log_mix (const fvar <T >&theta , const fvar <T >&lambda1 , const fvar <T >&lambda2 );

// ./stan/math/fwd/fun/log_mix.hpp:117
template <typename T,
    typename P,
    require_all_arithmetic_t<P>*>
inline fvar <T >log_mix (const fvar <T >&theta , const fvar <T >&lambda1 , P lambda2 );

// ./stan/math/fwd/fun/log_mix.hpp:135
template <typename T,
    typename P,
    require_all_arithmetic_t<P>*>
inline fvar <T >log_mix (const fvar <T >&theta , P lambda1 , const fvar <T >&lambda2 );

// ./stan/math/fwd/fun/log_mix.hpp:153
template <typename T,
    typename P,
    require_all_arithmetic_t<P>*>
inline fvar <T >log_mix (P theta , const fvar <T >&lambda1 , const fvar <T >&lambda2 );

// ./stan/math/fwd/fun/log_mix.hpp:171
template <typename T,
    typename P1,
    typename P2,
    require_all_arithmetic_t<P1,
    P2>*>
inline fvar <T >log_mix (const fvar <T >&theta , P1 lambda1 , P2 lambda2 );

// ./stan/math/fwd/fun/log_mix.hpp:187
template <typename T,
    typename P1,
    typename P2,
    require_all_arithmetic_t<P1,
    P2>*>
inline fvar <T >log_mix (P1 theta , const fvar <T >&lambda1 , P2 lambda2 );

// ./stan/math/fwd/fun/log_mix.hpp:203
template <typename T,
    typename P1,
    typename P2,
    require_all_arithmetic_t<P1,
    P2>*>
inline fvar <T >log_mix (P1 theta , P2 lambda1 , const fvar <T >&lambda2 );

// ./stan/math/prim/fun/log_rising_factorial.hpp:52
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline return_type_t <T1 ,T2 >log_rising_factorial (const T1 &x , const T2 &n );

// ./stan/math/prim/fun/log_rising_factorial.hpp:72
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto log_rising_factorial (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/log_rising_factorial.hpp:12
template <typename T>
inline fvar <T >log_rising_factorial (const fvar <T >&x , const fvar <T >&n );

// ./stan/math/fwd/fun/log_rising_factorial.hpp:19
template <typename T>
inline fvar <T >log_rising_factorial (const fvar <T >&x , double n );

// ./stan/math/fwd/fun/log_rising_factorial.hpp:25
template <typename T>
inline fvar <T >log_rising_factorial (double x , const fvar <T >&n );

// ./stan/math/prim/fun/softmax.hpp:45
template <typename ColVec,
    require_eigen_col_vector_vt<std::is_arithmetic,
    ColVec>*>
inline plain_type_t <ColVec >softmax (const ColVec &v );

// ./stan/math/fwd/fun/softmax.hpp:14
template <typename ColVec,
    require_eigen_col_vector_vt<is_fvar,
    ColVec>*>
inline auto softmax (const ColVec &alpha );

// ./stan/math/prim/fun/log_softmax.hpp:42
template <typename Container,
    require_st_arithmetic<Container>*>
inline auto log_softmax (const Container &x );

// ./stan/math/fwd/fun/log_softmax.hpp:22
template <typename T,
    require_vector_st<is_fvar,
    T>*>
inline auto log_softmax (const T &x );

// ./stan/math/fwd/fun/log_sum_exp.hpp:17
template <typename T>
inline fvar <T >log_sum_exp (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/log_sum_exp.hpp:25
template <typename T>
inline fvar <T >log_sum_exp (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/log_sum_exp.hpp:34
template <typename T>
inline fvar <T >log_sum_exp (const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/fun/log_sum_exp.hpp:54
template <typename T,
    require_container_st<is_fvar,
    T>*>
inline auto log_sum_exp (const T &x );

// ./stan/math/prim/fun/logit.hpp:46

inline double logit (double u );

// ./stan/math/prim/fun/logit.hpp:57

inline double logit (int u );

// ./stan/math/prim/fun/logit.hpp:86
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto logit (const Container &x );

// ./stan/math/prim/fun/logit.hpp:106
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto logit (const Container &x );

// ./stan/math/fwd/fun/logit.hpp:13
template <typename T>
inline fvar <T >logit (const fvar <T >&x );

// ./stan/math/prim/fun/mdivide_left.hpp:23
template <typename T1,
    typename T2,
    require_all_eigen_vt<std::is_arithmetic,
    T1,
    T2>*>
inline Eigen ::Matrix <return_type_t <T1 ,T2 >,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left (const T1 &A , const T2 &b );

// ./stan/math/fwd/fun/mdivide_left.hpp:16
template <typename T1,
    typename T2,
    require_all_eigen_vt<is_fvar,
    T1,
    T2>*>
inline Eigen ::Matrix <value_type_t <T1 >,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left (const T1 &A , const T2 &b );

// ./stan/math/fwd/fun/mdivide_left.hpp:57
template <typename T1,
    typename T2,
    require_eigen_vt<std::is_arithmetic,
    T1>*>
inline Eigen ::Matrix <value_type_t <T2 >,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left (const T1 &A , const T2 &b );

// ./stan/math/fwd/fun/mdivide_left.hpp:74
template <typename T1,
    typename T2,
    require_eigen_vt<is_fvar,
    T1>*>
inline Eigen ::Matrix <value_type_t <T1 >,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left (const T1 &A , const T2 &b );

// ./stan/math/fwd/fun/mdivide_left_ldlt.hpp:23
template <typename T,
    typename EigMat,
    require_eigen_vt<std::is_arithmetic,
    T>*>
inline Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,EigMat ::ColsAtCompileTime >mdivide_left_ldlt (LDLT_factor <T >&A , const EigMat &b );

// ./stan/math/fwd/fun/mdivide_left_tri_low.hpp:15
template <typename T1,
    typename T2,
    require_all_eigen_vt<is_fvar,
    T1,
    T2>*>
inline Eigen ::Matrix <value_type_t <T1 >,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left_tri_low (const T1 &A , const T2 &b );

// ./stan/math/fwd/fun/mdivide_left_tri_low.hpp:52
template <typename T1,
    typename T2,
    require_eigen_t<T1>*>
inline Eigen ::Matrix <value_type_t <T2 >,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left_tri_low (const T1 &A , const T2 &b );

// ./stan/math/fwd/fun/mdivide_left_tri_low.hpp:81
template <typename T1,
    typename T2,
    require_eigen_vt<is_fvar,
    T1>*>
inline Eigen ::Matrix <value_type_t <T1 >,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left_tri_low (const T1 &A , const T2 &b );

// ./stan/math/prim/fun/mdivide_right.hpp:23
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline Eigen ::Matrix <return_type_t <EigMat1 ,EigMat2 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_right (const EigMat1 &b , const EigMat2 &A );

// ./stan/math/fwd/fun/mdivide_right.hpp:16
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_vt<is_fvar,
    EigMat1,
    EigMat2>*>
inline Eigen ::Matrix <value_type_t <EigMat1 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_right (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/fwd/fun/mdivide_right.hpp:62
template <typename EigMat1,
    typename EigMat2,
    require_eigen_vt<is_fvar,
    EigMat1>*>
inline Eigen ::Matrix <value_type_t <EigMat1 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_right (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/fwd/fun/mdivide_right.hpp:92
template <typename EigMat1,
    typename EigMat2,
    require_eigen_vt<std::is_arithmetic,
    EigMat1>*>
inline Eigen ::Matrix <value_type_t <EigMat2 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_right (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/fwd/fun/mdivide_right_tri_low.hpp:14
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_vt<is_fvar,
    EigMat1,
    EigMat2>*>
inline Eigen ::Matrix <value_type_t <EigMat1 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_right_tri_low (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/fwd/fun/mdivide_right_tri_low.hpp:61
template <typename EigMat1,
    typename EigMat2,
    require_eigen_vt<is_fvar,
    EigMat1>*>
inline Eigen ::Matrix <value_type_t <EigMat1 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_right_tri_low (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/fwd/fun/mdivide_right_tri_low.hpp:93
template <typename EigMat1,
    typename EigMat2,
    require_eigen_vt<std::is_arithmetic,
    EigMat1>*>
inline Eigen ::Matrix <value_type_t <EigMat2 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_right_tri_low (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/prim/fun/modified_bessel_first_kind.hpp:37
template <typename T2,
    require_arithmetic_t<T2>*>
inline T2 modified_bessel_first_kind (int v , const T2 z );

// ./stan/math/prim/fun/modified_bessel_first_kind.hpp:49

inline double modified_bessel_first_kind (int v , int z );

// ./stan/math/prim/fun/modified_bessel_first_kind.hpp:65
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto modified_bessel_first_kind (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/modified_bessel_first_kind.hpp:11
template <typename T>
inline fvar <T >modified_bessel_first_kind (int v , const fvar <T >&z );

// ./stan/math/prim/fun/modified_bessel_second_kind.hpp:42
template <typename T2,
    require_arithmetic_t<T2>*>
inline T2 modified_bessel_second_kind (int v , const T2 z );

// ./stan/math/prim/fun/modified_bessel_second_kind.hpp:57
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto modified_bessel_second_kind (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/modified_bessel_second_kind.hpp:11
template <typename T>
inline fvar <T >modified_bessel_second_kind (int v , const fvar <T >&z );

// ./stan/math/fwd/fun/multiply_log.hpp:12
template <typename T>
inline fvar <T >multiply_log (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/multiply_log.hpp:19
template <typename T>
inline fvar <T >multiply_log (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/multiply_log.hpp:25
template <typename T>
inline fvar <T >multiply_log (const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/fun/multiply_lower_tri_self_transpose.hpp:12
template <typename EigMat,
    require_eigen_vt<is_fvar,
    EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,EigMat ::RowsAtCompileTime ,EigMat ::RowsAtCompileTime >multiply_lower_tri_self_transpose (const EigMat &m );

// ./stan/math/prim/fun/norm.hpp:17
template <typename V>
inline V norm (const std ::complex <V >&z );

// ./stan/math/fwd/fun/norm.hpp:18
template <typename T>
inline fvar <T >norm (const std ::complex <fvar <T >>&z );

// ./stan/math/prim/fun/norm1.hpp:19
template <typename T,
    require_eigen_vt<std::is_arithmetic,
    T>*>
inline double norm1 (const T &v );

// ./stan/math/prim/fun/norm1.hpp:33
template <typename Container,
    require_std_vector_t<Container>*>
inline auto norm1 (const Container &x );

// ./stan/math/fwd/fun/norm1.hpp:23
template <typename Container,
    require_eigen_vt<is_fvar,
    Container>*>
inline auto norm1 (const Container &x );

// ./stan/math/prim/fun/norm2.hpp:19
template <typename T,
    require_eigen_vt<std::is_arithmetic,
    T>*>
inline double norm2 (const T &v );

// ./stan/math/prim/fun/norm2.hpp:33
template <typename Container,
    require_std_vector_t<Container>*>
inline auto norm2 (const Container &x );

// ./stan/math/fwd/fun/norm2.hpp:21
template <typename Container,
    require_eigen_vt<is_fvar,
    Container>*>
inline auto norm2 (const Container &x );

// ./stan/math/prim/fun/owens_t.hpp:58

inline double owens_t (double h , double a );

// ./stan/math/prim/fun/owens_t.hpp:70
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto owens_t (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/owens_t.hpp:24
template <typename T>
inline fvar <T >owens_t (const fvar <T >&x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/owens_t.hpp:46
template <typename T>
inline fvar <T >owens_t (double x1 , const fvar <T >&x2 );

// ./stan/math/fwd/fun/owens_t.hpp:65
template <typename T>
inline fvar <T >owens_t (const fvar <T >&x1 , double x2 );

// ./stan/math/fwd/fun/Phi.hpp:13
template <typename T>
inline fvar <T >Phi (const fvar <T >&x );

// ./stan/math/fwd/fun/Phi_approx.hpp:22
template <typename T>
inline fvar <T >Phi_approx (const fvar <T >&x );

// ./stan/math/fwd/fun/polar.hpp:21
template <typename T>
inline std ::complex <fvar <T >>polar (const fvar <T >&r , const fvar <T >&theta );

// ./stan/math/fwd/fun/polar.hpp:35
template <typename T,
    typename U>
inline std ::complex <fvar <T >>polar (const fvar <T >&r , U theta );

// ./stan/math/fwd/fun/polar.hpp:49
template <typename T,
    typename U>
inline std ::complex <fvar <T >>polar (U r , const fvar <T >&theta );

// ./stan/math/fwd/fun/primitive_value.hpp:20
template <typename T>
inline double primitive_value (const fvar <T >&v );

// ./stan/math/prim/fun/proj.hpp:18
template <typename V>
inline std ::complex <V >proj (const std ::complex <V >&z );

// ./stan/math/fwd/fun/proj.hpp:19
template <typename T>
inline std ::complex <fvar <T >>proj (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/quad_form.hpp:27
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline promote_scalar_t <return_type_t <EigMat1 ,EigMat2 >,EigMat2 >quad_form (const EigMat1 &A , const EigMat2 &B );

// ./stan/math/fwd/fun/quad_form.hpp:51
template <typename EigMat,
    typename ColVec,
    require_eigen_t<EigMat>*>
inline return_type_t <EigMat ,ColVec >quad_form (const EigMat &A , const ColVec &B );

// ./stan/math/fwd/fun/quad_form_sym.hpp:26
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline promote_scalar_t <return_type_t <EigMat1 ,EigMat2 >,EigMat2 >quad_form_sym (const EigMat1 &A , const EigMat2 &B );

// ./stan/math/fwd/fun/quad_form_sym.hpp:53
template <typename EigMat,
    typename ColVec,
    require_eigen_t<EigMat>*>
inline return_type_t <EigMat ,ColVec >quad_form_sym (const EigMat &A , const ColVec &B );

// ./stan/math/prim/fun/rising_factorial.hpp:62
template <typename T,
    require_arithmetic_t<T>*>
inline return_type_t <T >rising_factorial (const T &x , int n );

// ./stan/math/prim/fun/rising_factorial.hpp:80
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto rising_factorial (const T1 &a , const T2 &b );

// ./stan/math/fwd/fun/rising_factorial.hpp:25
template <typename T>
inline fvar <T >rising_factorial (const fvar <T >&x , int n );

// ./stan/math/prim/fun/round.hpp:35
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto round (const Container &x );

// ./stan/math/prim/fun/round.hpp:51
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto round (const Container &x );

// ./stan/math/fwd/fun/round.hpp:23
template <typename T>
inline fvar <T >round (const fvar <T >&x );

// ./stan/math/fwd/fun/sin.hpp:13
template <typename T>
inline fvar <T >sin (const fvar <T >&x );

// ./stan/math/fwd/fun/sin.hpp:27
template <typename T>
inline std ::complex <fvar <T >>sin (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/sinh.hpp:12
template <typename T>
inline fvar <T >sinh (const fvar <T >&x );

// ./stan/math/fwd/fun/sinh.hpp:26
template <typename T>
inline std ::complex <fvar <T >>sinh (const std ::complex <fvar <T >>&z );

// ./stan/math/prim/fun/tanh.hpp:38
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto tanh (const Container &x );

// ./stan/math/prim/fun/tanh.hpp:55
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto tanh (const Container &x );

// ./stan/math/prim/fun/tan.hpp:38
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto tan (const Container &x );

// ./stan/math/prim/fun/tan.hpp:55
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto tan (const Container &x );

// ./stan/math/fwd/fun/tan.hpp:13
template <typename T>
inline fvar <T >tan (const fvar <T >&x );

// ./stan/math/fwd/fun/tan.hpp:26
template <typename T>
inline std ::complex <fvar <T >>tan (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/tanh.hpp:14
template <typename T>
inline fvar <T >tanh (const fvar <T >&x );

// ./stan/math/fwd/fun/tanh.hpp:28
template <typename T>
inline std ::complex <fvar <T >>tanh (const std ::complex <fvar <T >>&z );

// ./stan/math/fwd/fun/tgamma.hpp:20
template <typename T>
inline fvar <T >tgamma (const fvar <T >&x );

// ./stan/math/prim/fun/trace.hpp:20
template <typename T,
    require_eigen_t<T>*>
inline value_type_t <T >trace (const T &m );

// ./stan/math/fwd/fun/trace_quad_form.hpp:15
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline return_type_t <EigMat1 ,EigMat2 >trace_quad_form (const EigMat1 &A , const EigMat2 &B );

// ./stan/math/fwd/fun/trigamma.hpp:20
template <typename T>
inline fvar <T >trigamma (const fvar <T >&u );

// ./stan/math/prim/fun/trunc.hpp:40
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto trunc (const T &x );

// ./stan/math/fwd/fun/trunc.hpp:19
template <typename T>
inline fvar <T >trunc (const fvar <T >&x );

// ./stan/math/fwd/fun/value_of_rec.hpp:19
template <typename T>
inline double value_of_rec (const fvar <T >&v );

// ./stan/math/fwd/functor/gradient.hpp:39
template <typename T,
    typename F>
void gradient (const F &f , const Eigen ::Matrix <T , Eigen ::Dynamic , 1>&x , T &fx , Eigen ::Matrix <T , Eigen ::Dynamic , 1>&grad_fx );

// ./stan/math/prim/fun/to_array_1d.hpp:14
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline std ::vector <value_type_t <EigMat >>to_array_1d (const EigMat &matrix );

// ./stan/math/prim/fun/to_array_1d.hpp:26
template <typename T>
inline std ::vector <T >to_array_1d (const std ::vector <T >&x );

// ./stan/math/prim/fun/to_array_1d.hpp:32
template <typename T>
inline std ::vector <typename scalar_type <T >::type >to_array_1d (const std ::vector <std ::vector <T >>&x );

// ./stan/math/prim/fun/serializer.hpp:271
template <typename T>
deserializer <T >to_deserializer (const std ::vector <T >&vals );

// ./stan/math/prim/fun/serializer.hpp:283
template <typename T,
    require_eigen_vector_t<T>*>
deserializer <scalar_type_t <T >>to_deserializer (const T &vals );

// ./stan/math/prim/fun/serializer.hpp:288
template <typename T>
deserializer <T >to_deserializer (const std ::complex <T >&vals );

// ./stan/math/prim/fun/serializer.hpp:293
template <typename U>
void serialize_helper (serializer <U >&s );

// ./stan/math/prim/fun/serializer.hpp:296
template <typename U,
    typename T,
    typename ...Ts>
void serialize_helper (serializer <U >&s , const T &x , const Ts ...xs );

// ./stan/math/prim/fun/serializer.hpp:311
template <typename U,
    typename ...Ts>
std ::vector <U >serialize (const Ts ...xs );

// ./stan/math/prim/fun/serializer.hpp:325
template <typename T>
std ::vector <real_return_t <T >>serialize_return (const T &x );

// ./stan/math/prim/fun/serializer.hpp:338
template <typename ...Ts>
Eigen ::VectorXd serialize_args (const Ts ...xs );

// ./stan/math/fwd/functor/finite_diff.hpp:67
template <typename F,
    typename ...TArgs,
    require_any_st_fvar<TArgs...>*>
inline auto finite_diff (const F &func , const TArgs &...args );

// ./stan/math/fwd/functor/finite_diff.hpp:111
template <typename F,
    typename ...TArgs,
    require_all_not_st_fvar<TArgs...>*>
inline auto finite_diff (const F &func , const TArgs &...args );

// ./stan/math/fwd/functor/hessian.hpp:40
template <typename T,
    typename F>
void hessian (const F &f , const Eigen ::Matrix <T , Eigen ::Dynamic , 1>&x , T &fx , Eigen ::Matrix <T , Eigen ::Dynamic , 1>&grad , Eigen ::Matrix <T , Eigen ::Dynamic , Eigen ::Dynamic >&H );

// ./stan/math/prim/functor/integrate_1d.hpp:51
template <typename F>
inline double integrate (const F &f , double a , double b , double relative_tolerance );

// ./stan/math/prim/functor/integrate_1d.hpp:173
template <typename F,
    typename ...Args,
    require_all_st_arithmetic<Args...>*>
inline double integrate_1d_impl (const F &f , double a , double b , double relative_tolerance , std ::ostream *msgs , const Args &...args );

// ./stan/math/prim/functor/integrate_1d.hpp:238
template <typename F>
inline double integrate_1d (const F &f , double a , double b , const std ::vector <double >&theta , const std ::vector <double >&x_r , const std ::vector <int >&x_i , std ::ostream *msgs , const double relative_tolerance =std ::sqrt (EPSILON ));

// ./stan/math/prim/functor/apply.hpp:51
template <class F,
    class Tuple,
    typename ...PreArgs>
constexpr decltype (auto )apply (F &&f , Tuple &&t , PreArgs &&...pre_args );

// ./stan/math/fwd/functor/integrate_1d.hpp:29
template <typename F,
    typename T_a,
    typename T_b,
    typename ...Args,
    require_any_st_fvar<T_a,
    T_b,
    Args...>*>
inline return_type_t <T_a ,T_b ,Args ...>integrate_1d_impl (const F &f , const T_a &a , const T_b &b , double relative_tolerance , std ::ostream *msgs , const Args &...args );

// ./stan/math/fwd/functor/integrate_1d.hpp:89
template <typename F,
    typename T_a,
    typename T_b,
    typename T_theta,
    require_any_fvar_t<T_a,
    T_b,
    T_theta>*>
inline return_type_t <T_a ,T_b ,T_theta >integrate_1d (const F &f , const T_a &a , const T_b &b , const std ::vector <T_theta >&theta , const std ::vector <double >&x_r , const std ::vector <int >&x_i , std ::ostream *msgs , const double relative_tolerance );

// ./stan/math/fwd/functor/jacobian.hpp:10
template <typename T,
    typename F>
void jacobian (const F &f , const Eigen ::Matrix <T , Eigen ::Dynamic , 1>&x , Eigen ::Matrix <T , Eigen ::Dynamic , 1>&fx , Eigen ::Matrix <T , Eigen ::Dynamic , Eigen ::Dynamic >&J );

// ./stan/math/prim/fun/elt_multiply.hpp:23
template <typename Mat1,
    typename Mat2,
    require_all_eigen_t<Mat1,
    Mat2>*>
auto elt_multiply (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/prim/fun/elt_multiply.hpp:42
template <typename Scalar1,
    typename Scalar2,
    require_all_stan_scalar_t<Scalar1,
    Scalar2>*>
auto elt_multiply (const Scalar1 &a , const Scalar2 &b );

// ./stan/math/prim/fun/elt_multiply.hpp:60
template <typename T1,
    typename T2,
    require_any_matrix_t<T1,
    T2>*>
inline auto elt_multiply (const T1 &A , const T2 &B );

// ./stan/math/prim/functor/for_each.hpp:65
template <typename F,
    typename T>
constexpr inline auto for_each (F &&f , T &&t );

// ./stan/math/prim/functor/for_each.hpp:81
template <typename F,
    typename T1,
    typename T2>
constexpr inline auto for_each (F &&f , T1 &&t1 , T2 &&t2 );

// ./stan/math/prim/functor/for_each.hpp:103
template <typename F,
    typename T1,
    typename T2,
    typename T3>
constexpr inline auto for_each (F &&f , T1 &&t1 , T2 &&t2 , T3 &&t3 );

// ./stan/math/prim/prob/std_normal_log_qf.hpp:26

inline double std_normal_log_qf (double log_p );

// ./stan/math/prim/prob/std_normal_log_qf.hpp:149
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto std_normal_log_qf (const T &x );

// ./stan/math/fwd/prob/std_normal_log_qf.hpp:16
template <typename T>
inline fvar <T >std_normal_log_qf (const fvar <T >&p );

// ./stan/math/prim/constraint/corr_constrain.hpp:26
template <typename T>
inline plain_type_t <T >corr_constrain (const T &x );

// ./stan/math/prim/constraint/corr_constrain.hpp:45
template <typename T_x,
    typename T_lp>
inline auto corr_constrain (const T_x &x , T_lp &lp );

// ./stan/math/prim/constraint/corr_constrain.hpp:67
template <bool Jacobian,
    typename T_x,
    typename T_lp>
inline auto corr_constrain (const T_x &x , T_lp &lp );

// ./stan/math/rev/constraint/cholesky_corr_constrain.hpp:27
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >cholesky_corr_constrain (const T &y , int K );

// ./stan/math/rev/constraint/cholesky_corr_constrain.hpp:91
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >cholesky_corr_constrain (const T &y , int K , scalar_type_t <T >&lp );

// ./stan/math/rev/constraint/cholesky_factor_constrain.hpp:29
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >cholesky_factor_constrain (const T &x , int M , int N );

// ./stan/math/rev/constraint/cholesky_factor_constrain.hpp:88
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >cholesky_factor_constrain (const T &x , int M , int N , scalar_type_t <T >&lp );

// ./stan/math/prim/fun/read_corr_L.hpp:36
template <typename T,
    require_eigen_vector_t<T>*>
Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >read_corr_L (const T &CPCs , // on (-1, 1)size_t K );

// ./stan/math/prim/fun/read_corr_L.hpp:100
template <typename T,
    require_eigen_vector_t<T>*>
Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >read_corr_L (const T &CPCs , size_t K , value_type_t <T >&log_prob );

// ./stan/math/prim/fun/multiply_lower_tri_self_transpose.hpp:18
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
inline matrix_d multiply_lower_tri_self_transpose (const EigMat &L );

// ./stan/math/prim/fun/read_corr_matrix.hpp:25
template <typename T_CPCs,
    require_eigen_vector_t<T_CPCs>*>
Eigen ::Matrix <value_type_t <T_CPCs >,Eigen ::Dynamic ,Eigen ::Dynamic >read_corr_matrix (const T_CPCs &CPCs , size_t K );

// ./stan/math/prim/fun/read_corr_matrix.hpp:54
template <typename T_CPCs,
    require_eigen_vector_t<T_CPCs>*>
Eigen ::Matrix <value_type_t <T_CPCs >,Eigen ::Dynamic ,Eigen ::Dynamic >read_corr_matrix (const T_CPCs &CPCs , size_t K , value_type_t <T_CPCs >&log_prob );

// ./stan/math/rev/constraint/corr_matrix_constrain.hpp:40
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >corr_matrix_constrain (const T &x , Eigen ::Index k );

// ./stan/math/rev/constraint/corr_matrix_constrain.hpp:68
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >corr_matrix_constrain (const T &x , Eigen ::Index k , scalar_type_t <T >&lp );

// ./stan/math/rev/fun/dot_product.hpp:34
template <typename T1,
    typename T2,
    require_all_vector_t<T1,
    T2>*
;

// ./stan/math/rev/fun/dot_self.hpp:23
template <typename T,
    require_eigen_vector_vt<is_var,
    T>*>
inline var dot_self (const T &v );

// ./stan/math/rev/fun/dot_self.hpp:47
template <typename T,
    require_var_matrix_t<T>*>
inline var dot_self (const T &v );

// ./stan/math/rev/fun/multiply_lower_tri_self_transpose.hpp:16
template <typename T,
    require_rev_matrix_t<T>*>
inline auto multiply_lower_tri_self_transpose (const T &L );

// ./stan/math/rev/constraint/cov_matrix_constrain.hpp:31
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >cov_matrix_constrain (const T &x , Eigen ::Index K );

// ./stan/math/rev/constraint/cov_matrix_constrain.hpp:76
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >cov_matrix_constrain (const T &x , Eigen ::Index K , scalar_type_t <T >&lp );

// ./stan/math/prim/constraint/positive_constrain.hpp:22
template <typename T>
inline auto positive_constrain (const T &x );

// ./stan/math/prim/constraint/positive_constrain.hpp:43
template <typename T,
    typename S>
inline auto positive_constrain (const T &x , S &lp );

// ./stan/math/prim/constraint/positive_constrain.hpp:64
template <bool Jacobian,
    typename T,
    require_not_std_vector_t<T>*>
inline auto positive_constrain (const T &x , return_type_t <T >&lp );

// ./stan/math/prim/constraint/positive_constrain.hpp:89
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto positive_constrain (const T &x , return_type_t <T >&lp );

// ./stan/math/prim/fun/read_cov_L.hpp:27
template <typename T_CPCs,
    typename T_sds,
    require_all_eigen_vector_t<T_CPCs,
    T_sds>*>
Eigen ::Matrix <value_type_t <T_CPCs >,Eigen ::Dynamic ,Eigen ::Dynamic >read_cov_L (const T_CPCs &CPCs , const T_sds &sds , value_type_t <T_CPCs >&log_prob );

// ./stan/math/prim/fun/read_cov_matrix.hpp:25
template <typename T_CPCs,
    typename T_sds,
    require_all_eigen_vector_t<T_CPCs,
    T_sds>*>
Eigen ::Matrix <value_type_t <T_CPCs >,Eigen ::Dynamic ,Eigen ::Dynamic >read_cov_matrix (const T_CPCs &CPCs , const T_sds &sds , value_type_t <T_CPCs >&log_prob );

// ./stan/math/prim/fun/read_cov_matrix.hpp:46
template <typename T_CPCs,
    typename T_sds,
    require_all_eigen_vector_t<T_CPCs,
    T_sds>*>
Eigen ::Matrix <value_type_t <T_CPCs >,Eigen ::Dynamic ,Eigen ::Dynamic >read_cov_matrix (const T_CPCs &CPCs , const T_sds &sds );

// ./stan/math/rev/constraint/cov_matrix_constrain_lkj.hpp:34
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >cov_matrix_constrain_lkj (const T &x , size_t k );

// ./stan/math/rev/constraint/cov_matrix_constrain_lkj.hpp:66
template <typename T,
    require_var_vector_t<T>*>
var_value <Eigen ::MatrixXd >cov_matrix_constrain_lkj (const T &x , size_t k , scalar_type_t <T >&lp );

// ./stan/math/prim/constraint/identity_constrain.hpp:21
template <bool Jacobian>
inline auto identity_constrain (T &&x , Types &&.../* args */);

// ./stan/math/rev/fun/to_var_value.hpp:20
template <typename T,
    require_eigen_vt<is_var,
    T>*>
inline var_value <Eigen ::Matrix <double ,T ::RowsAtCompileTime ,T ::ColsAtCompileTime >>to_var_value (const T &a );

// ./stan/math/rev/fun/to_var_value.hpp:37
template <typename T,
    require_var_t<T>*>
inline T to_var_value (T &&a );

// ./stan/math/rev/fun/to_var_value.hpp:49
template <typename T>
inline auto to_var_value (const std ::vector <T >&a );

// ./stan/math/rev/constraint/identity_constrain.hpp:22
template <typename T,
    typename ...Types,
    require_eigen_vt<std::is_arithmetic,
    T>*>
inline auto identity_constrain (T &&x , Types &&.../* args */);

// ./stan/math/rev/constraint/identity_constrain.hpp:29
template <typename T,
    typename ...Types,
    require_eigen_vt<is_var,
    T>*>
inline auto identity_constrain (T &&x , Types &&.../* args */);

// ./stan/math/rev/constraint/identity_constrain.hpp:35
template <typename T,
    typename ...Types,
    require_var_matrix_t<T>*>
inline auto identity_constrain (T &&x , Types &&.../* args */);

// ./stan/math/prim/constraint/identity_free.hpp:20
template <typename T,
    typename ...Types,
    require_all_not_var_matrix_t<T,
    Types...>*>
inline auto identity_free (T &&x , Types &&.../* args */);

// ./stan/math/prim/constraint/lb_constrain.hpp:35
template <typename T,
    typename L,
    require_all_stan_scalar_t<T,
    L>*>
inline auto lb_constrain (const T &x , const L &lb );

// ./stan/math/prim/constraint/lb_constrain.hpp:58
template <typename T,
    typename L,
    require_all_stan_scalar_t<T,
    L>*>
inline auto lb_constrain (const T &x , const L &lb , return_type_t <T , L >&lp );

// ./stan/math/prim/constraint/lb_constrain.hpp:79
template <typename T,
    typename L,
    require_eigen_t<T>*>
inline auto lb_constrain (T &&x , L &&lb );

// ./stan/math/prim/constraint/lb_constrain.hpp:97
template <typename T,
    typename L,
    require_eigen_t<T>*>
inline auto lb_constrain (const T &x , const L &lb , return_type_t <T , L >&lp );

// ./stan/math/prim/constraint/lb_constrain.hpp:115
template <typename T,
    typename L,
    require_all_eigen_t<T,
    L>*>
inline auto lb_constrain (T &&x , L &&lb );

// ./stan/math/prim/constraint/lb_constrain.hpp:134
template <typename T,
    typename L,
    require_all_eigen_t<T,
    L>*>
inline auto lb_constrain (const T &x , const L &lb , return_type_t <T , L >&lp );

// ./stan/math/prim/constraint/lb_constrain.hpp:152
template <typename T,
    typename L,
    require_not_std_vector_t<L>*>
inline auto lb_constrain (const std ::vector <T >&x , const L &lb );

// ./stan/math/prim/constraint/lb_constrain.hpp:172
template <typename T,
    typename L,
    require_not_std_vector_t<L>*>
inline auto lb_constrain (const std ::vector <T >&x , const L &lb , return_type_t <T , L >&lp );

// ./stan/math/prim/constraint/lb_constrain.hpp:192
template <typename T,
    typename L>
inline auto lb_constrain (const std ::vector <T >&x , const std ::vector <L >&lb );

// ./stan/math/prim/constraint/lb_constrain.hpp:213
template <typename T,
    typename L>
inline auto lb_constrain (const std ::vector <T >&x , const std ::vector <L >&lb , return_type_t <T , L >&lp );

// ./stan/math/prim/constraint/lb_constrain.hpp:242
template <bool Jacobian,
    typename T,
    typename L>
inline auto lb_constrain (const T &x , const L &lb , return_type_t <T , L >&lp );

// ./stan/math/rev/fun/sum.hpp:25
template <typename Alloc>
inline var sum (const std ::vector <var , Alloc >&m );

// ./stan/math/rev/fun/sum.hpp:46
template <typename T,
    require_rev_matrix_t<T>*>
inline var sum (T &&x );

// ./stan/math/rev/constraint/lb_constrain.hpp:40
template <typename T,
    typename L,
    require_all_stan_scalar_t<T,
    L>*>
inline auto lb_constrain (const T &x , const L &lb );

// ./stan/math/rev/constraint/lb_constrain.hpp:90
template <typename T,
    typename L,
    require_all_stan_scalar_t<T,
    L>*>
inline auto lb_constrain (const T &x , const L &lb , var &lp );

// ./stan/math/rev/constraint/lb_constrain.hpp:132
template <typename T,
    typename L,
    require_matrix_t<T>*>
inline auto lb_constrain (const T &x , const L &lb );

// ./stan/math/rev/constraint/lb_constrain.hpp:181
template <typename T,
    typename L,
    require_matrix_t<T>*>
inline auto lb_constrain (const T &x , const L &lb , return_type_t <T , L >&lp );

// ./stan/math/rev/constraint/lb_constrain.hpp:234
template <typename T,
    typename L,
    require_all_matrix_t<T,
    L>*>
inline auto lb_constrain (const T &x , const L &lb );

// ./stan/math/rev/constraint/lb_constrain.hpp:298
template <typename T,
    typename L,
    require_all_matrix_t<T,
    L>*>
inline auto lb_constrain (const T &x , const L &lb , return_type_t <T , L >&lp );

// ./stan/math/prim/fun/subtract.hpp:21
template <typename ScalarA,
    typename ScalarB,
    require_all_stan_scalar_t<ScalarA,
    ScalarB>*>
inline return_type_t <ScalarA ,ScalarB >subtract (const ScalarA &a , const ScalarB &b );

// ./stan/math/prim/fun/subtract.hpp:40
template <typename Mat1,
    typename Mat2,
    require_all_eigen_t<Mat1,
    Mat2>*>
inline auto subtract (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/prim/fun/subtract.hpp:58
template <typename Scal,
    typename Mat,
    require_stan_scalar_t<Scal>*>
inline auto subtract (const Scal c , const Mat &m );

// ./stan/math/prim/fun/subtract.hpp:75
template <typename Mat,
    typename Scal,
    require_eigen_t<Mat>*>
inline auto subtract (const Mat &m , const Scal c );

// ./stan/math/prim/constraint/ub_constrain.hpp:34
template <typename T,
    typename U,
    require_all_stan_scalar_t<T,
    U>*>
inline auto ub_constrain (const T &x , const U &ub );

// ./stan/math/prim/constraint/ub_constrain.hpp:65
template <typename T,
    typename U,
    require_all_stan_scalar_t<T,
    U>*>
inline auto ub_constrain (const T &x , const U &ub , std ::decay_t <return_type_t <T , U >>&lp );

// ./stan/math/prim/constraint/ub_constrain.hpp:87
template <typename T,
    typename U,
    require_eigen_t<T>*>
inline auto ub_constrain (const T &x , const U &ub );

// ./stan/math/prim/constraint/ub_constrain.hpp:105
template <typename T,
    typename U,
    require_eigen_t<T>*>
inline auto ub_constrain (const T &x , const U &ub , std ::decay_t <return_type_t <T , U >>&lp );

// ./stan/math/prim/constraint/ub_constrain.hpp:124
template <typename T,
    typename U,
    require_all_eigen_t<T,
    U>*>
inline auto ub_constrain (const T &x , const U &ub );

// ./stan/math/prim/constraint/ub_constrain.hpp:143
template <typename T,
    typename U,
    require_all_eigen_t<T,
    U>*>
inline auto ub_constrain (const T &x , const U &ub , std ::decay_t <return_type_t <T , U >>&lp );

// ./stan/math/prim/constraint/ub_constrain.hpp:162
template <typename T,
    typename U,
    require_not_std_vector_t<U>*>
inline auto ub_constrain (const std ::vector <T >&x , const U &ub );

// ./stan/math/prim/constraint/ub_constrain.hpp:182
template <typename T,
    typename U,
    require_not_std_vector_t<U>*>
inline auto ub_constrain (const std ::vector <T >&x , const U &ub , return_type_t <T , U >&lp );

// ./stan/math/prim/constraint/ub_constrain.hpp:202
template <typename T,
    typename U>
inline auto ub_constrain (const std ::vector <T >&x , const std ::vector <U >&ub );

// ./stan/math/prim/constraint/ub_constrain.hpp:223
template <typename T,
    typename U>
inline auto ub_constrain (const std ::vector <T >&x , const std ::vector <U >&ub , return_type_t <T , U >&lp );

// ./stan/math/prim/constraint/ub_constrain.hpp:252
template <bool Jacobian,
    typename T,
    typename U>
inline auto ub_constrain (const T &x , const U &ub , return_type_t <T , U >&lp );

// ./stan/math/prim/constraint/lub_constrain.hpp:44
template <typename T,
    typename L,
    typename U,
    require_all_stan_scalar_t<T,
    L,
    U>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:95
template <typename T,
    typename L,
    typename U,
    require_all_stan_scalar_t<T,
    L,
    U>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:118
template <typename T,
    typename L,
    typename U,
    require_eigen_t<T>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:129
template <typename T,
    typename L,
    typename U,
    require_eigen_t<T>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:142
template <typename T,
    typename L,
    typename U,
    require_all_eigen_t<T,
    L>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:156
template <typename T,
    typename L,
    typename U,
    require_all_eigen_t<T,
    L>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:172
template <typename T,
    typename L,
    typename U,
    require_all_eigen_t<T,
    U>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:186
template <typename T,
    typename L,
    typename U,
    require_all_eigen_t<T,
    U>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:201
template <typename T,
    typename L,
    typename U,
    require_all_eigen_t<T,
    L,
    U>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:223
template <typename T,
    typename L,
    typename U,
    require_all_eigen_t<T,
    L,
    U>*
;

// ./stan/math/prim/constraint/lub_constrain.hpp:246
template <typename T,
    typename L,
    typename U,
    require_all_not_std_vector_t<L,
    U>*>
inline auto lub_constrain (const std ::vector <T >&x , const L &lb , const U &ub );

// ./stan/math/prim/constraint/lub_constrain.hpp:260
template <typename T,
    typename L,
    typename U,
    require_all_not_std_vector_t<L,
    U>*>
inline auto lub_constrain (const std ::vector <T >&x , const L &lb , const U &ub , return_type_t <T , L , U >&lp );

// ./stan/math/prim/constraint/lub_constrain.hpp:275
template <typename T,
    typename L,
    typename U,
    require_not_std_vector_t<L>*>
inline auto lub_constrain (const std ::vector <T >&x , const L &lb , const std ::vector <U >&ub );

// ./stan/math/prim/constraint/lub_constrain.hpp:291
template <typename T,
    typename L,
    typename U,
    require_not_std_vector_t<L>*>
inline auto lub_constrain (const std ::vector <T >&x , const L &lb , const std ::vector <U >&ub , return_type_t <T , L , U >&lp );

// ./stan/math/prim/constraint/lub_constrain.hpp:308
template <typename T,
    typename L,
    typename U,
    require_not_std_vector_t<U>*>
inline auto lub_constrain (const std ::vector <T >&x , const std ::vector <L >&lb , const U &ub );

// ./stan/math/prim/constraint/lub_constrain.hpp:324
template <typename T,
    typename L,
    typename U,
    require_not_std_vector_t<U>*>
inline auto lub_constrain (const std ::vector <T >&x , const std ::vector <L >&lb , const U &ub , return_type_t <T , L , U >&lp );

// ./stan/math/prim/constraint/lub_constrain.hpp:340
template <typename T,
    typename L,
    typename U>
inline auto lub_constrain (const std ::vector <T >&x , const std ::vector <L >&lb , const std ::vector <U >&ub );

// ./stan/math/prim/constraint/lub_constrain.hpp:356
template <typename T,
    typename L,
    typename U>
inline auto lub_constrain (const std ::vector <T >&x , const std ::vector <L >&lb , const std ::vector <U >&ub , return_type_t <T , L , U >&lp );

// ./stan/math/prim/constraint/lub_constrain.hpp:394
template <bool Jacobian,
    typename T,
    typename L,
    typename U>
inline auto lub_constrain (const T &x , const L &lb , const U &ub , return_type_t <T , L , U >&lp );

// ./stan/math/prim/constraint/lub_constrain.hpp:407
template <typename T,
    typename L,
    typename U>
inline auto lub_constrain (const T &x , const std ::tuple <L , U >&bounds );

// ./stan/math/prim/constraint/lub_constrain.hpp:415
template <typename T,
    typename L,
    typename U>
inline auto lub_constrain (const T &x , const std ::tuple <L , U >&bounds , return_type_t <T , L , U >&lp );

// ./stan/math/prim/constraint/lub_constrain.hpp:424
template <bool Jacobian,
    typename T,
    typename L,
    typename U>
inline auto lub_constrain (const T &x , const std ::tuple <L , U >&bounds , return_type_t <T , L , U >&lp );

// ./stan/math/rev/constraint/ub_constrain.hpp:29
template <typename T,
    typename U,
    require_all_stan_scalar_t<T,
    U>*>
inline auto ub_constrain (const T &x , const U &ub );

// ./stan/math/rev/constraint/ub_constrain.hpp:80
template <typename T,
    typename U,
    require_all_stan_scalar_t<T,
    U>*>
inline auto ub_constrain (const T &x , const U &ub , return_type_t <T , U >&lp );

// ./stan/math/rev/constraint/ub_constrain.hpp:135
template <typename T,
    typename U,
    require_matrix_t<T>*>
inline auto ub_constrain (const T &x , const U &ub );

// ./stan/math/rev/constraint/ub_constrain.hpp:184
template <typename T,
    typename U,
    require_matrix_t<T>*>
inline auto ub_constrain (const T &x , const U &ub , return_type_t <T , U >&lp );

// ./stan/math/rev/constraint/ub_constrain.hpp:237
template <typename T,
    typename U,
    require_all_matrix_t<T,
    U>*>
inline auto ub_constrain (const T &x , const U &ub );

// ./stan/math/rev/constraint/ub_constrain.hpp:300
template <typename T,
    typename U,
    require_all_matrix_t<T,
    U>*>
inline auto ub_constrain (const T &x , const U &ub , return_type_t <T , U >&lp );

// ./stan/math/rev/constraint/lub_constrain.hpp:33
template <typename T,
    typename L,
    typename U,
    require_all_stan_scalar_t<T,
    L,
    U>*
;

// ./stan/math/rev/constraint/lub_constrain.hpp:102
template <typename T,
    typename L,
    typename U,
    require_all_stan_scalar_t<T,
    L,
    U>*
;

// ./stan/math/rev/constraint/lub_constrain.hpp:152
template <typename T,
    typename L,
    typename U,
    require_matrix_t<T>*
;

// ./stan/math/rev/constraint/lub_constrain.hpp:195
template <typename T,
    typename L,
    typename U,
    require_matrix_t<T>*
;

// ./stan/math/rev/constraint/lub_constrain.hpp:252
template <typename T,
    typename L,
    typename U,
    require_all_matrix_t<T,
    L>*
;

// ./stan/math/rev/constraint/lub_constrain.hpp:301
template <typename T,
    typename L,
    typename U,
    require_all_matrix_t<T,
    L>*
;

// ./stan/math/rev/constraint/lub_constrain.hpp:362
template <typename T,
    typename L,
    typename U,
    require_all_matrix_t<T,
    U>*
;

// ./stan/math/rev/constraint/lub_constrain.hpp:413
template <typename T,
    typename L,
    typename U,
    require_all_matrix_t<T,
    U>*
;

// ./stan/math/rev/constraint/lub_constrain.hpp:472
template <typename T,
    typename L,
    typename U,
    require_all_matrix_t<T,
    L,
    U>*
;

// ./stan/math/rev/constraint/lub_constrain.hpp:550
template <typename T,
    typename L,
    typename U,
    require_all_matrix_t<T,
    L,
    U>*
;

// ./stan/math/rev/constraint/ordered_constrain.hpp:24
template <typename T,
    require_rev_col_vector_t<T>*>
inline auto ordered_constrain (const T &x );

// ./stan/math/rev/constraint/ordered_constrain.hpp:72
template <typename VarVec,
    require_var_col_vector_t<VarVec>*>
auto ordered_constrain (const VarVec &x , scalar_type_t <VarVec >&lp );

// ./stan/math/rev/constraint/positive_ordered_constrain.hpp:24
template <typename T,
    require_rev_col_vector_t<T>*>
inline auto positive_ordered_constrain (const T &x );

// ./stan/math/rev/constraint/simplex_constrain.hpp:30
template <typename T,
    require_rev_col_vector_t<T>*>
inline auto simplex_constrain (const T &y );

// ./stan/math/rev/constraint/simplex_constrain.hpp:84
template <typename T,
    require_rev_col_vector_t<T>*>
auto simplex_constrain (const T &y , scalar_type_t <T >&lp );

// ./stan/math/rev/constraint/stochastic_column_constrain.hpp:26
template <typename T,
    require_rev_matrix_t<T>*>
inline plain_type_t <T >stochastic_column_constrain (const T &y );

// ./stan/math/rev/constraint/stochastic_column_constrain.hpp:82
template <typename T,
    require_rev_matrix_t<T>*>
inline plain_type_t <T >stochastic_column_constrain (const T &y , scalar_type_t <T >&lp );

// ./stan/math/rev/constraint/stochastic_row_constrain.hpp:25
template <typename T,
    require_rev_matrix_t<T>*>
inline plain_type_t <T >stochastic_row_constrain (const T &y );

// ./stan/math/rev/constraint/stochastic_row_constrain.hpp:78
template <typename T,
    require_rev_matrix_t<T>*>
inline plain_type_t <T >stochastic_row_constrain (const T &y , scalar_type_t <T >&lp );

// ./stan/math/prim/constraint/sum_to_zero_constrain.hpp:38
template <typename Vec,
    require_eigen_col_vector_t<Vec>*>
inline plain_type_t <Vec >sum_to_zero_constrain (const Vec &y );

// ./stan/math/prim/constraint/sum_to_zero_constrain.hpp:88
template <typename Vec,
    require_eigen_col_vector_t<Vec>*>
inline plain_type_t <Vec >sum_to_zero_constrain (const Vec &y , value_type_t <Vec >&lp );

// ./stan/math/prim/constraint/sum_to_zero_constrain.hpp:123
template <bool Jacobian,
    typename Vec,
    require_not_std_vector_t<Vec>*>
inline plain_type_t <Vec >sum_to_zero_constrain (const Vec &y , return_type_t <Vec >&lp );

// ./stan/math/prim/constraint/sum_to_zero_constrain.hpp:157
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto sum_to_zero_constrain (const T &y , return_type_t <T >&lp );

// ./stan/math/rev/constraint/sum_to_zero_constrain.hpp:41
template <typename T,
    require_rev_col_vector_t<T>*>
inline auto sum_to_zero_constrain (T &&y );

// ./stan/math/rev/constraint/sum_to_zero_constrain.hpp:97
template <typename T,
    require_rev_col_vector_t<T>*>
inline auto sum_to_zero_constrain (T &&y , scalar_type_t <T >&lp );

// ./stan/math/rev/constraint/unit_vector_constrain.hpp:29
template <typename T,
    require_rev_col_vector_t<T>*>
inline auto unit_vector_constrain (const T &y );

// ./stan/math/rev/constraint/unit_vector_constrain.hpp:60
template <typename T,
    require_eigen_col_vector_vt<is_var,
    T>*>
inline auto unit_vector_constrain (const T &y , var &lp );

// ./stan/math/rev/constraint/unit_vector_constrain.hpp:78
template <typename T,
    require_var_col_vector_t<T>*>
inline auto unit_vector_constrain (const T &y , var &lp );

// ./stan/math/prim/fun/add_diag.hpp:25
template <typename T_m,
    typename T_a,
    typename 
;

// ./stan/math/prim/fun/all.hpp:35
template <typename ContainerT,
    require_eigen_st<std::is_integral,
    ContainerT>*>
inline bool all (const ContainerT &x );

// ./stan/math/prim/fun/all.hpp:42
template <typename ...Types>
inline bool all (const std ::tuple <Types ...>&x );

// ./stan/math/prim/fun/all.hpp:56
template <typename InnerT>
inline bool all (const std ::vector <InnerT >&x );

// ./stan/math/prim/fun/all.hpp:70
template <typename ...Types>
inline bool all (const std ::tuple <Types ...>&x );

// ./stan/math/prim/fun/any.hpp:35
template <typename ContainerT,
    require_eigen_st<std::is_integral,
    ContainerT>*>
inline bool any (const ContainerT &x );

// ./stan/math/prim/fun/any.hpp:42
template <typename ...Types>
inline bool any (const std ::tuple <Types ...>&x );

// ./stan/math/prim/fun/any.hpp:56
template <typename InnerT>
inline bool any (const std ::vector <InnerT >&x );

// ./stan/math/prim/fun/any.hpp:70
template <typename ...Types>
inline bool any (const std ::tuple <Types ...>&x );

// ./stan/math/prim/fun/resize.hpp:42
template <typename T>
inline void resize (T &x , std ::vector <int >dims );

// ./stan/math/prim/fun/assign.hpp:21
template <int N>
inline void print_mat_size (std ::ostream &o );

// ./stan/math/prim/fun/assign.hpp:44
template <typename T_lhs,
    typename T_rhs,
    require_all_stan_scalar_t<T_lhs,
    T_rhs>*>
inline void assign (T_lhs &x , const T_rhs &y );

// ./stan/math/prim/fun/assign.hpp:65
template <typename T_lhs,
    typename T_rhs,
    require_all_eigen_t<T_lhs,
    T_rhs>*>
inline void assign (T_lhs &&x , const T_rhs &y );

// ./stan/math/prim/fun/append_array.hpp:58
template <typename T1>
inline std ::vector <T1 >append_array (const std ::vector <T1 >&x , const std ::vector <T1 >&y );

// ./stan/math/prim/fun/append_col.hpp:32
template <typename T1,
    typename T2,
    typename >
inline auto append_col (const T1 &A , const T2 &B );

// ./stan/math/prim/fun/append_col.hpp:68
template <typename Scal,
    typename RowVec,
    require_stan_scalar_t<Scal>*
;

// ./stan/math/prim/fun/append_col.hpp:96
template <typename RowVec,
    typename Scal,
    require_t<is_eigen_row_vector<RowVec>>*>
inline Eigen ::Matrix <return_type_t <RowVec ,Scal >,1,Eigen ::Dynamic >append_col (const RowVec &A , const Scal &B );

// ./stan/math/prim/fun/as_bool.hpp:16
template <typename T>
inline bool as_bool (const T &x );

// ./stan/math/prim/fun/as_value_column_vector_or_scalar.hpp:22
template <typename T>
inline auto as_value_column_vector_or_scalar (T &&a );

// ./stan/math/prim/fun/as_value_column_array_or_scalar.hpp:21
template <typename T>
inline auto as_value_column_array_or_scalar (T &&a );

// ./stan/math/prim/fun/as_value_array_or_scalar.hpp:20
template <typename T>
inline auto as_value_array_or_scalar (T &&v );

// ./stan/math/prim/fun/atan2.hpp:22
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
double atan2 (T1 y , T2 x );

// ./stan/math/prim/fun/atan2.hpp:37
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto atan2 (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/mean.hpp:21
template <typename T,
    require_container_t<T>*>
inline return_type_t <T >mean (const T &m );

// ./stan/math/prim/fun/autocorrelation.hpp:61
template <typename T>
void autocorrelation (const std ::vector <T >&y , std ::vector <T >&ac , Eigen ::FFT <T >&fft );

// ./stan/math/prim/fun/autocorrelation.hpp:118
template <typename T,
    typename DerivedA,
    typename DerivedB>
void autocorrelation (const Eigen ::MatrixBase <DerivedA >&y , Eigen ::MatrixBase <DerivedB >&ac , Eigen ::FFT <T >&fft );

// ./stan/math/prim/fun/autocorrelation.hpp:163
template <typename T>
void autocorrelation (const std ::vector <T >&y , std ::vector <T >&ac );

// ./stan/math/prim/fun/autocorrelation.hpp:190
template <typename T,
    typename DerivedA,
    typename DerivedB>
void autocorrelation (const Eigen ::MatrixBase <DerivedA >&y , Eigen ::MatrixBase <DerivedB >&ac );

// ./stan/math/prim/fun/variance.hpp:23
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline value_type_t <EigMat >variance (const EigMat &m );

// ./stan/math/prim/fun/variance.hpp:44
template <typename StdVec,
    require_std_vector_t<StdVec>*>
inline value_type_t <StdVec >variance (const StdVec &v );

// ./stan/math/prim/fun/autocovariance.hpp:32
template <typename T>
void autocovariance (const std ::vector <T >&y , std ::vector <T >&acov , Eigen ::FFT <T >&fft );

// ./stan/math/prim/fun/autocovariance.hpp:65
template <typename T,
    typename DerivedA,
    typename DerivedB>
void autocovariance (const Eigen ::MatrixBase <DerivedA >&y , Eigen ::MatrixBase <DerivedB >&acov , Eigen ::FFT <T >&fft );

// ./stan/math/prim/fun/autocovariance.hpp:88
template <typename T>
void autocovariance (const std ::vector <T >&y , std ::vector <T >&acov );

// ./stan/math/prim/fun/autocovariance.hpp:117
template <typename T,
    typename DerivedA,
    typename DerivedB>
void autocovariance (const Eigen ::MatrixBase <DerivedA >&y , Eigen ::MatrixBase <DerivedB >&acov );

// ./stan/math/prim/fun/binomial_coefficient_log.hpp:79
template <typename T_n,
    typename T_k,
    require_all_stan_scalar_t<T_n,
    T_k>*>
inline return_type_t <T_n ,T_k >binomial_coefficient_log (const T_n n , const T_k k );

// ./stan/math/prim/fun/binomial_coefficient_log.hpp:163
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto binomial_coefficient_log (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/block.hpp:21
template <typename T,
    require_matrix_t<T>*>
inline auto block (const T &m , size_t i , size_t j , size_t nrows , size_t ncols );

// ./stan/math/prim/fun/ceil.hpp:36
template <typename Container,
    require_not_container_st<std::is_arithmetic,
    Container>*>
inline auto ceil (const Container &x );

// ./stan/math/prim/fun/ceil.hpp:52
template <typename Container,
    require_container_st<std::is_arithmetic,
    Container>*>
inline auto ceil (const Container &x );

// ./stan/math/prim/fun/mdivide_left_tri.hpp:26
template <Eigen::UpLoTypeTriView,
    typename T1,
    typename T2,
    require_all_eigen_t<T1,
    T2>*>
inline auto mdivide_left_tri (const T1 &A , const T2 &b );

// ./stan/math/prim/fun/mdivide_left_tri.hpp:53
template <Eigen::UpLoTypeTriView,
    typename T,
    require_eigen_t<T>*>
inline plain_type_t <T >mdivide_left_tri (const T &A );

// ./stan/math/prim/fun/chol2inv.hpp:23
template <typename T,
    require_eigen_t<T>*>
plain_type_t <T >chol2inv (const T &L );

// ./stan/math/prim/fun/cholesky_decompose.hpp:29
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,EigMat ::RowsAtCompileTime ,EigMat ::ColsAtCompileTime >cholesky_decompose (const EigMat &m );

// ./stan/math/prim/fun/choose.hpp:29

inline int choose (int n , int k );

// ./stan/math/prim/fun/choose.hpp:51
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto choose (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/col.hpp:23
template <typename T,
    require_matrix_t<T>*>
inline auto col (const T &m , size_t j );

// ./stan/math/prim/fun/cols.hpp:19
template <typename T,
    require_matrix_t<T>*>
inline int64_t cols (const T &m );

// ./stan/math/prim/fun/columns_dot_product.hpp:24
template <typename Mat1,
    typename Mat2,
    require_all_eigen_t<Mat1,
    Mat2>*>
inline Eigen ::Matrix <return_type_t <Mat1 ,Mat2 >,1,Mat1 ::ColsAtCompileTime >columns_dot_product (const Mat1 &v1 , const Mat2 &v2 );

// ./stan/math/prim/fun/columns_dot_self.hpp:19
template <typename T,
    require_eigen_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,1,T ::ColsAtCompileTime >columns_dot_self (const T &x );

// ./stan/math/prim/fun/complex_schur_decompose.hpp:27
template <typename M,
    require_eigen_dense_dynamic_t<M>*>
inline Eigen ::Matrix <complex_return_t <scalar_type_t <M >>,-1,-1>complex_schur_decompose_u (const M &m );

// ./stan/math/prim/fun/complex_schur_decompose.hpp:51
template <typename M,
    require_eigen_dense_dynamic_t<M>*>
inline Eigen ::Matrix <complex_return_t <scalar_type_t <M >>,-1,-1>complex_schur_decompose_t (const M &m );

// ./stan/math/prim/fun/squared_distance.hpp:23
template <typename Scal1,
    typename Scal2,
    require_all_stan_scalar_t<Scal1,
    Scal2>*>
inline return_type_t <Scal1 ,Scal2 >squared_distance (const Scal1 &x1 , const Scal2 &x2 );

// ./stan/math/prim/fun/squared_distance.hpp:47
template <typename EigVec1,
    typename EigVec2,
    require_all_eigen_vector_t<EigVec1,
    EigVec2>*>
inline return_type_t <EigVec1 ,EigVec2 >squared_distance (const EigVec1 &v1 , const EigVec2 &v2 );

// ./stan/math/prim/fun/gp_exp_quad_cov.hpp:123
template <typename T_x,
    typename T_sigma,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_sigma ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_exp_quad_cov (const std ::vector <T_x >&x , const T_sigma &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/gp_exp_quad_cov.hpp:163
template <typename T_x,
    typename T_sigma,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_sigma ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_exp_quad_cov (const std ::vector <Eigen ::Matrix <T_x , -1, 1>>&x , const T_sigma &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/gp_exp_quad_cov.hpp:206
template <typename T_x1,
    typename T_x2,
    typename T_sigma,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_sigma ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_exp_quad_cov (const std ::vector <T_x1 >&x1 , const std ::vector <T_x2 >&x2 , const T_sigma &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/gp_exp_quad_cov.hpp:255
template <typename T_x1,
    typename T_x2,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_exp_quad_cov (const std ::vector <Eigen ::Matrix <T_x1 , -1, 1>>&x1 , const std ::vector <Eigen ::Matrix <T_x2 , -1, 1>>&x2 , const T_s &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/cov_exp_quad.hpp:17
template <typename T_x,
    typename T_sigma,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_sigma ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >cov_exp_quad (const std ::vector <T_x >&x , const T_sigma &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/cov_exp_quad.hpp:28
template <typename T_x,
    typename T_sigma,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_sigma ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >cov_exp_quad (const std ::vector <T_x >&x , const T_sigma &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/cov_exp_quad.hpp:39
template <typename T_x1,
    typename T_x2,
    typename T_sigma,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_sigma ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >cov_exp_quad (const std ::vector <T_x1 >&x1 , const std ::vector <T_x2 >&x2 , const T_sigma &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/cov_exp_quad.hpp:50
template <typename T_x1,
    typename T_x2,
    typename T_sigma,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_sigma ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >cov_exp_quad (const std ::vector <T_x1 >&x1 , const std ::vector <T_x2 >&x2 , const T_sigma &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/crossprod.hpp:19
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline auto crossprod (const EigMat &M );

// ./stan/math/prim/fun/csr_extract_u.hpp:24
template <typename T>
const std ::vector <int >csr_extract_u (const Eigen ::SparseMatrix <T , Eigen ::RowMajor >&A );

// ./stan/math/prim/fun/csr_extract_u.hpp:43
template <typename T,
    require_eigen_dense_base_t<T>*>
const std ::vector <int >csr_extract_u (const T &A );

// ./stan/math/prim/fun/csr_extract_v.hpp:25
template <typename T>
const std ::vector <int >csr_extract_v (const Eigen ::SparseMatrix <T , Eigen ::RowMajor >&A );

// ./stan/math/prim/fun/csr_extract_v.hpp:47
template <typename T,
    require_eigen_dense_base_t<T>*>
const std ::vector <int >csr_extract_v (const T &A );

// ./stan/math/prim/fun/csr_extract_w.hpp:20
template <typename T>
const Eigen ::Matrix <T ,Eigen ::Dynamic ,1>csr_extract_w (const Eigen ::SparseMatrix <T , Eigen ::RowMajor >&A );

// ./stan/math/prim/fun/csr_extract_w.hpp:42
template <typename T,
    require_eigen_dense_base_t<T>*>
const Eigen ::Matrix <scalar_type_t <T >,Eigen ::Dynamic ,1>csr_extract_w (const T &A );

// ./stan/math/prim/fun/csr_u_to_z.hpp:25

inline int csr_u_to_z (const std ::vector <int >&u , int i );

// ./stan/math/prim/fun/csr_matrix_times_vector.hpp:72
template <typename T1,
    typename T2,
    require_all_not_rev_matrix_t<T1,
    T2>*>
inline Eigen ::Matrix <return_type_t <T1 ,T2 >,Eigen ::Dynamic ,1>csr_matrix_times_vector (int m , int n , const T1 &w , const std ::vector <int >&v , const std ::vector <int >&u , const T2 &b );

// ./stan/math/prim/fun/csr_to_dense_matrix.hpp:33
template <typename T>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >csr_to_dense_matrix (int m , int n , const T &w , const std ::vector <int >&v , const std ::vector <int >&u );

// ./stan/math/prim/fun/cumulative_sum.hpp:24
template <typename T>
inline std ::vector <T >cumulative_sum (const std ::vector <T >&x );

// ./stan/math/prim/fun/cumulative_sum.hpp:48
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
inline auto cumulative_sum (const EigVec &m );

// ./stan/math/prim/fun/determinant.hpp:20
template <typename T,
    require_eigen_vt<std::is_arithmetic,
    T>*>
inline value_type_t <T >determinant (const T &m );

// ./stan/math/prim/fun/diag_matrix.hpp:19
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
inline Eigen ::Matrix <value_type_t <EigVec >,Eigen ::Dynamic ,Eigen ::Dynamic >diag_matrix (const EigVec &v );

// ./stan/math/prim/fun/diag_post_multiply.hpp:22
template <typename T1,
    typename T2,
    require_eigen_t<T1>*>
auto diag_post_multiply (const T1 &m1 , const T2 &m2 );

// ./stan/math/prim/fun/diag_pre_multiply.hpp:22
template <typename T1,
    typename T2,
    require_eigen_vector_t<T1>*>
auto diag_pre_multiply (const T1 &m1 , const T2 &m2 );

// ./stan/math/prim/fun/diagonal.hpp:18
template <typename T,
    require_matrix_t<T>*>
inline auto diagonal (const T &m );

// ./stan/math/prim/fun/distance.hpp:25
template <typename T1,
    typename T2,
    require_all_stan_scalar_t<T1,
    T2>*>
inline return_type_t <T1 ,T2 >distance (const T1 &x1 , const T2 &x2 );

// ./stan/math/prim/fun/distance.hpp:46
template <typename T1,
    typename T2,
    require_all_vector_t<T1,
    T2>*>
inline return_type_t <T1 ,T2 >distance (const T1 &x1 , const T2 &x2 );

// ./stan/math/prim/fun/dot.hpp:11

inline double dot (const std ::vector <double >&x , const std ::vector <double >&y );

// ./stan/math/prim/fun/eigen_comparisons.hpp:27
template <typename T_a,
    typename T_b,
    require_any_eigen_t<T_a,
    T_b>*>
auto OPERATOR (const T_a &a , const T_b &b );

// ./stan/math/prim/fun/eigen_comparisons.hpp:28
template <typename T_a,
    typename T_b,
    require_any_eigen_t<T_a,
    T_b>*>
auto OPERATOR (const T_a &a , const T_b &b );

// ./stan/math/prim/fun/eigen_comparisons.hpp:29
template <typename T_a,
    typename T_b,
    require_any_eigen_t<T_a,
    T_b>*>
auto OPERATOR (const T_a &a , const T_b &b );

// ./stan/math/prim/fun/eigen_comparisons.hpp:30
template <typename T_a,
    typename T_b,
    require_any_eigen_t<T_a,
    T_b>*>
auto OPERATOR (const T_a &a , const T_b &b );

// ./stan/math/prim/fun/eigen_comparisons.hpp:31
template <typename T_a,
    typename T_b,
    require_any_eigen_t<T_a,
    T_b>*>
auto OPERATOR (const T_a &a , const T_b &b );

// ./stan/math/prim/fun/eigen_comparisons.hpp:32
template <typename T_a,
    typename T_b,
    require_any_eigen_t<T_a,
    T_b>*>
auto OPERATOR (const T_a &a , const T_b &b );

// ./stan/math/prim/fun/eigenvalues.hpp:18
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
inline Eigen ::Matrix <complex_return_t <value_type_t <EigMat >>,-1,1>eigenvalues (const EigMat &m );

// ./stan/math/prim/fun/eigenvalues.hpp:41
template <typename EigCplxMat,
    require_eigen_matrix_dynamic_vt<is_complex,
    EigCplxMat>*>
inline Eigen ::Matrix <complex_return_t <value_type_t <EigCplxMat >>,-1,1>eigenvalues (const EigCplxMat &m );

// ./stan/math/prim/fun/eigenvalues_sym.hpp:21
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,-1,1>eigenvalues_sym (const EigMat &m );

// ./stan/math/prim/fun/eigenvectors.hpp:19
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
inline Eigen ::Matrix <complex_return_t <value_type_t <EigMat >>,-1,-1>eigenvectors (const EigMat &m );

// ./stan/math/prim/fun/eigenvectors.hpp:43
template <typename EigCplxMat,
    require_eigen_matrix_dynamic_vt<is_complex,
    EigCplxMat>*>
inline Eigen ::Matrix <complex_return_t <value_type_t <EigCplxMat >>,-1,-1>eigenvectors (const EigCplxMat &m );

// ./stan/math/prim/fun/eigenvectors_sym.hpp:11
template <typename EigMat,
    require_eigen_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >eigenvectors_sym (const EigMat &m );

// ./stan/math/prim/fun/elt_divide.hpp:22
template <typename Mat1,
    typename Mat2,
    require_all_eigen_t<Mat1,
    Mat2>*>
auto elt_divide (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/prim/fun/elt_divide.hpp:41
template <typename Mat,
    typename Scal,
    require_matrix_t<Mat>*>
auto elt_divide (const Mat &m , Scal s );

// ./stan/math/prim/fun/elt_divide.hpp:58
template <typename Scal,
    typename Mat,
    require_stan_scalar_t<Scal>*>
auto elt_divide (Scal s , const Mat &m );

// ./stan/math/prim/fun/elt_divide.hpp:64
template <typename Scal1,
    typename Scal2,
    require_all_stan_scalar_t<Scal1,
    Scal2>*>
auto elt_divide (Scal1 s1 , Scal2 s2 );

// ./stan/math/prim/fun/factor_U.hpp:22
template <typename T_U,
    typename T_CPCs,
    require_eigen_t<T_U>*>
void factor_U (const T_U &U , T_CPCs &&CPCs );

// ./stan/math/prim/fun/factor_cov_matrix.hpp:24
template <typename T_Sigma,
    typename T_CPCs,
    typename T_sds,
    require_eigen_t<T_Sigma>*>
bool factor_cov_matrix (const T_Sigma &Sigma , T_CPCs &&CPCs , T_sds &&sds );

// ./stan/math/prim/fun/fft.hpp:33
template <typename V,
    require_eigen_vector_vt<is_complex,
    V>*
;

// ./stan/math/prim/fun/fft.hpp:64
template <typename V,
    require_eigen_vector_vt<is_complex,
    V>*
;

// ./stan/math/prim/fun/fft.hpp:85
template <typename M,
    require_eigen_dense_dynamic_vt<is_complex,
    M>*
;

// ./stan/math/prim/fun/fft.hpp:107
template <typename M,
    require_eigen_dense_dynamic_vt<is_complex,
    M>*
;

// ./stan/math/prim/fun/fill.hpp:22
template <typename EigMat,
    typename S,
    require_eigen_t<EigMat>*>
inline void fill (EigMat &x , const S &y );

// ./stan/math/prim/fun/fill.hpp:56
template <typename Vec,
    typename S,
    require_std_vector_t<Vec>*>
inline void fill (Vec &x , S &&y );

// ./stan/math/prim/fun/fmod.hpp:21
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
inline double fmod (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/fmod.hpp:36
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto fmod (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/generalized_inverse.hpp:29
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,EigMat ::ColsAtCompileTime ,EigMat ::RowsAtCompileTime >generalized_inverse (const EigMat &G );

// ./stan/math/prim/fun/get_base1.hpp:27
template <typename T>
inline const T &get_base1 (const std ::vector <T >&x , size_t i , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:50
template <typename T>
inline const T &get_base1 (const std ::vector <std ::vector <T >>&x , size_t i1 , size_t i2 , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:74
template <typename T>
inline const T &get_base1 (const std ::vector <std ::vector <std ::vector <T >>>&x , size_t i1 , size_t i2 , size_t i3 , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:100
template <typename T>
inline const T &get_base1 (const std ::vector <std ::vector <std ::vector <std ::vector <T >>>>&x , size_t i1 , size_t i2 , size_t i3 , size_t i4 , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:128
template <typename T>
inline const T &get_base1 (const std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <T >>>>>&x , size_t i1 , size_t i2 , size_t i3 , size_t i4 , size_t i5 , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:158
template <typename T>
inline const T &get_base1 (const std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <T >>>>>>&x , size_t i1 , size_t i2 , size_t i3 , size_t i4 , size_t i5 , size_t i6 , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:189
template <typename T>
inline const T &get_base1 (const std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <T >>>>>>>&x , size_t i1 , size_t i2 , size_t i3 , size_t i4 , size_t i5 , size_t i6 , size_t i7 , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:221
template <typename T>
inline const T &get_base1 (const std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <std ::vector <T >>>>>>>>&x , size_t i1 , size_t i2 , size_t i3 , size_t i4 , size_t i5 , size_t i6 , size_t i7 , size_t i8 , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:251
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,1,Eigen ::Dynamic >get_base1 (const EigMat &x , size_t m , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:276
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline const value_type_t <EigMat >&get_base1 (const EigMat &x , size_t m , size_t n , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1.hpp:300
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
inline const value_type_t <EigVec >&get_base1 (const EigVec &x , size_t m , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1_lhs.hpp:249
template <typename T>
inline Eigen ::Block <Eigen ::Matrix <T ,Eigen ::Dynamic ,Eigen ::Dynamic >>get_base1_lhs (Eigen ::Matrix <T , Eigen ::Dynamic , Eigen ::Dynamic >&x , size_t m , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1_lhs.hpp:274
template <typename T>
inline T &get_base1_lhs (Eigen ::Matrix <T , Eigen ::Dynamic , Eigen ::Dynamic >&x , size_t m , size_t n , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1_lhs.hpp:297
template <typename T>
inline T &get_base1_lhs (Eigen ::Matrix <T , Eigen ::Dynamic , 1>&x , size_t m , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_base1_lhs.hpp:319
template <typename T>
inline T &get_base1_lhs (Eigen ::Matrix <T , 1, Eigen ::Dynamic >&x , size_t n , const char *error_msg , size_t idx );

// ./stan/math/prim/fun/get_imag.hpp:17
template <typename T>
inline T get_imag (const std ::complex <T >&z );

// ./stan/math/prim/fun/get_imag.hpp:29
template <typename Eig,
    require_eigen_vt<is_complex,
    Eig>*>
inline auto get_imag (const Eig &z );

// ./stan/math/prim/fun/get_imag.hpp:41
template <typename StdVec,
    require_std_vector_st<is_complex,
    StdVec>*>
inline auto get_imag (const StdVec &z );

// ./stan/math/prim/fun/get_lp.hpp:9
template <typename T_lp,
    typename T_lp_accum>
inline return_type_t <T_lp ,T_lp_accum >get_lp (const T_lp &lp , const accumulator <T_lp_accum >&lp_accum );

// ./stan/math/prim/fun/get_real.hpp:17
template <typename T>
inline T get_real (const std ::complex <T >&z );

// ./stan/math/prim/fun/get_real.hpp:29
template <typename Eig,
    require_eigen_vt<is_complex,
    Eig>*>
inline auto get_real (const Eig &z );

// ./stan/math/prim/fun/get_real.hpp:41
template <typename StdVec,
    require_std_vector_st<is_complex,
    StdVec>*>
inline auto get_real (const StdVec &z );

// ./stan/math/prim/fun/gp_dot_prod_cov.hpp:36
template <typename T_x,
    typename T_sigma>
Eigen ::Matrix <return_type_t <T_x ,T_sigma >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_dot_prod_cov (const std ::vector <Eigen ::Matrix <T_x , Eigen ::Dynamic , 1>>&x , const T_sigma &sigma );

// ./stan/math/prim/fun/gp_dot_prod_cov.hpp:98
template <typename T_x,
    typename T_sigma>
Eigen ::Matrix <return_type_t <T_x ,T_sigma >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_dot_prod_cov (const std ::vector <T_x >&x , const T_sigma &sigma );

// ./stan/math/prim/fun/gp_dot_prod_cov.hpp:154
template <typename T_x1,
    typename T_x2,
    typename T_sigma>
Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_sigma >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_dot_prod_cov (const std ::vector <Eigen ::Matrix <T_x1 , Eigen ::Dynamic , 1>>&x1 , const std ::vector <Eigen ::Matrix <T_x2 , Eigen ::Dynamic , 1>>&x2 , const T_sigma &sigma );

// ./stan/math/prim/fun/gp_dot_prod_cov.hpp:218
template <typename T_x1,
    typename T_x2,
    typename T_sigma>
Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_sigma >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_dot_prod_cov (const std ::vector <T_x1 >&x1 , const std ::vector <T_x2 >&x2 , const T_sigma &sigma );

// ./stan/math/prim/fun/gp_exponential_cov.hpp:35
template <typename T_x,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_exponential_cov (const std ::vector <T_x >&x , const T_s &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/gp_exponential_cov.hpp:100
template <typename T_x,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_exponential_cov (const std ::vector <Eigen ::Matrix <T_x , -1, 1>>&x , const T_s &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/gp_exponential_cov.hpp:168
template <typename T_x1,
    typename T_x2,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_exponential_cov (const std ::vector <T_x1 >&x1 , const std ::vector <T_x2 >&x2 , const T_s &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/gp_exponential_cov.hpp:245
template <typename T_x1,
    typename T_x2,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_exponential_cov (const std ::vector <Eigen ::Matrix <T_x1 , -1, 1>>&x1 , const std ::vector <Eigen ::Matrix <T_x2 , -1, 1>>&x2 , const T_s &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/gp_matern32_cov.hpp:36
template <typename T_x,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_matern32_cov (const std ::vector <T_x >&x , const T_s &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/gp_matern32_cov.hpp:106
template <typename T_x,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_matern32_cov (const std ::vector <Eigen ::Matrix <T_x , -1, 1>>&x , const T_s &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/gp_matern32_cov.hpp:183
template <typename T_x1,
    typename T_x2,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_matern32_cov (const std ::vector <T_x1 >&x1 , const std ::vector <T_x2 >&x2 , const T_s &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/gp_matern32_cov.hpp:265
template <typename T_x1,
    typename T_x2,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_matern32_cov (const std ::vector <Eigen ::Matrix <T_x1 , -1, 1>>&x1 , const std ::vector <Eigen ::Matrix <T_x2 , -1, 1>>&x2 , const T_s &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/gp_matern52_cov.hpp:38
template <typename T_x,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_matern52_cov (const std ::vector <T_x >&x , const T_s &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/gp_matern52_cov.hpp:109
template <typename T_x,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_matern52_cov (const std ::vector <Eigen ::Matrix <T_x , Eigen ::Dynamic , 1>>&x , const T_s &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/gp_matern52_cov.hpp:185
template <typename T_x1,
    typename T_x2,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_matern52_cov (const std ::vector <T_x1 >&x1 , const std ::vector <T_x2 >&x2 , const T_s &sigma , const T_l &length_scale );

// ./stan/math/prim/fun/gp_matern52_cov.hpp:267
template <typename T_x1,
    typename T_x2,
    typename T_s,
    typename T_l>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_s ,T_l >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_matern52_cov (const std ::vector <Eigen ::Matrix <T_x1 , Eigen ::Dynamic , 1>>&x1 , const std ::vector <Eigen ::Matrix <T_x2 , Eigen ::Dynamic , 1>>&x2 , const T_s &sigma , const std ::vector <T_l >&length_scale );

// ./stan/math/prim/fun/gp_periodic_cov.hpp:44
template <typename T_x,
    typename T_sigma,
    typename T_l,
    typename T_p>
inline typename Eigen ::Matrix <return_type_t <T_x ,T_sigma ,T_l ,T_p >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_periodic_cov (const std ::vector <T_x >&x , const T_sigma &sigma , const T_l &l , const T_p &p );

// ./stan/math/prim/fun/gp_periodic_cov.hpp:123
template <typename T_x1,
    typename T_x2,
    typename T_sigma,
    typename T_l,
    typename T_p>
inline typename Eigen ::Matrix <return_type_t <T_x1 ,T_x2 ,T_sigma ,T_l ,T_p >,Eigen ::Dynamic ,Eigen ::Dynamic >gp_periodic_cov (const std ::vector <T_x1 >&x1 , const std ::vector <T_x2 >&x2 , const T_sigma &sigma , const T_l &l , const T_p &p );

// ./stan/math/prim/fun/grad_F32.hpp:38
template <typename T>
void grad_F32 (T *g , const T &a1 , const T &a2 , const T &a3 , const T &b1 , const T &b2 , const T &z , const T &precision =1e-6, int max_steps =1e5);

// ./stan/math/prim/fun/head.hpp:21
template <typename T,
    require_vector_t<T>*>
inline auto head (const T &v , size_t n );

// ./stan/math/prim/fun/head.hpp:39
template <typename T>
std ::vector <T >head (const std ::vector <T >&sv , size_t n );

// ./stan/math/prim/fun/hypergeometric_2F2.hpp:22
template <typename Ta,
    typename Tb,
    typename Tz,
    require_all_eigen_t<Ta,
    Tb>*>
return_type_t <Ta ,Tb ,Tz >hypergeometric_2F2 (const Ta &a , const Tb &b , const Tz &z );

// ./stan/math/prim/fun/identity_matrix.hpp:18
template <typename T>
inline auto identity_matrix (int K );

// ./stan/math/prim/fun/if_else.hpp:25
template <typename T_true,
    typename T_false>
inline return_type_t <T_true ,T_false >if_else (const bool c , const T_true y_true , const T_false y_false );

// ./stan/math/prim/fun/imag.hpp:17
template <typename T,
    require_autodiff_t<T>>
T imag (const std ::complex <T >&z );

// ./stan/math/prim/fun/initialize.hpp:15
template <typename T,
    typename V,
    require_all_stan_scalar_t<T,
    V>*>
inline void initialize (T &x , V v );

// ./stan/math/prim/fun/initialize.hpp:21
template <typename T,
    int R,
    int C,
    typename V>
inline void initialize (Eigen ::Matrix <T , R , C >&x , const V &v );

// ./stan/math/prim/fun/initialize_fill.hpp:22
template <typename EigMat,
    typename S,
    require_eigen_t<EigMat>*>
inline void initialize_fill (EigMat &x , const S &y );

// ./stan/math/prim/fun/initialize_fill.hpp:56
template <typename Vec,
    typename S,
    require_std_vector_t<Vec>*>
inline void initialize_fill (Vec &x , S &&y );

// ./stan/math/prim/fun/int_step.hpp:26
template <typename T>
inline int int_step (const T &y );

// ./stan/math/prim/fun/inverse_softmax.hpp:36
template <typename Vector,
    require_vector_t<Vector>*>
void inverse_softmax (const Vector &simplex , Vector &y );

// ./stan/math/prim/fun/inverse_spd.hpp:21
template <typename EigMat>
inline Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >inverse_spd (const EigMat &m );

// ./stan/math/prim/fun/isnormal.hpp:20
template <typename ADType,
    require_autodiff_t<ADType>*>
inline bool isnormal (ADType &&v );

// ./stan/math/prim/fun/is_uninitialized.hpp:18
template <typename T>
inline bool is_uninitialized (T x );

// ./stan/math/prim/fun/linspaced_array.hpp:26

inline std ::vector <double >linspaced_array (int K , double low , double high );

// ./stan/math/prim/fun/linspaced_int_array.hpp:33

inline std ::vector <int >linspaced_int_array (int K , int low , int high );

// ./stan/math/prim/fun/linspaced_row_vector.hpp:25

inline auto linspaced_row_vector (int K , double low , double high );

// ./stan/math/prim/fun/linspaced_vector.hpp:25

inline auto linspaced_vector (int K , double low , double high );

// ./stan/math/prim/fun/logb.hpp:23
template <typename T,
    typename >
double logb (const T &x );

// ./stan/math/prim/fun/log_determinant.hpp:20
template <typename EigMat,
    require_eigen_vt<std::is_arithmetic,
    EigMat>*>
inline value_type_t <EigMat >log_determinant (const EigMat &m );

// ./stan/math/prim/fun/log_determinant_ldlt.hpp:16
template <typename T,
    require_not_rev_matrix_t<T>*>
inline value_type_t <T >log_determinant_ldlt (LDLT_factor <T >&A );

// ./stan/math/prim/fun/log_determinant_spd.hpp:23
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline value_type_t <EigMat >log_determinant_spd (const EigMat &m );

// ./stan/math/prim/fun/log_modified_bessel_first_kind.hpp:45
template <typename T1,
    typename T2,
    require_all_stan_scalar_t<T1,
    T2>*>
inline return_type_t <T1 ,T2 ,double >log_modified_bessel_first_kind (const T1 v , const T2 z );

// ./stan/math/prim/fun/log_modified_bessel_first_kind.hpp:234
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto log_modified_bessel_first_kind (const T1 &a , const T2 &b );

// ./stan/math/prim/fun/log_sum_exp_signed.hpp:26
template <typename T1,
    typename T2,
    require_all_stan_scalar_t<T1,
    T2>*>
inline std ::tuple <return_type_t <T1 ,T2 >,int >log_sum_exp_signed (const T1 &a , int a_sign , const T2 &b , int b_sign );

// ./stan/math/prim/fun/logical_and.hpp:29
template <typename T1,
    typename T2>
inline int logical_and (const T1 x1 , const T2 x2 );

// ./stan/math/prim/fun/logical_eq.hpp:18
template <typename T1,
    typename T2>
inline int logical_eq (const T1 x1 , const T2 x2 );

// ./stan/math/prim/fun/logical_gt.hpp:17
template <typename T1,
    typename T2>
inline bool logical_gt (const T1 x1 , const T2 x2 );

// ./stan/math/prim/fun/logical_gte.hpp:17
template <typename T1,
    typename T2>
inline bool logical_gte (const T1 x1 , const T2 x2 );

// ./stan/math/prim/fun/logical_lt.hpp:17
template <typename T1,
    typename T2>
inline bool logical_lt (T1 x1 , T2 x2 );

// ./stan/math/prim/fun/logical_lte.hpp:17
template <typename T1,
    typename T2>
inline bool logical_lte (const T1 x1 , const T2 x2 );

// ./stan/math/prim/fun/logical_negation.hpp:16
template <typename T>
inline int logical_negation (const T &x );

// ./stan/math/prim/fun/logical_neq.hpp:18
template <typename T1,
    typename T2>
inline int logical_neq (const T1 x1 , const T2 x2 );

// ./stan/math/prim/fun/logical_or.hpp:28
template <typename T1,
    typename T2>
inline int logical_or (T1 x1 , T2 x2 );

// ./stan/math/prim/fun/make_nu.hpp:18
template <typename T>
Eigen ::Array <T ,Eigen ::Dynamic ,1>make_nu (const T &eta , size_t K );

// ./stan/math/prim/fun/matrix_exp_pade.hpp:20
template <typename EigMat,
    require_eigen_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,EigMat ::RowsAtCompileTime ,EigMat ::ColsAtCompileTime >matrix_exp_pade (const EigMat &arg );

// ./stan/math/prim/fun/matrix_exp_2x2.hpp:24
template <typename EigMat,
    require_eigen_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >matrix_exp_2x2 (const EigMat &A );

// ./stan/math/prim/fun/matrix_exp.hpp:24
template <typename T,
    typename >
inline plain_type_t <T >matrix_exp (const T &A_in );

// ./stan/math/prim/fun/matrix_exp_multiply.hpp:20
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline Eigen ::Matrix <double ,Eigen ::Dynamic ,EigMat2 ::ColsAtCompileTime >matrix_exp_multiply (const EigMat1 &A , const EigMat2 &B );

// ./stan/math/prim/fun/matrix_power.hpp:22
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,EigMat ::RowsAtCompileTime ,EigMat ::ColsAtCompileTime >matrix_power (const EigMat &M , const int n );

// ./stan/math/prim/fun/matrix_power.hpp:48
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,EigMat ::RowsAtCompileTime ,EigMat ::ColsAtCompileTime >operator ^(const EigMat &M , const int n );

// ./stan/math/prim/fun/size_mvt.hpp:24
template <typename ScalarT,
    require_stan_scalar_t<ScalarT>*>
int64_t size_mvt (const ScalarT &/* unused */);

// ./stan/math/prim/fun/size_mvt.hpp:29
template <typename MatrixT,
    require_matrix_t<MatrixT>*>
int64_t size_mvt (const MatrixT &/* unused */);

// ./stan/math/prim/fun/size_mvt.hpp:34
template <typename MatrixT,
    require_matrix_t<MatrixT>*>
int64_t size_mvt (const std ::vector <MatrixT >&x );

// ./stan/math/prim/fun/max_size_mvt.hpp:24
template <typename T1,
    typename ...Ts>
inline int64_t max_size_mvt (const T1 &x1 , const Ts &...xs );

// ./stan/math/prim/fun/mdivide_left_spd.hpp:25
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline Eigen ::Matrix <return_type_t <EigMat1 ,EigMat2 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_left_spd (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/prim/fun/mdivide_left_tri_low.hpp:29
template <typename T1,
    typename T2,
    require_all_eigen_t<T1,
    T2>*>
inline Eigen ::Matrix <return_type_t <T1 ,T2 >,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left_tri_low (const T1 &A , const T2 &b );

// ./stan/math/prim/fun/mdivide_left_tri_low.hpp:43
template <typename T,
    require_eigen_t<T>*>
inline plain_type_t <T >mdivide_left_tri_low (const T &A );

// ./stan/math/prim/fun/mdivide_right_ldlt.hpp:24
template <typename EigMat,
    typename T,
    require_all_matrix_t<EigMat,
    T>*>
inline auto mdivide_right_ldlt (const EigMat &b , LDLT_factor <T >&A );

// ./stan/math/prim/fun/mdivide_right_ldlt.hpp:47
template <typename EigMat,
    typename T,
    require_all_matrix_t<EigMat,
    T>*>
inline Eigen ::Matrix <double ,EigMat ::RowsAtCompileTime ,T ::ColsAtCompileTime >mdivide_right_ldlt (const EigMat &b , LDLT_factor <T >&A );

// ./stan/math/prim/fun/mdivide_right_spd.hpp:26
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline Eigen ::Matrix <return_type_t <EigMat1 ,EigMat2 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_right_spd (const EigMat1 &b , const EigMat2 &A );

// ./stan/math/prim/fun/mdivide_right_tri.hpp:27
template <Eigen::UpLoTypeTriView,
    typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline auto mdivide_right_tri (const EigMat1 &b , const EigMat2 &A );

// ./stan/math/prim/fun/mdivide_right_tri_low.hpp:24
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline Eigen ::Matrix <return_type_t <EigMat1 ,EigMat2 >,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_right_tri_low (const EigMat1 &b , const EigMat2 &A );

// ./stan/math/prim/fun/min.hpp:23
template <typename T1,
    typename T2,
    require_all_arithmetic_t<T1,
    T2>*>
auto min (T1 x , T2 y );

// ./stan/math/prim/fun/min.hpp:40
template <typename T,
    require_container_t<T>*>
inline value_type_t <T >min (const T &m );

// ./stan/math/prim/fun/minus.hpp:17
template <typename T,
    require_not_std_vector_t<T>*>
inline auto minus (const T &x );

// ./stan/math/prim/fun/minus.hpp:29
template <typename T>
inline auto minus (const std ::vector <T >&x );

// ./stan/math/prim/fun/modulus.hpp:12

inline int modulus (int x , int y );

// ./stan/math/prim/fun/one_hot_array.hpp:20

inline std ::vector <double >one_hot_array (int K , int k );

// ./stan/math/prim/fun/one_hot_int_array.hpp:20

inline std ::vector <int >one_hot_int_array (int K , int k );

// ./stan/math/prim/fun/one_hot_row_vector.hpp:20

inline Eigen ::RowVectorXd one_hot_row_vector (int K , int k );

// ./stan/math/prim/fun/one_hot_vector.hpp:20

inline Eigen ::VectorXd one_hot_vector (int K , int k );

// ./stan/math/prim/fun/ones_array.hpp:17

inline std ::vector <double >ones_array (int K );

// ./stan/math/prim/fun/ones_int_array.hpp:17

inline std ::vector <int >ones_int_array (int K );

// ./stan/math/prim/fun/ones_row_vector.hpp:17

inline auto ones_row_vector (int K );

// ./stan/math/prim/fun/ones_vector.hpp:17

inline auto ones_vector (int K );

// ./stan/math/prim/fun/Phi_approx.hpp:23

inline double Phi_approx (double x );

// ./stan/math/prim/fun/Phi_approx.hpp:35

inline double Phi_approx (int x );

// ./stan/math/prim/fun/Phi_approx.hpp:65
template <typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>*>
inline auto Phi_approx (const T &x );

// ./stan/math/prim/fun/plus.hpp:16
template <typename T>
inline T plus (T &&x );

// ./stan/math/prim/fun/poisson_binomial_log_probs.hpp:24
template <typename T_theta,
    typename T_scalar>
plain_type_t <T_theta >poisson_binomial_log_probs (int y , const T_theta &theta );

// ./stan/math/prim/fun/poisson_binomial_log_probs.hpp:55
template <typename T_y,
    typename T_theta,
    require_vt_integral<T_y>*>
auto poisson_binomial_log_probs (const T_y &y , const T_theta &theta );

// ./stan/math/prim/fun/promote_scalar.hpp:21
template <typename PromotionScalar,
    typename UnPromotedType,
    require_constructible_t<PromotionScalar,
    UnPromotedType>*>
inline constexpr auto promote_scalar (UnPromotedType &&x );

// ./stan/math/prim/fun/promote_scalar.hpp:37
template <typename PromotionScalar,
    typename UnPromotedType,
    require_same_t<PromotionScalar,
    scalar_type_t<UnPromotedType>>*>
inline constexpr auto promote_scalar (UnPromotedType &&x )noexcept ;

// ./stan/math/prim/fun/promote_scalar.hpp:52
template <typename PromotionScalar,
    typename UnPromotedType,
    require_eigen_t<UnPromotedType>*
;

// ./stan/math/prim/fun/promote_scalar.hpp:61
template <typename PromotionScalars,
    typename UnPromotedTypes,
    require_all_tuple_t<PromotionScalars,
    UnPromotedTypes>*>
inline constexpr promote_scalar_t <PromotionScalars ,UnPromotedTypes >promote_scalar (UnPromotedTypes &&x );

// ./stan/math/prim/fun/promote_scalar.hpp:75
template <typename PromotionScalar,
    typename UnPromotedType,
    require_std_vector_t<UnPromotedType>*
;

// ./stan/math/prim/fun/promote_scalar.hpp:97
template <typename PromotionScalars,
    typename UnPromotedTypes,
    require_all_tuple_t<PromotionScalars,
    UnPromotedTypes>*,
    require_not_same_t<PromotionScalars,
    UnPromotedTypes>*>
inline constexpr promote_scalar_t <PromotionScalars ,UnPromotedTypes >promote_scalar (UnPromotedTypes &&x );

// ./stan/math/prim/fun/pseudo_eigenvalues.hpp:10
template <typename T>
Eigen ::Matrix <T ,-1,-1>pseudo_eigenvalues (const Eigen ::Matrix <T , -1, -1>&m );

// ./stan/math/prim/fun/pseudo_eigenvectors.hpp:10
template <typename T>
Eigen ::Matrix <T ,-1,-1>pseudo_eigenvectors (const Eigen ::Matrix <T , -1, -1>&m );

// ./stan/math/prim/fun/qr_Q.hpp:19
template <typename EigMat,
    require_eigen_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >qr_Q (const EigMat &m );

// ./stan/math/prim/fun/qr_R.hpp:18
template <typename EigMat,
    require_eigen_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >qr_R (const EigMat &m );

// ./stan/math/prim/fun/qr_thin_Q.hpp:18
template <typename EigMat,
    require_eigen_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >qr_thin_Q (const EigMat &m );

// ./stan/math/prim/fun/qr_thin_R.hpp:18
template <typename EigMat,
    require_eigen_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >qr_thin_R (const EigMat &m );

// ./stan/math/prim/fun/quad_form.hpp:26
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline auto quad_form (const EigMat1 &A , const EigMat2 &B );

// ./stan/math/prim/fun/quad_form.hpp:51
template <typename EigMat,
    typename ColVec,
    require_eigen_t<EigMat>*>
inline value_type_t <EigMat >quad_form (const EigMat &A , const ColVec &B );

// ./stan/math/prim/fun/quad_form_diag.hpp:11
template <typename EigMat,
    typename EigVec,
    require_eigen_t<EigMat>*>
inline auto quad_form_diag (const EigMat &mat , const EigVec &vec );

// ./stan/math/prim/fun/quad_form_sym.hpp:25
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline plain_type_t <EigMat2 >quad_form_sym (const EigMat1 &A , const EigMat2 &B );

// ./stan/math/prim/fun/quad_form_sym.hpp:52
template <typename EigMat,
    typename ColVec,
    require_eigen_t<EigMat>*>
inline value_type_t <EigMat >quad_form_sym (const EigMat &A , const ColVec &B );

// ./stan/math/prim/fun/quantile.hpp:30
template <typename T,
    require_vector_t<T>*>
inline double quantile (const T &samples_vec , const double p );

// ./stan/math/prim/fun/quantile.hpp:78
template <typename T,
    typename Tp,
    require_all_vector_t<T,
    Tp>*>
inline std ::vector <double >quantile (const T &samples_vec , const Tp &ps );

// ./stan/math/prim/fun/rank.hpp:19
template <typename C,
    require_container_t<C>*>
inline int rank (const C &v , int s );

// ./stan/math/prim/fun/real.hpp:17
template <typename T,
    require_autodiff_t<T>>
T real (const std ::complex <T >&z );

// ./stan/math/prim/fun/rep_array.hpp:11
template <typename T_ret,
    typename In,
    require_std_vector_t<T_ret>*>
inline std ::vector <plain_type_t <In >>rep_array (const In &x , int n );

// ./stan/math/prim/fun/rep_array.hpp:17
template <typename In>
inline std ::vector <plain_type_t <In >>rep_array (const In &x , int n );

// ./stan/math/prim/fun/rep_matrix.hpp:21
template <typename Ret,
    typename T,
    require_eigen_matrix_dynamic_vt<is_stan_scalar,
    Ret>*>
inline auto rep_matrix (const T &x , int m , int n );

// ./stan/math/prim/fun/rep_matrix.hpp:39
template <typename T,
    require_stan_scalar_t<T>*>
inline auto rep_matrix (const T &x , int m , int n );

// ./stan/math/prim/fun/rep_matrix.hpp:53
template <typename Vec,
    require_eigen_vector_t<Vec>*>
inline auto rep_matrix (const Vec &x , int n );

// ./stan/math/prim/fun/rep_row_vector.hpp:10
template <typename T_ret,
    typename T,
    require_eigen_row_vector_t<T_ret>*>
inline auto rep_row_vector (const T &x , int n );

// ./stan/math/prim/fun/rep_row_vector.hpp:17
template <typename T,
    require_stan_scalar_t<T>*>
inline auto rep_row_vector (const T &x , int n );

// ./stan/math/prim/fun/rep_vector.hpp:10
template <typename T_ret,
    typename T,
    require_eigen_col_vector_t<T_ret>*>
inline auto rep_vector (const T &x , int n );

// ./stan/math/prim/fun/rep_vector.hpp:17
template <typename T,
    require_stan_scalar_t<T>*>
inline auto rep_vector (const T &x , int n );

// ./stan/math/prim/fun/reverse.hpp:19
template <typename T>
inline std ::vector <T >reverse (const std ::vector <T >&x );

// ./stan/math/prim/fun/reverse.hpp:34
template <typename T,
    typename >
inline auto reverse (const T &x );

// ./stan/math/prim/fun/row.hpp:23
template <typename T,
    require_matrix_t<T>*>
inline auto row (const T &m , size_t i );

// ./stan/math/prim/fun/rows.hpp:19
template <typename T,
    require_matrix_t<T>*>
inline int64_t rows (const T &m );

// ./stan/math/prim/fun/rows_dot_product.hpp:24
template <typename Mat1,
    typename Mat2,
    require_all_eigen_t<Mat1,
    Mat2>*>
inline Eigen ::Matrix <return_type_t <Mat1 ,Mat2 >,Mat1 ::RowsAtCompileTime ,1>rows_dot_product (const Mat1 &v1 , const Mat2 &v2 );

// ./stan/math/prim/fun/rows_dot_self.hpp:18
template <typename T,
    require_eigen_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,T ::RowsAtCompileTime ,1>rows_dot_self (const T &x );

// ./stan/math/prim/fun/scalbn.hpp:11
template <typename T,
    typename >
double scalbn (const T &x , int n );

// ./stan/math/prim/fun/scale_matrix_exp_multiply.hpp:27
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_vt<std::is_arithmetic,
    EigMat1,
    EigMat2>*>
inline Eigen ::Matrix <double ,Eigen ::Dynamic ,EigMat2 ::ColsAtCompileTime >scale_matrix_exp_multiply (const double &t , const EigMat1 &A , const EigMat2 &B );

// ./stan/math/prim/fun/sd.hpp:24
template <typename T,
    require_container_t<T>*>
inline auto sd (const T &m );

// ./stan/math/prim/fun/segment.hpp:17
template <typename Vec,
    require_vector_t<Vec>*>
inline auto segment (const Vec &v , size_t i , size_t n );

// ./stan/math/prim/fun/segment.hpp:29
template <typename T>
std ::vector <T >segment (const std ::vector <T >&sv , size_t i , size_t n );

// ./stan/math/prim/fun/signbit.hpp:20
template <typename ADType,
    require_autodiff_t<ADType>*>
inline bool signbit (ADType &&v );

// ./stan/math/prim/fun/singular_values.hpp:21
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
auto singular_values (const EigMat &m );

// ./stan/math/prim/fun/sort_asc.hpp:37
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
inline plain_type_t <EigVec >sort_asc (EigVec &&xs );

// ./stan/math/prim/fun/sort_desc.hpp:38
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
inline plain_type_t <EigVec >sort_desc (EigVec &&xs );

// ./stan/math/prim/fun/sort_indices.hpp:63
template <bool ascending,
    typename C>
std ::vector <int >sort_indices (const C &xs );

// ./stan/math/prim/fun/sort_indices_asc.hpp:20
template <typename C>
std ::vector <int >sort_indices_asc (const C &xs );

// ./stan/math/prim/fun/sort_indices_desc.hpp:20
template <typename C>
std ::vector <int >sort_indices_desc (const C &xs );

// ./stan/math/prim/fun/stan_print.hpp:13
template <typename T,
    require_not_container_t<T>*>
void stan_print (std ::ostream *o , const T &x );

// ./stan/math/prim/fun/stan_print.hpp:19
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
void stan_print (std ::ostream *o , const EigVec &x );

// ./stan/math/prim/fun/stan_print.hpp:33
template <typename EigMat,
    require_eigen_t<EigMat>*>
void stan_print (std ::ostream *o , const EigMat &x );

// ./stan/math/prim/fun/stan_print.hpp:56
template <typename T,
    require_tuple_t<T>*>
void stan_print (std ::ostream *o , const T &x );

// ./stan/math/prim/fun/stan_print.hpp:59
template <typename T,
    require_std_vector_t<T>*>
void stan_print (std ::ostream *o , const T &x );

// ./stan/math/prim/fun/stan_print.hpp:71
template <typename T,
    require_tuple_t<T>*>
void stan_print (std ::ostream *o , const T &x );

// ./stan/math/prim/fun/step.hpp:30
template <typename T,
    require_stan_scalar_t<T>*>
inline T step (const T &y );

// ./stan/math/prim/fun/step.hpp:60
template <typename T,
    require_container_t<T>*>
inline auto step (const T &x );

// ./stan/math/prim/fun/sub_col.hpp:20
template <typename T,
    require_matrix_t<T>*>
inline auto sub_col (const T &m , size_t i , size_t j , size_t nrows );

// ./stan/math/prim/fun/sub_row.hpp:20
template <typename T,
    require_matrix_t<T>*>
inline auto sub_row (const T &m , size_t i , size_t j , size_t ncols );

// ./stan/math/prim/fun/svd_U.hpp:17
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >svd_U (const EigMat &m );

// ./stan/math/prim/fun/svd_V.hpp:17
template <typename EigMat,
    require_eigen_matrix_dynamic_t<EigMat>*>
Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >svd_V (const EigMat &m );

// ./stan/math/prim/fun/symmetrize_from_lower_tri.hpp:18
template <typename T,
    require_eigen_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >symmetrize_from_lower_tri (const T &m );

// ./stan/math/prim/fun/symmetrize_from_upper_tri.hpp:18
template <typename T,
    require_eigen_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >symmetrize_from_upper_tri (const T &m );

// ./stan/math/prim/fun/tail.hpp:22
template <typename T,
    require_vector_t<T>*>
inline auto tail (const T &v , size_t n );

// ./stan/math/prim/fun/tail.hpp:40
template <typename T>
std ::vector <T >tail (const std ::vector <T >&sv , size_t n );

// ./stan/math/prim/fun/to_complex.hpp:21
template <typename T>
constexpr inline std ::complex <stan ::real_return_t <T ,S >>to_complex (const T &re =0, const S &im =0);

// ./stan/math/prim/fun/to_complex.hpp:38
template <typename T1,
    typename T2,
    require_any_container_t<T1,
    T2>*>
inline auto to_complex (const T1 &re , const T2 &im );

// ./stan/math/prim/fun/to_int.hpp:18
template <typename T,
    require_integral_t<T>*>
inline T to_int (T x );

// ./stan/math/prim/fun/to_int.hpp:41
template <typename T,
    require_floating_point_t<T>*>
inline int to_int (T x );

// ./stan/math/prim/fun/to_int.hpp:72
template <typename Container,
    require_std_vector_st<std::is_arithmetic,
    Container>*>
inline auto to_int (const Container &x );

// ./stan/math/prim/fun/to_matrix.hpp:21
template <typename EigMat,
    require_eigen_dense_dynamic_t<EigMat>*>
inline EigMat to_matrix (EigMat &&x );

// ./stan/math/prim/fun/to_matrix.hpp:36
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
inline auto to_matrix (EigVec &&matrix );

// ./stan/math/prim/fun/to_matrix.hpp:50
template <typename T>
inline Eigen ::Matrix <T ,Eigen ::Dynamic ,Eigen ::Dynamic >to_matrix (const std ::vector <Eigen ::Matrix <T , 1, Eigen ::Dynamic >>&x );

// ./stan/math/prim/fun/to_matrix.hpp:75
template <typename T>
inline Eigen ::Matrix <return_type_t <T ,double >,Eigen ::Dynamic ,Eigen ::Dynamic >to_matrix (const std ::vector <std ::vector <T >>&x );

// ./stan/math/prim/fun/to_matrix.hpp:106
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >to_matrix (EigMat &&x , int m , int n );

// ./stan/math/prim/fun/to_matrix.hpp:129
template <typename T>
inline auto to_matrix (const std ::vector <T >&x , int m , int n );

// ./stan/math/prim/fun/to_matrix.hpp:148

inline Eigen ::Matrix <double ,Eigen ::Dynamic ,Eigen ::Dynamic >to_matrix (const std ::vector <int >&x , int m , int n );

// ./stan/math/prim/fun/to_matrix.hpp:177
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,Eigen ::Dynamic ,Eigen ::Dynamic >to_matrix (EigMat &&x , int m , int n , bool col_major );

// ./stan/math/prim/fun/to_matrix.hpp:206
template <typename T>
inline Eigen ::Matrix <return_type_t <T ,double >,Eigen ::Dynamic ,Eigen ::Dynamic >to_matrix (const std ::vector <T >&x , int m , int n , bool col_major );

// ./stan/math/prim/fun/to_row_vector.hpp:14
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline Eigen ::Matrix <value_type_t <EigMat >,1,Eigen ::Dynamic >to_row_vector (const EigMat &matrix );

// ./stan/math/prim/fun/to_row_vector.hpp:27
template <typename T>
inline Eigen ::Matrix <T ,1,Eigen ::Dynamic >to_row_vector (const std ::vector <T >&vec );

// ./stan/math/prim/fun/to_row_vector.hpp:34

inline Eigen ::Matrix <double ,1,Eigen ::Dynamic >to_row_vector (const std ::vector <int >&vec );

// ./stan/math/prim/fun/trace_gen_inv_quad_form_ldlt.hpp:33
template <typename EigMat1,
    typename T2,
    typename EigMat3,
    require_not_col_vector_t<EigMat1>*>
inline return_type_t <EigMat1 ,T2 ,EigMat3 >trace_gen_inv_quad_form_ldlt (const EigMat1 &D , LDLT_factor <T2 >&A , const EigMat3 &B );

// ./stan/math/prim/fun/trace_gen_inv_quad_form_ldlt.hpp:66
template <typename EigVec,
    typename T,
    typename EigMat,
    require_col_vector_t<EigVec>*>
inline return_type_t <EigVec ,T ,EigMat >trace_gen_inv_quad_form_ldlt (const EigVec &D , LDLT_factor <T >&A , const EigMat &B );

// ./stan/math/prim/fun/trace_gen_quad_form.hpp:32
template <typename TD,
    typename TA,
    typename TB,
    typename >
inline auto trace_gen_quad_form (const TD &D , const TA &A , const TB &B );

// ./stan/math/prim/fun/trace_gen_quad_form.hpp:63
template <typename EigMatD,
    typename EigMatA,
    typename EigMatB,
    require_all_eigen_vt<std::is_arithmetic,
    EigMatD,
    EigMatA,
    EigMatB>*>
inline double trace_gen_quad_form (const EigMatD &D , const EigMatA &A , const EigMatB &B );

// ./stan/math/prim/fun/trace_inv_quad_form_ldlt.hpp:26
template <typename T,
    typename EigMat2,
    typename >
inline return_type_t <T ,EigMat2 >trace_inv_quad_form_ldlt (LDLT_factor <T >&A , const EigMat2 &B );

// ./stan/math/prim/fun/trace_quad_form.hpp:24
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_vt<std::is_arithmetic,
    EigMat1,
    EigMat2>*>
inline return_type_t <EigMat1 ,EigMat2 >trace_quad_form (const EigMat1 &A , const EigMat2 &B );

// ./stan/math/prim/fun/uniform_simplex.hpp:18

inline auto uniform_simplex (int K );

// ./stan/math/prim/fun/unitspaced_array.hpp:23

inline std ::vector <int >unitspaced_array (int low , int high );

// ./stan/math/prim/fun/vec_concat.hpp:56
template <typename Vec,
    typename ...Args>
inline auto vec_concat (const Vec &v1 , const Args &...args );

// ./stan/math/prim/fun/zeros_array.hpp:17

inline std ::vector <double >zeros_array (int K );

// ./stan/math/prim/fun/zeros_int_array.hpp:17

inline std ::vector <int >zeros_int_array (int K );

// ./stan/math/prim/fun/zeros_row_vector.hpp:17

inline auto zeros_row_vector (int K );

// ./stan/math/prim/fun/zeros_vector.hpp:17

inline auto zeros_vector (int K );

// ./stan/math/prim/constraint/cholesky_corr_constrain.hpp:15
template <typename EigVec,
    require_eigen_col_vector_t<EigVec>*>
inline Eigen ::Matrix <value_type_t <EigVec >,Eigen ::Dynamic ,Eigen ::Dynamic >cholesky_corr_constrain (const EigVec &y , int K );

// ./stan/math/prim/constraint/cholesky_corr_constrain.hpp:46
template <typename EigVec,
    require_eigen_vector_t<EigVec>*>
inline Eigen ::Matrix <value_type_t <EigVec >,Eigen ::Dynamic ,Eigen ::Dynamic >cholesky_corr_constrain (const EigVec &y , int K , return_type_t <EigVec >&lp );

// ./stan/math/prim/constraint/cholesky_corr_constrain.hpp:93
template <bool Jacobian,
    typename T,
    require_not_std_vector_t<T>*>
inline auto cholesky_corr_constrain (const T &y , int K , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cholesky_corr_constrain.hpp:118
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto cholesky_corr_constrain (const T &y , int K , return_type_t <T >&lp );

// ./stan/math/prim/constraint/corr_free.hpp:27
template <typename T>
inline T corr_free (const T &y );

// ./stan/math/prim/constraint/cholesky_corr_free.hpp:13
template <typename T,
    require_eigen_t<T>*>
auto cholesky_corr_free (const T &x );

// ./stan/math/prim/constraint/cholesky_corr_free.hpp:43
template <typename T,
    require_std_vector_t<T>*>
auto cholesky_corr_free (const T &x );

// ./stan/math/prim/constraint/cholesky_factor_constrain.hpp:28
template <typename T,
    require_eigen_col_vector_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >cholesky_factor_constrain (const T &x , int M , int N );

// ./stan/math/prim/constraint/cholesky_factor_constrain.hpp:73
template <typename T,
    require_eigen_vector_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >cholesky_factor_constrain (const T &x , int M , int N , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cholesky_factor_constrain.hpp:108
template <bool Jacobian,
    typename T,
    require_not_std_vector_t<T>*>
inline auto cholesky_factor_constrain (const T &x , int M , int N , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cholesky_factor_constrain.hpp:138
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto cholesky_factor_constrain (const T &x , int M , int N , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cholesky_factor_free.hpp:25
template <typename T,
    require_eigen_t<T>*>
Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,1>cholesky_factor_free (const T &y );

// ./stan/math/prim/constraint/cholesky_factor_free.hpp:58
template <typename T,
    require_std_vector_t<T>*>
auto cholesky_factor_free (const T &x );

// ./stan/math/prim/constraint/corr_matrix_constrain.hpp:39
template <typename T,
    require_eigen_col_vector_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >corr_matrix_constrain (const T &x , Eigen ::Index k );

// ./stan/math/prim/constraint/corr_matrix_constrain.hpp:68
template <typename T,
    require_eigen_col_vector_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >corr_matrix_constrain (const T &x , Eigen ::Index k , return_type_t <T >&lp );

// ./stan/math/prim/constraint/corr_matrix_constrain.hpp:96
template <bool Jacobian,
    typename T,
    require_not_std_vector_t<T>*>
inline auto corr_matrix_constrain (const T &x , Eigen ::Index k , return_type_t <T >&lp );

// ./stan/math/prim/constraint/corr_matrix_constrain.hpp:125
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto corr_matrix_constrain (const T &y , int K , return_type_t <T >&lp );

// ./stan/math/prim/constraint/corr_matrix_free.hpp:33
template <typename T,
    require_eigen_t<T>*>
Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,1>corr_matrix_free (const T &y );

// ./stan/math/prim/constraint/corr_matrix_free.hpp:62
template <typename T,
    require_std_vector_t<T>*>
auto corr_matrix_free (const T &x );

// ./stan/math/prim/constraint/cov_matrix_constrain.hpp:30
template <typename T,
    require_eigen_col_vector_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >cov_matrix_constrain (const T &x , Eigen ::Index K );

// ./stan/math/prim/constraint/cov_matrix_constrain.hpp:63
template <typename T,
    require_eigen_col_vector_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >cov_matrix_constrain (const T &x , Eigen ::Index K , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cov_matrix_constrain.hpp:107
template <bool Jacobian,
    typename T,
    require_not_std_vector_t<T>*>
inline auto cov_matrix_constrain (const T &x , Eigen ::Index K , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cov_matrix_constrain.hpp:135
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto cov_matrix_constrain (const T &x , Eigen ::Index K , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cov_matrix_constrain_lkj.hpp:33
template <typename T,
    require_eigen_vector_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >cov_matrix_constrain_lkj (const T &x , size_t k );

// ./stan/math/prim/constraint/cov_matrix_constrain_lkj.hpp:57
template <typename T,
    require_eigen_vector_t<T>*>
inline Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,Eigen ::Dynamic >cov_matrix_constrain_lkj (const T &x , size_t k , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cov_matrix_constrain_lkj.hpp:86
template <bool Jacobian,
    typename T,
    require_not_std_vector_t<T>*>
inline auto cov_matrix_constrain_lkj (const T &x , size_t k , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cov_matrix_constrain_lkj.hpp:116
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto cov_matrix_constrain_lkj (const T &x , size_t k , return_type_t <T >&lp );

// ./stan/math/prim/constraint/cov_matrix_free.hpp:37
template <typename T,
    require_eigen_t<T>*>
Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,1>cov_matrix_free (const T &y );

// ./stan/math/prim/constraint/cov_matrix_free.hpp:68
template <typename T,
    require_std_vector_t<T>*>
auto cov_matrix_free (const T &x );

// ./stan/math/prim/constraint/cov_matrix_free_lkj.hpp:30
template <typename T,
    require_eigen_t<T>*>
Eigen ::Matrix <value_type_t <T >,Eigen ::Dynamic ,1>cov_matrix_free_lkj (const T &y );

// ./stan/math/prim/constraint/cov_matrix_free_lkj.hpp:59
template <typename T,
    require_std_vector_t<T>*>
auto cov_matrix_free_lkj (const T &x );

// ./stan/math/prim/constraint/lb_free.hpp:28
template <typename T,
    typename L,
    require_not_std_vector_t<T>*>
inline auto lb_free (T &&y , L &&lb );

// ./stan/math/prim/constraint/lb_free.hpp:56
template <typename T,
    typename L,
    require_all_eigen_t<T,
    L>*>
inline auto lb_free (T &&y , L &&lb );

// ./stan/math/prim/constraint/lb_free.hpp:83
template <typename T,
    typename L,
    require_not_std_vector_t<L>*>
inline auto lb_free (const std ::vector <T >y , const L &lb );

// ./stan/math/prim/constraint/lb_free.hpp:106
template <typename T,
    typename L>
inline auto lb_free (const std ::vector <T >y , const std ::vector <L >&lb );

// ./stan/math/prim/constraint/ub_free.hpp:33
template <typename T,
    typename U,
    require_not_std_vector_t<T>*>
inline auto ub_free (T &&y , U &&ub );

// ./stan/math/prim/constraint/ub_free.hpp:62
template <typename T,
    typename U,
    require_all_eigen_t<T,
    U>*>
inline auto ub_free (T &&y , U &&ub );

// ./stan/math/prim/constraint/ub_free.hpp:90
template <typename T,
    typename U,
    require_not_std_vector_t<U>*>
inline auto ub_free (const std ::vector <T >y , const U &ub );

// ./stan/math/prim/constraint/ub_free.hpp:114
template <typename T,
    typename U>
inline auto ub_free (const std ::vector <T >y , const std ::vector <U >&ub );

// ./stan/math/prim/constraint/lub_free.hpp:43
template <typename T,
    typename L,
    typename U,
    require_not_std_vector_t<T>*>
inline auto lub_free (T &&y , L &&lb , U &&ub );

// ./stan/math/prim/constraint/lub_free.hpp:68
template <typename T,
    typename L,
    typename U,
    require_all_eigen_t<T,
    L>*>
inline auto lub_free (T &&y , L &&lb , U &&ub );

// ./stan/math/prim/constraint/lub_free.hpp:88
template <typename T,
    typename L,
    typename U,
    require_all_eigen_t<T,
    U>*>
inline auto lub_free (T &&y , L &&lb , U &&ub );

// ./stan/math/prim/constraint/lub_free.hpp:108
template <typename T,
    typename L,
    typename U,
    require_all_eigen_t<T,
    L,
    U>*>
inline auto lub_free (T &&y , L &&lb , U &&ub );

// ./stan/math/prim/constraint/lub_free.hpp:129
template <typename T,
    typename L,
    typename U,
    require_all_not_std_vector_t<L,
    U>*>
inline auto lub_free (const std ::vector <T >y , const L &lb , const U &ub );

// ./stan/math/prim/constraint/lub_free.hpp:142
template <typename T,
    typename L,
    typename U,
    require_not_std_vector_t<L>*>
inline auto lub_free (const std ::vector <T >y , const L &lb , const std ::vector <U >&ub );

// ./stan/math/prim/constraint/lub_free.hpp:157
template <typename T,
    typename L,
    typename U,
    require_not_std_vector_t<U>*>
inline auto lub_free (const std ::vector <T >y , const std ::vector <L >&lb , const U &ub );

// ./stan/math/prim/constraint/lub_free.hpp:172
template <typename T,
    typename L,
    typename U>
inline auto lub_free (const std ::vector <T >y , const std ::vector <L >&lb , const std ::vector <U >&ub );

// ./stan/math/prim/constraint/lub_free.hpp:187
template <typename T,
    typename L,
    typename U>
inline auto lub_free (T &&y , const std ::tuple <L , U >&bounds );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:41
template <typename T,
    typename M,
    typename S,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T,
    M,
    S>*>
inline auto offset_multiplier_constrain (const T &x , const M &mu , const S &sigma );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:90
template <typename T,
    typename M,
    typename S,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T,
    M,
    S>*>
inline auto offset_multiplier_constrain (const T &x , const M &mu , const S &sigma , return_type_t <T , M , S >&lp );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:121
template <typename T,
    typename M,
    typename S,
    require_all_not_std_vector_t<M,
    S>*>
inline auto offset_multiplier_constrain (const std ::vector <T >&x , const M &mu , const S &sigma );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:140
template <typename T,
    typename M,
    typename S,
    require_all_not_std_vector_t<M,
    S>*>
inline auto offset_multiplier_constrain (const std ::vector <T >&x , const M &mu , const S &sigma , return_type_t <T , M , S >&lp );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:160
template <typename T,
    typename M,
    typename S,
    require_not_std_vector_t<M>*>
inline auto offset_multiplier_constrain (const std ::vector <T >&x , const M &mu , const std ::vector <S >&sigma );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:179
template <typename T,
    typename M,
    typename S,
    require_not_std_vector_t<M>*>
inline auto offset_multiplier_constrain (const std ::vector <T >&x , const M &mu , const std ::vector <S >&sigma , return_type_t <T , M , S >&lp );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:199
template <typename T,
    typename M,
    typename S,
    require_not_std_vector_t<S>*>
inline auto offset_multiplier_constrain (const std ::vector <T >&x , const std ::vector <M >&mu , const S &sigma );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:219
template <typename T,
    typename M,
    typename S,
    require_not_std_vector_t<S>*>
inline auto offset_multiplier_constrain (const std ::vector <T >&x , const std ::vector <M >&mu , const S &sigma , return_type_t <T , M , S >&lp );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:240
template <typename T,
    typename M,
    typename S>
inline auto offset_multiplier_constrain (const std ::vector <T >&x , const std ::vector <M >&mu , const std ::vector <S >&sigma );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:259
template <typename T,
    typename M,
    typename S>
inline auto offset_multiplier_constrain (const std ::vector <T >&x , const std ::vector <M >&mu , const std ::vector <S >&sigma , return_type_t <T , M , S >&lp );

// ./stan/math/prim/constraint/offset_multiplier_constrain.hpp:299
template <bool Jacobian,
    typename T,
    typename M,
    typename S>
inline auto offset_multiplier_constrain (const T &x , const M &mu , const S &sigma , return_type_t <T , M , S >&lp );

// ./stan/math/prim/constraint/offset_multiplier_free.hpp:41
template <typename T,
    typename M,
    typename S>
inline auto offset_multiplier_free (const T &y , const M &mu , const S &sigma );

// ./stan/math/prim/constraint/offset_multiplier_free.hpp:64
template <typename T,
    typename M,
    typename S,
    require_all_not_std_vector_t<M,
    S>*>
inline auto offset_multiplier_free (const std ::vector <T >&x , const M &mu , const S &sigma );

// ./stan/math/prim/constraint/offset_multiplier_free.hpp:82
template <typename T,
    typename M,
    typename S,
    require_not_std_vector_t<M>*>
inline auto offset_multiplier_free (const std ::vector <T >&x , const M &mu , const std ::vector <S >&sigma );

// ./stan/math/prim/constraint/offset_multiplier_free.hpp:101
template <typename T,
    typename M,
    typename S,
    require_not_std_vector_t<S>*>
inline auto offset_multiplier_free (const std ::vector <T >&x , const std ::vector <M >&mu , const S &sigma );

// ./stan/math/prim/constraint/offset_multiplier_free.hpp:120
template <typename T,
    typename M,
    typename S>
inline auto offset_multiplier_free (const std ::vector <T >&x , const std ::vector <M >&mu , const std ::vector <S >&sigma );

// ./stan/math/prim/constraint/ordered_constrain.hpp:24
template <typename EigVec,
    require_eigen_col_vector_t<EigVec>*>
inline plain_type_t <EigVec >ordered_constrain (const EigVec &x );

// ./stan/math/prim/constraint/ordered_constrain.hpp:53
template <typename EigVec,
    require_eigen_col_vector_t<EigVec>*>
inline auto ordered_constrain (const EigVec &x , value_type_t <EigVec >&lp );

// ./stan/math/prim/constraint/ordered_constrain.hpp:80
template <bool Jacobian,
    typename T,
    require_not_std_vector_t<T>*>
inline auto ordered_constrain (const T &x , return_type_t <T >&lp );

// ./stan/math/prim/constraint/ordered_constrain.hpp:107
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto ordered_constrain (const T &x , return_type_t <T >&lp );

// ./stan/math/prim/constraint/ordered_free.hpp:27
template <typename EigVec,
    require_eigen_col_vector_t<EigVec>*>
plain_type_t <EigVec >ordered_free (const EigVec &y );

// ./stan/math/prim/constraint/ordered_free.hpp:51
template <typename T,
    require_std_vector_t<T>*>
auto ordered_free (const T &x );

// ./stan/math/prim/constraint/positive_free.hpp:28
template <typename T>
inline T positive_free (const T &y );

// ./stan/math/prim/constraint/positive_ordered_constrain.hpp:22
template <typename EigVec,
    require_eigen_col_vector_t<EigVec>*>
inline auto positive_ordered_constrain (const EigVec &x );

// ./stan/math/prim/constraint/positive_ordered_constrain.hpp:51
template <typename Vec,
    require_col_vector_t<Vec>*>
inline auto positive_ordered_constrain (const Vec &x , return_type_t <Vec >&lp );

// ./stan/math/prim/constraint/positive_ordered_constrain.hpp:75
template <bool Jacobian,
    typename Vec,
    require_not_std_vector_t<Vec>*>
inline auto positive_ordered_constrain (const Vec &x , return_type_t <Vec >&lp );

// ./stan/math/prim/constraint/positive_ordered_constrain.hpp:102
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto positive_ordered_constrain (const T &x , return_type_t <T >&lp );

// ./stan/math/prim/constraint/positive_ordered_free.hpp:27
template <typename EigVec,
    require_eigen_col_vector_t<EigVec>*>
auto positive_ordered_free (const EigVec &y );

// ./stan/math/prim/constraint/positive_ordered_free.hpp:51
template <typename T,
    require_std_vector_t<T>*>
auto positive_ordered_free (const T &x );

// ./stan/math/prim/constraint/prob_constrain.hpp:25
template <typename T>
inline T prob_constrain (const T &x );

// ./stan/math/prim/constraint/prob_constrain.hpp:50
template <typename T>
inline T prob_constrain (const T &x , T &lp );

// ./stan/math/prim/constraint/prob_constrain.hpp:71
template <bool Jacobian,
    typename T>
inline auto prob_constrain (const T &x , T &lp );

// ./stan/math/prim/constraint/prob_free.hpp:26
template <typename T>
inline T prob_free (const T &y );

// ./stan/math/prim/constraint/simplex_constrain.hpp:27
template <typename Vec,
    require_eigen_vector_t<Vec>*>
inline plain_type_t <Vec >simplex_constrain (const Vec &y );

// ./stan/math/prim/constraint/simplex_constrain.hpp:59
template <typename Vec,
    require_eigen_vector_t<Vec>*>
inline plain_type_t <Vec >simplex_constrain (const Vec &y , value_type_t <Vec >&lp );

// ./stan/math/prim/constraint/simplex_constrain.hpp:101
template <bool Jacobian,
    typename Vec,
    require_not_std_vector_t<Vec>*>
inline plain_type_t <Vec >simplex_constrain (const Vec &y , return_type_t <Vec >&lp );

// ./stan/math/prim/constraint/simplex_constrain.hpp:127
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto simplex_constrain (const T &y , return_type_t <T >&lp );

// ./stan/math/prim/constraint/simplex_free.hpp:29
template <typename Vec,
    require_eigen_vector_t<Vec>*>
inline plain_type_t <Vec >simplex_free (const Vec &x );

// ./stan/math/prim/constraint/simplex_free.hpp:55
template <typename T,
    require_std_vector_t<T>*>
auto simplex_free (const T &x );

// ./stan/math/prim/constraint/stochastic_column_constrain.hpp:25
template <typename Mat,
    require_eigen_matrix_dynamic_t<Mat>*>
inline plain_type_t <Mat >stochastic_column_constrain (const Mat &y );

// ./stan/math/prim/constraint/stochastic_column_constrain.hpp:50
template <typename Mat,
    require_eigen_matrix_dynamic_t<Mat>*>
inline plain_type_t <Mat >stochastic_column_constrain (const Mat &y , value_type_t <Mat >&lp );

// ./stan/math/prim/constraint/stochastic_column_constrain.hpp:77
template <bool Jacobian,
    typename Mat,
    require_not_std_vector_t<Mat>*>
inline plain_type_t <Mat >stochastic_column_constrain (const Mat &y , return_type_t <Mat >&lp );

// ./stan/math/prim/constraint/stochastic_column_constrain.hpp:104
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto stochastic_column_constrain (const T &y , return_type_t <T >&lp );

// ./stan/math/prim/constraint/stochastic_column_free.hpp:20
template <typename Mat,
    require_eigen_matrix_dynamic_t<Mat>*>
inline plain_type_t <Mat >stochastic_column_free (const Mat &y );

// ./stan/math/prim/constraint/stochastic_column_free.hpp:39
template <typename T,
    require_std_vector_t<T>*>
inline auto stochastic_column_free (const T &y );

// ./stan/math/prim/constraint/stochastic_row_constrain.hpp:25
template <typename Mat,
    require_eigen_matrix_dynamic_t<Mat>*>
inline plain_type_t <Mat >stochastic_row_constrain (const Mat &y );

// ./stan/math/prim/constraint/stochastic_row_constrain.hpp:53
template <typename Mat,
    require_eigen_matrix_dynamic_t<Mat>*>
inline plain_type_t <Mat >stochastic_row_constrain (const Mat &y , value_type_t <Mat >&lp );

// ./stan/math/prim/constraint/stochastic_row_constrain.hpp:92
template <bool Jacobian,
    typename Mat,
    require_not_std_vector_t<Mat>*>
inline plain_type_t <Mat >stochastic_row_constrain (const Mat &y , return_type_t <Mat >&lp );

// ./stan/math/prim/constraint/stochastic_row_constrain.hpp:118
template <bool Jacobian,
    typename T,
    require_std_vector_t<T>*>
inline auto stochastic_row_constrain (const T &y , return_type_t <T >&lp );

// ./stan/math/prim/constraint/stochastic_row_free.hpp:19
template <typename Mat,
    require_eigen_matrix_dynamic_t<Mat>*>
inline plain_type_t <Mat >stochastic_row_free (const Mat &y );

// ./stan/math/prim/constraint/stochastic_row_free.hpp:38
template <typename T,
    require_std_vector_t<T>*>
inline auto stochastic_row_free (const T &y );

// ./stan/math/prim/constraint/sum_to_zero_free.hpp:37
template <typename Vec,
    require_eigen_vector_t<Vec>*>
inline plain_type_t <Vec >sum_to_zero_free (const Vec &z );

// ./stan/math/prim/constraint/sum_to_zero_free.hpp:71
template <typename T,
    require_std_vector_t<T>*>
auto sum_to_zero_free (const T &z );

// ./stan/math/prim/constraint/unit_vector_free.hpp:23
template <typename EigVec,
    require_eigen_col_vector_t<EigVec>*>
inline auto unit_vector_free (EigVec &&x );

// ./stan/math/prim/constraint/unit_vector_free.hpp:38
template <typename T,
    require_std_vector_t<T>*>
auto unit_vector_free (const T &x );

// ./stan/math/rev/fun/Phi.hpp:53

inline var Phi (const var &a );

// ./stan/math/rev/fun/Phi.hpp:66
template <typename T,
    require_var_matrix_t<T>*>
inline auto Phi (const T &a );

// ./stan/math/rev/fun/Phi_approx.hpp:46

inline var Phi_approx (const var &a );

// ./stan/math/rev/fun/Phi_approx.hpp:54
template <typename T,
    require_var_matrix_t<T>*>
inline auto Phi_approx (const T &a );

// ./stan/math/rev/fun/fabs.hpp:50

inline var fabs (const var &a );

// ./stan/math/rev/fun/fabs.hpp:70
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto fabs (const VarMat &a );

// ./stan/math/rev/fun/hypot.hpp:27

inline var hypot (const var &a , const var &b );

// ./stan/math/rev/fun/hypot.hpp:46

inline var hypot (const var &a , double b );

// ./stan/math/rev/fun/hypot.hpp:91

inline var hypot (double a , const var &b );

// ./stan/math/rev/fun/abs.hpp:40
template <typename T>
inline auto abs (const var_value <T >&a );

// ./stan/math/rev/fun/abs.hpp:51

inline var abs (const std ::complex <var >&z );

// ./stan/math/rev/fun/atan2.hpp:27

inline var atan2 (const var &a , const var &b );

// ./stan/math/rev/fun/atan2.hpp:48

inline var atan2 (const var &a , double b );

// ./stan/math/rev/fun/atan2.hpp:92

inline var atan2 (double a , const var &b );

// ./stan/math/rev/fun/atan2.hpp:100
template <typename Mat1,
    typename Mat2,
    require_any_var_matrix_t<Mat1,
    Mat2>*>
inline auto atan2 (const Mat1 &a , const Mat2 &b );

// ./stan/math/rev/fun/atan2.hpp:148
template <typename Scalar,
    typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto atan2 (const Scalar &a , const VarMat &b );

// ./stan/math/rev/fun/atan2.hpp:194
template <typename VarMat,
    typename Scalar,
    require_var_matrix_t<VarMat>*>
inline auto atan2 (const VarMat &a , const Scalar &b );

// ./stan/math/rev/fun/arg.hpp:18

inline var arg (const std ::complex <var >&z );

// ./stan/math/rev/fun/value_of_rec.hpp:16
template <typename T>
inline auto &value_of_rec (const var_value <T >&v );

// ./stan/math/rev/fun/is_inf.hpp:20

inline int is_inf (const var &v );

// ./stan/math/rev/fun/is_nan.hpp:20

inline bool is_nan (const var &v );

// ./stan/math/rev/fun/sinh.hpp:44

inline var sinh (const var &a );

// ./stan/math/rev/fun/sinh.hpp:58
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto sinh (const VarMat &a );

// ./stan/math/rev/fun/sinh.hpp:72

inline std ::complex <var >sinh (const std ::complex <var >&z );

// ./stan/math/rev/fun/cos.hpp:46

inline var cos (var a );

// ./stan/math/rev/fun/cos.hpp:60
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto cos (const VarMat &a );

// ./stan/math/rev/fun/cos.hpp:73

inline std ::complex <var >cos (const std ::complex <var >&z );

// ./stan/math/rev/fun/sin.hpp:46

inline var sin (const var &a );

// ./stan/math/rev/fun/sin.hpp:59
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto sin (const VarMat &a );

// ./stan/math/rev/fun/sin.hpp:73

inline std ::complex <var >sin (const std ::complex <var >&z );

// ./stan/math/rev/fun/exp.hpp:42

inline var exp (const var &a );

// ./stan/math/rev/fun/exp.hpp:53

inline std ::complex <var >exp (const std ::complex <var >&z );

// ./stan/math/rev/fun/exp.hpp:64
template <typename T,
    require_var_matrix_t<T>*>
inline auto exp (const T &x );

// ./stan/math/rev/fun/cosh.hpp:44

inline var cosh (const var &a );

// ./stan/math/rev/fun/cosh.hpp:57
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto cosh (const VarMat &a );

// ./stan/math/rev/fun/cosh.hpp:71

inline std ::complex <var >cosh (const std ::complex <var >&z );

// ./stan/math/rev/fun/square.hpp:35

inline var square (const var &x );

// ./stan/math/rev/fun/square.hpp:48
template <typename T,
    require_var_matrix_t<T>*>
inline auto square (const T &x );

// ./stan/math/rev/fun/norm.hpp:18

inline var norm (const std ::complex <var >&z );

// ./stan/math/rev/fun/sqrt.hpp:44

inline var sqrt (const var &a );

// ./stan/math/rev/fun/sqrt.hpp:59
template <typename T,
    require_var_matrix_t<T>*>
inline auto sqrt (const T &a );

// ./stan/math/rev/fun/sqrt.hpp:75

inline std ::complex <var >sqrt (const std ::complex <var >&z );

// ./stan/math/rev/fun/log.hpp:50
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto log (const var_value <T >&a );

// ./stan/math/rev/fun/log.hpp:64

inline std ::complex <var >log (const std ::complex <var >&z );

// ./stan/math/rev/fun/polar.hpp:26

inline std ::complex <var >polar (const var &r , const var &theta );

// ./stan/math/rev/fun/polar.hpp:38
template <typename T>
inline std ::complex <var >polar (T r , const var &theta );

// ./stan/math/rev/fun/polar.hpp:51
template <typename T>
inline std ::complex <var >polar (const var &r , T theta );

// ./stan/math/rev/fun/asinh.hpp:59

inline var asinh (const var &x );

// ./stan/math/rev/fun/asinh.hpp:72
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto asinh (const VarMat &x );

// ./stan/math/rev/fun/asinh.hpp:88

inline std ::complex <var >asinh (const std ::complex <var >&z );

// ./stan/math/rev/fun/asin.hpp:53

inline var asin (const var &x );

// ./stan/math/rev/fun/asin.hpp:67
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto asin (const VarMat &x );

// ./stan/math/rev/fun/asin.hpp:82

inline std ::complex <var >asin (const std ::complex <var >&z );

// ./stan/math/rev/fun/acos.hpp:55

inline var acos (const var &x );

// ./stan/math/rev/fun/acos.hpp:68
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto acos (const VarMat &x );

// ./stan/math/rev/fun/acos.hpp:83

inline std ::complex <var >acos (const std ::complex <var >&x );

// ./stan/math/rev/fun/acosh.hpp:63

inline var acosh (const var &x );

// ./stan/math/rev/fun/acosh.hpp:77
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto acosh (const VarMat &x );

// ./stan/math/rev/fun/acosh.hpp:93

inline std ::complex <var >acosh (const std ::complex <var >&z );

// ./stan/math/rev/fun/append_col.hpp:34
template <typename T1,
    typename T2,
    require_any_var_matrix_t<T1,
    T2>*>
inline auto append_col (const T1 &A , const T2 &B );

// ./stan/math/rev/fun/append_col.hpp:78
template <typename Scal,
    typename RowVec,
    require_stan_scalar_t<Scal>*
;

// ./stan/math/rev/fun/append_col.hpp:118
template <typename RowVec,
    typename Scal,
    require_t<is_eigen_row_vector<RowVec>>*>
inline auto append_col (const var_value <RowVec >&A , const Scal &B );

// ./stan/math/rev/fun/append_row.hpp:32
template <typename T1,
    typename T2,
    require_any_var_matrix_t<T1,
    T2>*>
inline auto append_row (const T1 &A , const T2 &B );

// ./stan/math/rev/fun/append_row.hpp:75
template <typename Scal,
    typename ColVec,
    require_stan_scalar_t<Scal>*
;

// ./stan/math/rev/fun/append_row.hpp:114
template <typename ColVec,
    typename Scal,
    require_t<is_eigen_col_vector<ColVec>>*>
inline auto append_row (const var_value <ColVec >&A , const Scal &B );

// ./stan/math/rev/fun/as_bool.hpp:16

inline int as_bool (const var &v );

// ./stan/math/rev/fun/as_array_or_scalar.hpp:19
template <typename T,
    require_var_matrix_t<T>*>
inline auto as_array_or_scalar (T &&v );

// ./stan/math/rev/fun/as_column_vector_or_scalar.hpp:18
template <typename T,
    require_var_col_vector_t<T>*>
inline auto as_column_vector_or_scalar (T &&a );

// ./stan/math/rev/fun/as_column_vector_or_scalar.hpp:30
template <typename T,
    require_var_row_vector_t<T>*>
inline auto as_column_vector_or_scalar (T &&a );

// ./stan/math/rev/fun/atanh.hpp:58

inline var atanh (const var &x );

// ./stan/math/rev/fun/atanh.hpp:72
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto atanh (const VarMat &x );

// ./stan/math/rev/fun/atanh.hpp:87

inline std ::complex <var >atanh (const std ::complex <var >&z );

// ./stan/math/rev/fun/atan.hpp:54

inline var atan (const var &x );

// ./stan/math/rev/fun/atan.hpp:69
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto atan (const VarMat &x );

// ./stan/math/rev/fun/atan.hpp:83

inline std ::complex <var >atan (const std ::complex <var >&z );

// ./stan/math/rev/fun/bessel_first_kind.hpp:11

inline var bessel_first_kind (int v , const var &a );

// ./stan/math/rev/fun/bessel_first_kind.hpp:25
template <typename T1,
    typename T2,
    require_st_integral<T1>*>
inline auto bessel_first_kind (const T1 &v , const var_value <T2 >&a );

// ./stan/math/rev/fun/bessel_second_kind.hpp:11

inline var bessel_second_kind (int v , const var &a );

// ./stan/math/rev/fun/bessel_second_kind.hpp:24
template <typename T1,
    typename T2,
    require_st_integral<T1>*>
inline auto bessel_second_kind (const T1 &v , const var_value <T2 >&a );

// ./stan/math/rev/fun/beta.hpp:37

inline var beta (const var &a , const var &b );

// ./stan/math/rev/fun/beta.hpp:67

inline var beta (const var &a , double b );

// ./stan/math/rev/fun/beta.hpp:92

inline var beta (double a , const var &b );

// ./stan/math/rev/fun/beta.hpp:100
template <typename Mat1,
    typename Mat2,
    require_any_var_matrix_t<Mat1,
    Mat2>*>
inline auto beta (const Mat1 &a , const Mat2 &b );

// ./stan/math/rev/fun/beta.hpp:146
template <typename Scalar,
    typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto beta (const Scalar &a , const VarMat &b );

// ./stan/math/rev/fun/beta.hpp:188
template <typename VarMat,
    typename Scalar,
    require_var_matrix_t<VarMat>*>
inline auto beta (const VarMat &a , const Scalar &b );

// ./stan/math/rev/fun/binary_log_loss.hpp:46

inline var binary_log_loss (int y , const var &y_hat );

// ./stan/math/rev/fun/binary_log_loss.hpp:61
template <typename Mat,
    require_eigen_t<Mat>*>
inline auto binary_log_loss (int y , const var_value <Mat >&y_hat );

// ./stan/math/rev/fun/binary_log_loss.hpp:80
template <typename StdVec,
    typename Mat,
    require_eigen_t<Mat>*>
inline auto binary_log_loss (const StdVec &y , const var_value <Mat >&y_hat );

// ./stan/math/rev/fun/cbrt.hpp:37

inline var cbrt (const var &a );

// ./stan/math/rev/fun/cbrt.hpp:49
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto cbrt (const VarMat &a );

// ./stan/math/rev/fun/ceil.hpp:47

inline var ceil (const var &a );

// ./stan/math/rev/fun/ceil.hpp:49
template <typename T,
    require_matrix_t<T>*>
inline auto ceil (const var_value <T >&a );

// ./stan/math/rev/fun/cholesky_decompose.hpp:134
template <typename EigMat,
    require_eigen_vt<is_var,
    EigMat>*>
inline auto cholesky_decompose (const EigMat &A );

// ./stan/math/rev/fun/cholesky_decompose.hpp:170
template <typename T,
    require_var_matrix_t<T>*>
inline auto cholesky_decompose (const T &A );

// ./stan/math/rev/fun/cumulative_sum.hpp:29
template <typename EigVec,
    require_rev_vector_t<EigVec>*>
inline auto cumulative_sum (const EigVec &x );

// ./stan/math/rev/fun/columns_dot_product.hpp:31
template <typename Mat1,
    typename Mat2,
    require_all_eigen_t<Mat1,
    Mat2>*>
inline Eigen ::Matrix <return_type_t <Mat1 ,Mat2 >,1,Mat1 ::ColsAtCompileTime >columns_dot_product (const Mat1 &v1 , const Mat2 &v2 );

// ./stan/math/rev/fun/columns_dot_product.hpp:61
template <typename Mat1,
    typename Mat2,
    require_all_matrix_t<Mat1,
    Mat2>*>
inline auto columns_dot_product (const Mat1 &v1 , const Mat2 &v2 );

// ./stan/math/rev/fun/columns_dot_self.hpp:19
template <typename Mat,
    require_eigen_vt<is_var,
    Mat>*>
inline Eigen ::Matrix <var ,1,Mat ::ColsAtCompileTime >columns_dot_self (const Mat &x );

// ./stan/math/rev/fun/columns_dot_self.hpp:35
template <typename Mat,
    require_var_matrix_t<Mat>*>
inline auto columns_dot_self (const Mat &x );

// ./stan/math/rev/fun/conj.hpp:17

inline std ::complex <var >conj (const std ::complex <var >&z );

// ./stan/math/rev/fun/adjoint_of.hpp:36
template <typename T,
    require_var_t<T>*>
auto &adjoint_of (const T &x );

// ./stan/math/rev/fun/adjoint_of.hpp:49
template <typename T,
    require_not_var_t<T>*>
internal ::nonexisting_adjoint adjoint_of (const T &x );

// ./stan/math/rev/fun/gp_exp_quad_cov.hpp:31
template <typename T_x,
    typename T_sigma,
    require_st_arithmetic<T_x>*>
inline Eigen ::Matrix <var ,-1,-1>gp_exp_quad_cov (const std ::vector <T_x >&x , const T_sigma sigma , const var length_scale );

// ./stan/math/rev/fun/cov_exp_quad.hpp:152
template <typename T_x,
    typename >
inline Eigen ::Matrix <var ,-1,-1>cov_exp_quad (const std ::vector <T_x >&x , const var &sigma , const var &l );

// ./stan/math/rev/fun/cov_exp_quad.hpp:162
template <typename T_x,
    typename >
inline Eigen ::Matrix <var ,-1,-1>cov_exp_quad (const std ::vector <T_x >&x , double sigma , const var &l );

// ./stan/math/rev/fun/to_soa_sparse_matrix.hpp:28
template <int Options
;

// ./stan/math/rev/fun/to_soa_sparse_matrix.hpp:62
template <int Options>
inline auto to_soa_sparse_matrix (int m , int n , MatrixVar &&w , Vec1 &&u , Vec2 &&v );

// ./stan/math/rev/fun/to_soa_sparse_matrix.hpp:109
template <int Options>
inline auto to_soa_sparse_matrix (int m , int n , Mat &&w , Vec1 &&u , Vec2 &&v );

// ./stan/math/rev/fun/csr_matrix_times_vector.hpp:157
template <typename T1,
    typename T2,
    require_any_rev_matrix_t<T1,
    T2>*>
inline auto csr_matrix_times_vector (int m , int n , const T1 &w , const std ::vector <int >&v , const std ::vector <int >&u , const T2 &b );

// ./stan/math/rev/fun/determinant.hpp:12
template <typename T,
    require_rev_matrix_t<T>*>
inline var determinant (const T &m );

// ./stan/math/rev/fun/diag_pre_multiply.hpp:23
template <typename T1,
    typename T2,
    require_vector_t<T1>*>
auto diag_pre_multiply (const T1 &m1 , const T2 &m2 );

// ./stan/math/rev/fun/diag_post_multiply.hpp:23
template <typename T1,
    typename T2,
    require_matrix_t<T1>*>
auto diag_post_multiply (const T1 &m1 , const T2 &m2 );

// ./stan/math/rev/fun/digamma.hpp:19

inline var digamma (const var &a );

// ./stan/math/rev/fun/digamma.hpp:33
template <typename T,
    require_var_matrix_t<T>*>
inline auto digamma (const T &a );

// ./stan/math/rev/fun/eigendecompose_sym.hpp:27
template <typename T,
    require_rev_matrix_t<T>*>
inline auto eigendecompose_sym (const T &m );

// ./stan/math/rev/fun/eigenvalues_sym.hpp:23
template <typename T,
    require_rev_matrix_t<T>*>
inline auto eigenvalues_sym (const T &m );

// ./stan/math/rev/fun/eigenvectors_sym.hpp:24
template <typename T,
    require_rev_matrix_t<T>*>
inline auto eigenvectors_sym (const T &m );

// ./stan/math/rev/fun/elt_divide.hpp:24
template <typename Mat1,
    typename Mat2,
    require_all_matrix_t<Mat1,
    Mat2>*>
auto elt_divide (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/rev/fun/elt_divide.hpp:78
template <typename Scal,
    typename Mat,
    require_stan_scalar_t<Scal>*>
auto elt_divide (Scal s , const Mat &m );

// ./stan/math/rev/fun/multiply.hpp:26
template <typename T1,
    typename T2,
    require_all_matrix_t<T1,
    T2>*>
inline auto multiply (T1 &&A , T2 &&B );

// ./stan/math/rev/fun/multiply.hpp:87
template <typename T1,
    typename T2,
    require_all_matrix_t<T1,
    T2>*>
inline var multiply (const T1 &A , const T2 &B );

// ./stan/math/rev/fun/multiply.hpp:136
template <typename T1,
    typename T2,
    require_not_matrix_t<T1>*>
inline auto multiply (const T1 &a , T2 &&B );

// ./stan/math/rev/fun/multiply.hpp:189
template <typename T1,
    typename T2,
    require_matrix_t<T1>*
;

// ./stan/math/rev/fun/multiply.hpp:209
template <typename T1,
    typename T2,
    require_any_var_matrix_t<T1,
    T2>*>
inline auto operator *(T1 &&a , T2 &&b );

// ./stan/math/rev/fun/elt_multiply.hpp:25
template <typename Mat1,
    typename Mat2,
    require_all_matrix_t<Mat1,
    Mat2>*>
auto elt_multiply (const Mat1 &m1 , const Mat2 &m2 );

// ./stan/math/rev/fun/erf.hpp:48

inline var erf (const var &a );

// ./stan/math/rev/fun/erf.hpp:55
template <typename T,
    require_matrix_t<T>*>
inline auto erf (const var_value <T >&a );

// ./stan/math/rev/fun/erfc.hpp:48

inline var erfc (const var &a );

// ./stan/math/rev/fun/erfc.hpp:55
template <typename T,
    require_matrix_t<T>*>
inline auto erfc (const var_value <T >&a );

// ./stan/math/rev/fun/exp2.hpp:39

inline var exp2 (const var &a );

// ./stan/math/rev/fun/exp2.hpp:45
template <typename T,
    require_eigen_t<T>*>
inline auto exp2 (const var_value <T >&a );

// ./stan/math/rev/fun/expm1.hpp:38

inline var expm1 (const var &a );

// ./stan/math/rev/fun/expm1.hpp:44
template <typename T,
    require_eigen_t<T>*>
inline auto expm1 (const var_value <T >&a );

// ./stan/math/rev/fun/falling_factorial.hpp:12

inline var falling_factorial (const var &a , int b );

// ./stan/math/rev/fun/falling_factorial.hpp:20
template <typename T1,
    typename T2,
    require_eigen_t<T1>*>
inline auto falling_factorial (const var_value <T1 >&a , const T2 &b );

// ./stan/math/rev/fun/fdim.hpp:91

inline var fdim (const var &a , const var &b );

// ./stan/math/rev/fun/fdim.hpp:110

inline var fdim (double a , const var &b );

// ./stan/math/rev/fun/fdim.hpp:127

inline var fdim (const var &a , double b );

// ./stan/math/rev/fun/fft.hpp:40
template <typename V,
    require_eigen_vector_vt<is_complex,
    V>*
;

// ./stan/math/rev/fun/fft.hpp:85
template <typename V,
    require_eigen_vector_vt<is_complex,
    V>*
;

// ./stan/math/rev/fun/fft.hpp:121
template <typename M,
    require_eigen_dense_dynamic_vt<is_complex,
    M>*
;

// ./stan/math/rev/fun/fft.hpp:153
template <typename M,
    require_eigen_dense_dynamic_vt<is_complex,
    M>*
;

// ./stan/math/rev/fun/fill.hpp:23
template <typename VarMat,
    typename S,
    require_var_matrix_t<VarMat>*>
inline void fill (VarMat &x , const S &y );

// ./stan/math/rev/fun/fill.hpp:46
template <typename VarMat,
    typename S,
    require_var_matrix_t<VarMat>*>
inline void fill (VarMat &x , const S &y );

// ./stan/math/rev/fun/floor.hpp:47

inline var floor (const var &a );

// ./stan/math/rev/fun/floor.hpp:49
template <typename T,
    require_eigen_t<T>*>
inline auto floor (const var_value <T >&a );

// ./stan/math/rev/fun/fma.hpp:31

inline var fma (const var &x , const var &y , const var &z );

// ./stan/math/rev/fun/fma.hpp:56
template <typename Tc,
    require_arithmetic_t<Tc>*>
inline var fma (const var &x , const var &y , Tc &&z );

// ./stan/math/rev/fun/fma.hpp:84
template <typename Tb,
    require_arithmetic_t<Tb>*>
inline var fma (const var &x , Tb &&y , const var &z );

// ./stan/math/rev/fun/fma.hpp:113
template <typename Tb,
    typename Tc,
    require_all_arithmetic_t<Tb,
    Tc>*>
inline var fma (const var &x , Tb &&y , Tc &&z );

// ./stan/math/rev/fun/fma.hpp:136
template <typename Ta,
    typename Tc,
    require_all_arithmetic_t<Ta,
    Tc>*>
inline var fma (Ta &&x , const var &y , Tc &&z );

// ./stan/math/rev/fun/fma.hpp:159
template <typename Ta,
    typename Tb,
    require_all_arithmetic_t<Ta,
    Tb>*>
inline var fma (Ta &&x , Tb &&y , const var &z );

// ./stan/math/rev/fun/fma.hpp:182
template <typename Ta,
    require_arithmetic_t<Ta>*>
inline var fma (Ta &&x , const var &y , const var &z );

// ./stan/math/rev/fun/fma.hpp:385
template <typename T1,
    typename T2,
    typename T3,
    require_any_matrix_t<T1,
    T2,
    T3>*
;

// ./stan/math/rev/fun/fmax.hpp:61

inline var fmax (const var &a , const var &b );

// ./stan/math/rev/fun/fmax.hpp:91

inline var fmax (const var &a , double b );

// ./stan/math/rev/fun/fmax.hpp:119

inline var fmax (double a , const var &b );

// ./stan/math/rev/fun/fmin.hpp:57

inline var fmin (const var &a , const var &b );

// ./stan/math/rev/fun/fmin.hpp:85

inline var fmin (const var &a , double b );

// ./stan/math/rev/fun/fmin.hpp:113

inline var fmin (double a , const var &b );

// ./stan/math/rev/fun/fmod.hpp:100

inline var fmod (const var &a , const var &b );

// ./stan/math/rev/fun/fmod.hpp:117

inline var fmod (const var &a , double b );

// ./stan/math/rev/fun/fmod.hpp:134

inline var fmod (double a , const var &b );

// ./stan/math/rev/fun/from_var_value.hpp:19
template <typename T,
    require_var_matrix_t<T>*>
Eigen ::Matrix <var ,T ::RowsAtCompileTime ,T ::ColsAtCompileTime >from_var_value (const T &a );

// ./stan/math/rev/fun/from_var_value.hpp:51
template <typename T>
auto from_var_value (const std ::vector <T >&a );

// ./stan/math/rev/fun/gamma_p.hpp:102

inline var gamma_p (const var &a , const var &b );

// ./stan/math/rev/fun/gamma_p.hpp:106

inline var gamma_p (const var &a , double b );

// ./stan/math/rev/fun/gamma_p.hpp:110

inline var gamma_p (double a , const var &b );

// ./stan/math/rev/fun/gamma_q.hpp:50

inline var gamma_q (const var &a , const var &b );

// ./stan/math/rev/fun/gamma_q.hpp:54

inline var gamma_q (const var &a , double b );

// ./stan/math/rev/fun/gamma_q.hpp:58

inline var gamma_q (double a , const var &b );

// ./stan/math/rev/fun/inverse.hpp:22
template <typename T,
    require_rev_matrix_t<T>*>
inline auto inverse (const T &m );

// ./stan/math/rev/fun/generalized_inverse.hpp:63
template <typename VarMat,
    require_rev_matrix_t<VarMat>*>
inline auto generalized_inverse (const VarMat &G );

// ./stan/math/rev/fun/gp_periodic_cov.hpp:40
template <typename T_x,
    typename T_sigma,
    require_st_arithmetic<T_x>*>
inline Eigen ::Matrix <var ,Eigen ::Dynamic ,Eigen ::Dynamic >gp_periodic_cov (const std ::vector <T_x >&x , const T_sigma sigma , const var l , const var p );

// ./stan/math/rev/fun/grad.hpp:24

inline void grad (var &v , Eigen ::Matrix <var , Eigen ::Dynamic , 1>&x , Eigen ::VectorXd &g );

// ./stan/math/rev/fun/inv.hpp:30
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto inv (const var_value <T >&a );

// ./stan/math/rev/fun/inv_sqrt.hpp:33
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto inv_sqrt (const var_value <T >&a );

// ./stan/math/rev/fun/inv_square.hpp:30

inline var inv_square (const var &a );

// ./stan/math/rev/fun/pow.hpp:68
template <typename Scal1,
    typename Scal2,
    require_any_st_var<Scal1,
    Scal2>*>
inline var pow (const Scal1 &base , const Scal2 &exponent );

// ./stan/math/rev/fun/pow.hpp:118
template <typename Mat1,
    typename Mat2,
    require_all_st_var_or_arithmetic<Mat1,
    Mat2>*>
inline auto pow (const Mat1 &base , const Mat2 &exponent );

// ./stan/math/rev/fun/pow.hpp:171
template <typename Mat1,
    typename Scal1,
    require_all_st_var_or_arithmetic<Mat1,
    Scal1>*>
inline auto pow (const Mat1 &base , const Scal1 &exponent );

// ./stan/math/rev/fun/pow.hpp:236
template <typename Scal1,
    typename Mat1,
    require_all_st_var_or_arithmetic<Scal1,
    Mat1>*>
inline auto pow (Scal1 base , const Mat1 &exponent );

// ./stan/math/rev/fun/pow.hpp:283

inline std ::complex <var >pow (const std ::complex <var >&x , const std ::complex <var >&y );

// ./stan/math/rev/fun/pow.hpp:296
template <typename T,
    typename >
inline std ::complex <var >pow (const std ::complex <var >&x , const std ::complex <T >y );

// ./stan/math/rev/fun/pow.hpp:309

inline std ::complex <var >pow (const std ::complex <var >&x , const var &y );

// ./stan/math/rev/fun/pow.hpp:321
template <typename T,
    typename >
inline std ::complex <var >pow (const std ::complex <var >&x , T y );

// ./stan/math/rev/fun/pow.hpp:334
template <typename T,
    typename >
inline std ::complex <var >pow (std ::complex <T >x , const std ::complex <var >&y );

// ./stan/math/rev/fun/pow.hpp:347
template <typename T,
    typename >
inline std ::complex <var >pow (std ::complex <T >x , const var &y );

// ./stan/math/rev/fun/pow.hpp:359

inline std ::complex <var >pow (const var &x , const std ::complex <var >&y );

// ./stan/math/rev/fun/pow.hpp:371
template <typename T,
    typename >
inline std ::complex <var >pow (const var &x , std ::complex <T >y );

// ./stan/math/rev/fun/pow.hpp:384
template <typename T,
    typename >
inline std ::complex <var >pow (T x , const std ::complex <var >&y );

// ./stan/math/rev/fun/pow.hpp:400

inline std ::complex <var >pow (const std ::complex <var >&x , int y );

// ./stan/math/rev/fun/inc_beta.hpp:38

inline var inc_beta (const var &a , const var &b , const var &c );

// ./stan/math/rev/fun/lgamma.hpp:23
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto lgamma (const var_value <T >&a );

// ./stan/math/rev/fun/log1m.hpp:22
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto log1m (const var_value <T >&a );

// ./stan/math/rev/fun/hypergeometric_2F1.hpp:29
template <typename Ta1,
    typename Ta2,
    typename Tb,
    typename Tz,
    require_all_stan_scalar_t<Ta1,
    Ta2,
    Tb,
    Tz>*>
inline return_type_t <Ta1 ,Ta1 ,Tb ,Tz >hypergeometric_2F1 (const Ta1 &a1 , const Ta2 &a2 , const Tb &b , const Tz &z );

// ./stan/math/rev/fun/grad_inc_beta.hpp:37

inline void grad_inc_beta (var &g1 , var &g2 , const var &a , const var &b , const var &z );

// ./stan/math/rev/fun/hypergeometric_1F0.hpp:31
template <typename Ta,
    typename Tz,
    require_all_stan_scalar_t<Ta,
    Tz>*>
var hypergeometric_1f0 (const Ta &a , const Tz &z );

// ./stan/math/rev/fun/hypergeometric_pFq.hpp:24
template <typename Ta,
    typename Tb,
    typename Tz,
    bool grad_a>
inline var hypergeometric_pFq (const Ta &a , const Tb &b , const Tz &z );

// ./stan/math/rev/fun/if_else.hpp:18

inline var if_else (bool c , const var &y_true , const var &y_false );

// ./stan/math/rev/fun/if_else.hpp:30

inline var if_else (bool c , double y_true , const var &y_false );

// ./stan/math/rev/fun/if_else.hpp:46

inline var if_else (bool c , const var &y_true , double y_false );

// ./stan/math/rev/fun/inv_inc_beta.hpp:54
template <typename T1,
    typename T2,
    typename T3,
    require_all_stan_scalar_t<T1,
    T2,
    T3>*>
inline var inv_inc_beta (const T1 &a , const T2 &b , const T3 &p );

// ./stan/math/rev/fun/initialize_fill.hpp:24
template <typename VarMat,
    typename S,
    require_var_matrix_t<VarMat>*>
inline void initialize_fill (VarMat &x , const S &y );

// ./stan/math/rev/fun/initialize_variable.hpp:16

inline void initialize_variable (var &variable , const var &value );

// ./stan/math/rev/fun/initialize_variable.hpp:26
template <int R,
    int C>
inline void initialize_variable (Eigen ::Matrix <var , R , C >&matrix , const var &value );

// ./stan/math/rev/fun/inv_Phi.hpp:21

inline var inv_Phi (const var &p );

// ./stan/math/rev/fun/inv_Phi.hpp:34
template <typename T,
    require_var_matrix_t<T>*>
inline auto inv_Phi (const T &p );

// ./stan/math/rev/fun/inv_cloglog.hpp:27
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto inv_cloglog (const var_value <T >&a );

// ./stan/math/rev/fun/inv_erfc.hpp:28

inline var inv_erfc (const var &a );

// ./stan/math/rev/fun/inv_erfc.hpp:37
template <typename T,
    require_matrix_t<T>*>
inline auto inv_erfc (const var_value <T >&a );

// ./stan/math/rev/fun/inv_logit.hpp:25
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto inv_logit (const var_value <T >&a );

// ./stan/math/rev/fun/is_uninitialized.hpp:23

inline bool is_uninitialized (var x );

// ./stan/math/rev/fun/lambert_w.hpp:21
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto lambert_w0 (const var_value <T >&a );

// ./stan/math/rev/fun/lambert_w.hpp:39
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto lambert_wm1 (const var_value <T >&a );

// ./stan/math/rev/fun/lbeta.hpp:64

inline var lbeta (const var &a , const var &b );

// ./stan/math/rev/fun/lbeta.hpp:83

inline var lbeta (const var &a , double b );

// ./stan/math/rev/fun/lbeta.hpp:102

inline var lbeta (double a , const var &b );

// ./stan/math/rev/fun/ldexp.hpp:19

inline var ldexp (const var &a , int b );

// ./stan/math/rev/fun/lmgamma.hpp:27

inline var lmgamma (int a , const var &b );

// ./stan/math/rev/fun/lmultiply.hpp:57

inline var lmultiply (const var &a , const var &b );

// ./stan/math/rev/fun/lmultiply.hpp:70

inline var lmultiply (const var &a , double b );

// ./stan/math/rev/fun/lmultiply.hpp:83

inline var lmultiply (double a , const var &b );

// ./stan/math/rev/fun/lmultiply.hpp:102
template <typename T1,
    typename T2,
    require_all_matrix_t<T1,
    T2>*>
inline auto lmultiply (const T1 &a , const T2 &b );

// ./stan/math/rev/fun/lmultiply.hpp:150
template <typename T1,
    typename T2,
    require_var_matrix_t<T1>*>
inline auto lmultiply (const T1 &a , const T2 &b );

// ./stan/math/rev/fun/lmultiply.hpp:196
template <typename T1,
    typename T2,
    require_stan_scalar_t<T1>*>
inline auto lmultiply (const T1 &a , const T2 &b );

// ./stan/math/rev/fun/log10.hpp:48
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto log10 (const var_value <T >&a );

// ./stan/math/rev/fun/log10.hpp:62

inline std ::complex <var >log10 (const std ::complex <var >&z );

// ./stan/math/rev/fun/log1m_exp.hpp:24
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto log1m_exp (const var_value <T >&x );

// ./stan/math/rev/fun/log1m_inv_logit.hpp:20
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto log1m_inv_logit (const var_value <T >&u );

// ./stan/math/rev/fun/log1p.hpp:22
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto log1p (const var_value <T >&a );

// ./stan/math/rev/fun/log1p_exp.hpp:18
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto log1p_exp (const var_value <T >&a );

// ./stan/math/rev/fun/log2.hpp:43
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto log2 (const var_value <T >&a );

// ./stan/math/rev/fun/log_determinant.hpp:13
template <typename T,
    require_rev_matrix_t<T>*>
inline var log_determinant (const T &m );

// ./stan/math/rev/fun/log_determinant_ldlt.hpp:20
template <typename T,
    require_rev_matrix_t<T>*>
var log_determinant_ldlt (LDLT_factor <T >&A );

// ./stan/math/rev/fun/log_determinant_spd.hpp:23
template <typename T,
    require_rev_matrix_t<T>*>
inline var log_determinant_spd (const T &M );

// ./stan/math/rev/fun/log_diff_exp.hpp:56

inline var log_diff_exp (const var &a , const var &b );

// ./stan/math/rev/fun/log_diff_exp.hpp:67

inline var log_diff_exp (const var &a , double b );

// ./stan/math/rev/fun/log_diff_exp.hpp:78

inline var log_diff_exp (double a , const var &b );

// ./stan/math/rev/fun/log_falling_factorial.hpp:61

inline var log_falling_factorial (const var &a , double b );

// ./stan/math/rev/fun/log_falling_factorial.hpp:65

inline var log_falling_factorial (const var &a , const var &b );

// ./stan/math/rev/fun/log_falling_factorial.hpp:69

inline var log_falling_factorial (double a , const var &b );

// ./stan/math/rev/fun/log_inv_logit.hpp:20
template <typename T,
    require_arithmetic_t<T>*>
inline auto log_inv_logit (const var_value <T >&u );

// ./stan/math/rev/fun/log_inv_logit.hpp:35
template <typename T,
    require_eigen_t<T>*>
inline auto log_inv_logit (const var_value <T >&u );

// ./stan/math/rev/fun/log_inv_logit_diff.hpp:71

inline var log_inv_logit_diff (const var &a , double b );

// ./stan/math/rev/fun/log_inv_logit_diff.hpp:75

inline var log_inv_logit_diff (const var &a , const var &b );

// ./stan/math/rev/fun/log_inv_logit_diff.hpp:79

inline var log_inv_logit_diff (double a , const var &b );

// ./stan/math/rev/fun/log_mix.hpp:24

inline void log_mix_partial_helper (double theta_val , double lambda1_val , double lambda2_val , double &one_m_exp_lam2_m_lam1 , double &one_m_t_prod_exp_lam2_m_lam1 , double &one_d_t_plus_one_m_t_prod_exp_lam2_m_lam1 );

// ./stan/math/rev/fun/log_mix.hpp:77
template <typename T_theta,
    typename T_lambda1,
    typename T_lambda2,
    require_any_var_t<T_theta,
    T_lambda1,
    T_lambda2>*>
inline return_type_t <T_theta ,T_lambda1 ,T_lambda2 >log_mix (const T_theta &theta , const T_lambda1 &lambda1 , const T_lambda2 &lambda2 );

// ./stan/math/rev/fun/log_rising_factorial.hpp:43

inline var log_rising_factorial (const var &a , double b );

// ./stan/math/rev/fun/log_rising_factorial.hpp:47

inline var log_rising_factorial (const var &a , const var &b );

// ./stan/math/rev/fun/log_rising_factorial.hpp:51

inline var log_rising_factorial (double a , const var &b );

// ./stan/math/rev/fun/log_softmax.hpp:56
template <typename T,
    require_eigen_st<is_var,
    T>*>
auto log_softmax (const T &x );

// ./stan/math/rev/fun/log_softmax.hpp:99
template <typename T,
    require_var_matrix_t<T>*>
inline auto log_softmax (const T &x );

// ./stan/math/rev/fun/log_softmax.hpp:122
template <typename T,
    require_std_vector_st<is_var,
    T>*>
inline auto log_softmax (const T &x );

// ./stan/math/rev/fun/log_sum_exp.hpp:45

inline var log_sum_exp (const var &a , const var &b );

// ./stan/math/rev/fun/log_sum_exp.hpp:51

inline var log_sum_exp (const var &a , double b );

// ./stan/math/rev/fun/log_sum_exp.hpp:57

inline var log_sum_exp (double a , const var &b );

// ./stan/math/rev/fun/log_sum_exp.hpp:67
template <typename T,
    require_eigen_st<is_var,
    T>*>
inline var log_sum_exp (const T &v );

// ./stan/math/rev/fun/log_sum_exp.hpp:88
template <typename T,
    require_var_matrix_t<T>*>
inline var log_sum_exp (const T &x );

// ./stan/math/rev/fun/log_sum_exp.hpp:101
template <typename T,
    require_std_vector_st<is_var,
    T>*>
inline auto log_sum_exp (const T &x );

// ./stan/math/rev/fun/logit.hpp:19
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto logit (const var_value <T >&u );

// ./stan/math/rev/fun/matrix_exp_multiply.hpp:24
template <typename Ta,
    typename Tb,
    require_all_eigen_t<Ta,
    Tb>*>
inline Eigen ::Matrix <return_type_t <Ta ,Tb >,-1,Tb ::ColsAtCompileTime >matrix_exp_multiply (const Ta &A , const Tb &B );

// ./stan/math/rev/fun/matrix_power.hpp:28
template <typename T,
    require_rev_matrix_t<T>*>
inline plain_type_t <T >matrix_power (const T &M , const int n );

// ./stan/math/rev/fun/mdivide_left.hpp:28
template <typename T1,
    typename T2,
    require_all_matrix_t<T1,
    T2>*>
inline auto mdivide_left (const T1 &A , const T2 &B );

// ./stan/math/rev/fun/mdivide_left_ldlt.hpp:25
template <typename T1,
    typename T2,
    require_all_matrix_t<T1,
    T2>*>
inline auto mdivide_left_ldlt (LDLT_factor <T1 >&A , const T2 &B );

// ./stan/math/rev/fun/mdivide_left_spd.hpp:147
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_matrix_base_vt<is_var,
    EigMat1,
    EigMat2>*>
inline Eigen ::Matrix <var ,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_left_spd (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/rev/fun/mdivide_left_spd.hpp:178
template <typename EigMat1,
    typename EigMat2,
    require_eigen_matrix_base_vt<is_var,
    EigMat1>*>
inline Eigen ::Matrix <var ,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_left_spd (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/rev/fun/mdivide_left_spd.hpp:209
template <typename EigMat1,
    typename EigMat2,
    require_eigen_matrix_base_vt<std::is_arithmetic,
    EigMat1>*>
inline Eigen ::Matrix <var ,EigMat1 ::RowsAtCompileTime ,EigMat2 ::ColsAtCompileTime >mdivide_left_spd (const EigMat1 &A , const EigMat2 &b );

// ./stan/math/rev/fun/mdivide_left_spd.hpp:259
template <typename T1,
    typename T2,
    require_all_matrix_t<T1,
    T2>*>
inline auto mdivide_left_spd (const T1 &A , const T2 &B );

// ./stan/math/rev/fun/mdivide_left_tri.hpp:245
template <Eigen::UpLoTypeTriView,
    typename T1,
    typename T2,
    require_all_eigen_vt<is_var,
    T1,
    T2>*>
inline Eigen ::Matrix <var ,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left_tri (const T1 &A , const T2 &b );

// ./stan/math/rev/fun/mdivide_left_tri.hpp:270
template <Eigen::UpLoTypeTriView,
    typename T1,
    typename T2,
    require_eigen_vt<std::is_arithmetic,
    T1>*>
inline Eigen ::Matrix <var ,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left_tri (const T1 &A , const T2 &b );

// ./stan/math/rev/fun/mdivide_left_tri.hpp:296
template <Eigen::UpLoTypeTriView,
    typename T1,
    typename T2,
    require_eigen_vt<is_var,
    T1>*>
inline Eigen ::Matrix <var ,T1 ::RowsAtCompileTime ,T2 ::ColsAtCompileTime >mdivide_left_tri (const T1 &A , const T2 &b );

// ./stan/math/rev/fun/mdivide_left_tri.hpp:342
template <Eigen::UpLoTypeTriView,
    typename T1,
    typename T2,
    require_all_matrix_t<T1,
    T2>*>
inline auto mdivide_left_tri (const T1 &A , const T2 &B );

// ./stan/math/rev/fun/modified_bessel_first_kind.hpp:26

inline var modified_bessel_first_kind (int v , const var &a );

// ./stan/math/rev/fun/modified_bessel_second_kind.hpp:26

inline var modified_bessel_second_kind (int v , const var &a );

// ./stan/math/rev/fun/multiply_log.hpp:56

inline var multiply_log (const var &a , const var &b );

// ./stan/math/rev/fun/multiply_log.hpp:69

inline var multiply_log (const var &a , double b );

// ./stan/math/rev/fun/multiply_log.hpp:82

inline var multiply_log (double a , const var &b );

// ./stan/math/rev/fun/multiply_log.hpp:101
template <typename T1,
    typename T2,
    require_all_matrix_t<T1,
    T2>*>
inline auto multiply_log (const T1 &a , const T2 &b );

// ./stan/math/rev/fun/multiply_log.hpp:149
template <typename T1,
    typename T2,
    require_var_matrix_t<T1>*>
inline auto multiply_log (const T1 &a , const T2 &b );

// ./stan/math/rev/fun/multiply_log.hpp:195
template <typename T1,
    typename T2,
    require_stan_scalar_t<T1>*>
inline auto multiply_log (const T1 &a , const T2 &b );

// ./stan/math/rev/fun/norm1.hpp:22
template <typename T,
    require_eigen_vector_vt<is_var,
    T>*>
inline var norm1 (const T &v );

// ./stan/math/rev/fun/norm1.hpp:40
template <typename T,
    require_var_matrix_t<T>*>
inline var norm1 (const T &v );

// ./stan/math/rev/fun/norm2.hpp:21
template <typename T,
    require_eigen_vector_vt<is_var,
    T>*>
inline var norm2 (const T &v );

// ./stan/math/rev/fun/norm2.hpp:38
template <typename T,
    require_var_matrix_t<T>*>
inline var norm2 (const T &v );

// ./stan/math/rev/fun/owens_t.hpp:28
template <typename Var1,
    typename Var2,
    require_all_st_var<Var1,
    Var2>*>
inline auto owens_t (const Var1 &h , const Var2 &a );

// ./stan/math/rev/fun/owens_t.hpp:65
template <typename Var,
    typename Arith,
    require_st_arithmetic<Arith>*>
inline auto owens_t (const Var &h , const Arith &a );

// ./stan/math/rev/fun/owens_t.hpp:97
template <typename Arith,
    typename Var,
    require_st_arithmetic<Arith>*>
inline auto owens_t (const Arith &h , const Var &a );

// ./stan/math/rev/fun/primitive_value.hpp:18

inline double primitive_value (const var &v );

// ./stan/math/rev/fun/proj.hpp:20

inline std ::complex <var >proj (const std ::complex <var >&z );

// ./stan/math/rev/fun/quad_form.hpp:240
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline promote_scalar_t <var ,EigMat2 >quad_form (const EigMat1 &A , const EigMat2 &B , bool symmetric =false );

// ./stan/math/rev/fun/quad_form.hpp:271
template <typename EigMat,
    typename ColVec,
    require_eigen_t<EigMat>*>
inline var quad_form (const EigMat &A , const ColVec &B , bool symmetric =false );

// ./stan/math/rev/fun/quad_form.hpp:307
template <typename Mat1,
    typename Mat2,
    require_all_matrix_t<Mat1,
    Mat2>*>
inline auto quad_form (const Mat1 &A , const Mat2 &B , bool symmetric =false );

// ./stan/math/rev/fun/quad_form.hpp:333
template <typename Mat,
    typename Vec,
    require_matrix_t<Mat>*>
inline var quad_form (const Mat &A , const Vec &B , bool symmetric =false );

// ./stan/math/rev/fun/quad_form_sym.hpp:30
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline auto quad_form_sym (const EigMat1 &A , const EigMat2 &B );

// ./stan/math/rev/fun/read_corr_L.hpp:35
template <typename T,
    require_var_vector_t<T>*>
auto read_corr_L (const T &CPCs , size_t K );

// ./stan/math/rev/fun/read_corr_L.hpp:133
template <typename T1,
    typename T2,
    require_var_vector_t<T1>*>
auto read_corr_L (const T1 &CPCs , size_t K , T2 &log_prob );

// ./stan/math/rev/fun/read_corr_matrix.hpp:26
template <typename T_CPCs,
    require_var_vector_t<T_CPCs>*>
inline var_value <Eigen ::MatrixXd >read_corr_matrix (const T_CPCs &CPCs , size_t K );

// ./stan/math/rev/fun/read_corr_matrix.hpp:55
template <typename T_CPCs,
    require_var_vector_t<T_CPCs>*>
inline var_value <Eigen ::MatrixXd >read_corr_matrix (const T_CPCs &CPCs , size_t K , scalar_type_t <T_CPCs >&log_prob );

// ./stan/math/rev/fun/read_cov_L.hpp:29
template <typename T_CPCs,
    typename T_sds,
    require_any_var_vector_t<T_CPCs,
    T_sds>*>
inline auto read_cov_L (const T_CPCs &CPCs , const T_sds &sds , scalar_type_t <T_CPCs >&log_prob );

// ./stan/math/rev/fun/rows_dot_product.hpp:30
template <typename Mat1,
    typename Mat2,
    require_all_eigen_t<Mat1,
    Mat2>*>
inline Eigen ::Matrix <var ,Mat1 ::RowsAtCompileTime ,1>rows_dot_product (const Mat1 &v1 , const Mat2 &v2 );

// ./stan/math/rev/fun/rows_dot_product.hpp:60
template <typename Mat1,
    typename Mat2,
    require_all_matrix_t<Mat1,
    Mat2>*>
inline auto rows_dot_product (const Mat1 &v1 , const Mat2 &v2 );

// ./stan/math/rev/fun/read_cov_matrix.hpp:27
template <typename T_CPCs,
    typename T_sds,
    require_all_var_vector_t<T_CPCs,
    T_sds>*>
var_value <Eigen ::MatrixXd >read_cov_matrix (const T_CPCs &CPCs , const T_sds &sds , scalar_type_t <T_CPCs >&log_prob );

// ./stan/math/rev/fun/read_cov_matrix.hpp:44
template <typename T_CPCs,
    typename T_sds,
    require_all_var_vector_t<T_CPCs,
    T_sds>*>
inline var_value <Eigen ::MatrixXd >read_cov_matrix (const T_CPCs &CPCs , const T_sds &sds );

// ./stan/math/rev/fun/rep_matrix.hpp:23
template <typename Ret,
    typename T,
    require_var_matrix_t<Ret>*>
inline auto rep_matrix (const T &x , int m , int n );

// ./stan/math/rev/fun/rep_matrix.hpp:43
template <typename Ret,
    typename Vec,
    require_var_matrix_t<Ret>*>
inline auto rep_matrix (const Vec &x , int n );

// ./stan/math/rev/fun/rep_row_vector.hpp:20
template <typename T_ret,
    require_var_matrix_t<T_ret>*
;

// ./stan/math/rev/fun/rep_vector.hpp:20
template <typename T_ret,
    require_var_matrix_t<T_ret>*
;

// ./stan/math/rev/fun/rising_factorial.hpp:25

inline var rising_factorial (const var &a , int b );

// ./stan/math/rev/fun/round.hpp:43

inline var round (const var &a );

// ./stan/math/rev/fun/rows_dot_self.hpp:19
template <typename Mat,
    require_eigen_vt<is_var,
    Mat>*>
inline Eigen ::Matrix <var ,Mat ::RowsAtCompileTime ,1>rows_dot_self (const Mat &x );

// ./stan/math/rev/fun/rows_dot_self.hpp:35
template <typename Mat,
    require_var_matrix_t<Mat>*>
inline auto rows_dot_self (const Mat &x );

// ./stan/math/rev/fun/sd.hpp:25
template <typename T,
    require_eigen_st<is_var,
    T>*>
var sd (const T &x );

// ./stan/math/rev/fun/sd.hpp:68
template <typename T,
    require_var_matrix_t<T>*>
var sd (const T &x );

// ./stan/math/rev/fun/sd.hpp:94
template <typename T,
    require_std_vector_st<is_var,
    T>*>
auto sd (const T &m );

// ./stan/math/rev/fun/singular_values.hpp:22
template <typename EigMat,
    require_rev_matrix_t<EigMat>*>
inline auto singular_values (const EigMat &m );

// ./stan/math/rev/fun/svd.hpp:26
template <typename EigMat,
    require_rev_matrix_t<EigMat>*>
inline auto svd (const EigMat &m );

// ./stan/math/rev/fun/svd_U.hpp:24
template <typename EigMat,
    require_rev_matrix_t<EigMat>*>
inline auto svd_U (const EigMat &m );

// ./stan/math/rev/fun/svd_V.hpp:24
template <typename EigMat,
    require_rev_matrix_t<EigMat>*>
inline auto svd_V (const EigMat &m );

// ./stan/math/rev/fun/softmax.hpp:27
template <typename Mat,
    require_rev_matrix_t<Mat>*>
inline auto softmax (const Mat &alpha );

// ./stan/math/rev/fun/squared_distance.hpp:21

inline var squared_distance (const var &a , const var &b );

// ./stan/math/rev/fun/squared_distance.hpp:35

inline var squared_distance (const var &a , double b );

// ./stan/math/rev/fun/squared_distance.hpp:47

inline var squared_distance (double a , const var &b );

// ./stan/math/rev/fun/squared_distance.hpp:116
template <typename EigVecVar1,
    typename EigVecVar2,
    require_all_eigen_vector_vt<is_var,
    EigVecVar1,
    EigVecVar2>*>
inline var squared_distance (const EigVecVar1 &v1 , const EigVecVar2 &v2 );

// ./stan/math/rev/fun/squared_distance.hpp:124
template <typename EigVecVar,
    typename EigVecArith,
    require_eigen_vector_vt<is_var,
    EigVecVar>*>
inline var squared_distance (const EigVecVar &v1 , const EigVecArith &v2 );

// ./stan/math/rev/fun/squared_distance.hpp:132
template <typename EigVecArith,
    typename EigVecVar,
    require_eigen_vector_vt<std::is_arithmetic,
    EigVecArith>*>
inline var squared_distance (const EigVecArith &v1 , const EigVecVar &v2 );

// ./stan/math/rev/fun/squared_distance.hpp:155
template <typename T1,
    typename T2,
    require_all_vector_t<T1,
    T2>*>
inline var squared_distance (const T1 &A , const T2 &B );

// ./stan/math/rev/fun/stan_print.hpp:11

inline void stan_print (std ::ostream *o , const var &x );

// ./stan/math/rev/fun/step.hpp:26

inline var step (const var &a );

// ./stan/math/rev/fun/tanh.hpp:42

inline var tanh (const var &a );

// ./stan/math/rev/fun/tanh.hpp:56
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto tanh (const VarMat &a );

// ./stan/math/rev/fun/tanh.hpp:70

inline std ::complex <var >tanh (const std ::complex <var >&z );

// ./stan/math/rev/fun/tan.hpp:47

inline var tan (const var &a );

// ./stan/math/rev/fun/tan.hpp:61
template <typename VarMat,
    require_var_matrix_t<VarMat>*>
inline auto tan (const VarMat &a );

// ./stan/math/rev/fun/tan.hpp:76

inline std ::complex <var >tan (const std ::complex <var >&z );

// ./stan/math/rev/fun/tcrossprod.hpp:23
template <typename T,
    require_rev_matrix_t<T>*>
inline auto tcrossprod (const T &M );

// ./stan/math/rev/fun/tgamma.hpp:48

inline var tgamma (const var &a );

// ./stan/math/rev/fun/tgamma.hpp:61
template <typename T,
    require_var_matrix_t<T>*>
inline auto tgamma (const T &a );

// ./stan/math/rev/fun/to_var.hpp:22

inline var to_var (double x );

// ./stan/math/rev/fun/to_var.hpp:30

inline var &to_var (var &x );

// ./stan/math/rev/fun/to_var.hpp:38

inline const var &to_var (const var &x );

// ./stan/math/rev/fun/to_var.hpp:48

inline std ::vector <var >to_var (const std ::vector <double >&v );

// ./stan/math/rev/fun/to_var.hpp:64

inline const std ::vector <var >&to_var (const std ::vector <var >&v );

// ./stan/math/rev/fun/to_var.hpp:84

inline matrix_v to_var (const matrix_d &m );

// ./stan/math/rev/fun/to_var.hpp:95

inline matrix_v &to_var (matrix_v &m );

// ./stan/math/rev/fun/to_var.hpp:103

inline const matrix_v &to_var (const matrix_v &m );

// ./stan/math/rev/fun/to_var.hpp:114

inline vector_v to_var (const vector_d &v );

// ./stan/math/rev/fun/to_var.hpp:125

inline const vector_v &to_var (const vector_v &v );

// ./stan/math/rev/fun/to_var.hpp:133

inline vector_v &to_var (vector_v &v );

// ./stan/math/rev/fun/to_var.hpp:144

inline row_vector_v to_var (const row_vector_d &rv );

// ./stan/math/rev/fun/to_var.hpp:155

inline const row_vector_v &to_var (const row_vector_v &rv );

// ./stan/math/rev/fun/to_var.hpp:163

inline row_vector_v &to_var (row_vector_v &rv );

// ./stan/math/rev/fun/to_vector.hpp:19
template <typename EigMat,
    require_eigen_t<EigMat>*>
inline auto to_vector (const var_value <EigMat >&x );

// ./stan/math/rev/fun/trace.hpp:22
template <typename T,
    require_rev_matrix_t<T>*>
inline auto trace (const T &m );

// ./stan/math/rev/fun/trace_inv_quad_form_ldlt.hpp:30
template <typename T1,
    typename T2,
    require_all_matrix_t<T1,
    T2>*>
inline var trace_inv_quad_form_ldlt (LDLT_factor <T1 >&A , const T2 &B );

// ./stan/math/rev/fun/trace_gen_inv_quad_form_ldlt.hpp:29
template <typename Td,
    typename Ta,
    typename Tb,
    require_not_col_vector_t<Td>*>
inline var trace_gen_inv_quad_form_ldlt (const Td &D , LDLT_factor <Ta >&A , const Tb &B );

// ./stan/math/rev/fun/trace_gen_inv_quad_form_ldlt.hpp:186
template <typename Td,
    typename Ta,
    typename Tb,
    require_col_vector_t<Td>*>
inline var trace_gen_inv_quad_form_ldlt (const Td &D , const LDLT_factor <Ta >&A , const Tb &B );

// ./stan/math/rev/fun/trace_gen_quad_form.hpp:136
template <typename Td,
    typename Ta,
    typename Tb,
    require_any_var_matrix_t<Td,
    Ta,
    Tb>*>
inline var trace_gen_quad_form (const Td &D , const Ta &A , const Tb &B );

// ./stan/math/rev/fun/trace_quad_form.hpp:77
template <typename EigMat1,
    typename EigMat2,
    require_all_eigen_t<EigMat1,
    EigMat2>*>
inline return_type_t <EigMat1 ,EigMat2 >trace_quad_form (const EigMat1 &A , const EigMat2 &B );

// ./stan/math/rev/fun/trace_quad_form.hpp:115
template <typename Mat1,
    typename Mat2,
    require_all_matrix_t<Mat1,
    Mat2>*>
inline var trace_quad_form (const Mat1 &A , const Mat2 &B );

// ./stan/math/rev/fun/trigamma.hpp:23

inline var trigamma (const var &u );

// ./stan/math/rev/fun/trunc.hpp:43

inline var trunc (const var &a );

// ./stan/math/rev/fun/variance.hpp:42

inline var variance (const std ::vector <var >&v );

// ./stan/math/rev/fun/variance.hpp:60
template <typename EigMat,
    require_eigen_vt<is_var,
    EigMat>*>
var variance (const EigMat &m );

// ./stan/math/rev/fun/variance.hpp:78
template <typename Mat,
    require_var_matrix_t<Mat>*>
inline var variance (const Mat &x );

// ./stan/math/prim/functor/finite_diff_gradient.hpp:42
template <typename F>
void finite_diff_gradient (const F &f , const Eigen ::VectorXd &x , double &fx , Eigen ::VectorXd &grad_fx , double epsilon =1e-03);

// ./stan/math/prim/functor/hcubature.hpp:403
template <typename F,
    typename T_a,
    typename T_b,
    typename ParsTuple,
    typename TAbsErr,
    typename TRelErr>
inline auto hcubature (const F &integrand , const ParsTuple &pars , const int dim , const Eigen ::Matrix <T_a , Eigen ::Dynamic , 1>&a , const Eigen ::Matrix <T_b , Eigen ::Dynamic , 1>&b , const int max_eval , const TAbsErr reqAbsError , const TRelErr reqRelError );

// ./stan/math/prim/functor/ode_store_sensitivities.hpp:31
template <typename F,
    typename T_y0_t0,
    typename T_t0,
    typename T_t,
    typename ...Args,
    typename >
Eigen ::VectorXd ode_store_sensitivities (const F &f , const std ::vector <double >&coupled_state , const Eigen ::Matrix <T_y0_t0 , Eigen ::Dynamic , 1>&y0 , T_t0 t0 , T_t t , std ::ostream *msgs , const Args &...args );

// ./stan/math/prim/functor/ode_rk45.hpp:55
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...Args,
    require_eigen_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,Args ...>,Eigen ::Dynamic ,1>>ode_rk45_tol_impl (const char *function_name , const F &f , const T_y0 &y0_arg , T_t0 t0 , const std ::vector <T_ts >&ts , double relative_tolerance , double absolute_tolerance , long int max_num_steps , // NOLINT(runtime/int)std ::ostream *msgs , const Args &...args );

// ./stan/math/prim/functor/ode_rk45.hpp:197
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...Args,
    require_eigen_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,Args ...>,Eigen ::Dynamic ,1>>ode_rk45_tol (const F &f , const T_y0 &y0_arg , T_t0 t0 , const std ::vector <T_ts >&ts , double relative_tolerance , double absolute_tolerance , long int max_num_steps , // NOLINT(runtime/int)std ::ostream *msgs , const Args &...args );

// ./stan/math/prim/functor/ode_rk45.hpp:243
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...Args,
    require_eigen_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,Args ...>,Eigen ::Dynamic ,1>>ode_rk45 (const F &f , const T_y0 &y0 , T_t0 t0 , const std ::vector <T_ts >&ts , std ::ostream *msgs , const Args &...args );

// ./stan/math/prim/functor/integrate_ode_rk45.hpp:16
template <typename F,
    typename T_y0,
    typename T_param,
    typename T_t0,
    typename T_ts>
inline auto integrate_ode_rk45 (const F &f , const std ::vector <T_y0 >&y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , const std ::vector <T_param >&theta , const std ::vector <double >&x , const std ::vector <int >&x_int , std ::ostream *msgs =nullptr , double relative_tolerance =1e-6, double absolute_tolerance =1e-6, int max_num_steps =1e6);

// ./stan/math/prim/functor/ode_ckrk.hpp:54
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...Args,
    require_eigen_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,Args ...>,Eigen ::Dynamic ,1>>ode_ckrk_tol_impl (const char *function_name , const F &f , const T_y0 &y0_arg , T_t0 t0 , const std ::vector <T_ts >&ts , double relative_tolerance , double absolute_tolerance , long int max_num_steps , // NOLINT(runtime/int)std ::ostream *msgs , const Args &...args );

// ./stan/math/prim/functor/ode_ckrk.hpp:195
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...Args,
    require_eigen_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,Args ...>,Eigen ::Dynamic ,1>>ode_ckrk_tol (const F &f , const T_y0 &y0_arg , T_t0 t0 , const std ::vector <T_ts >&ts , double relative_tolerance , double absolute_tolerance , long int max_num_steps , // NOLINT(runtime/int)std ::ostream *msgs , const Args &...args );

// ./stan/math/prim/functor/ode_ckrk.hpp:241
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...Args,
    require_eigen_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,Args ...>,Eigen ::Dynamic ,1>>ode_ckrk (const F &f , const T_y0 &y0 , T_t0 t0 , const std ::vector <T_ts >&ts , std ::ostream *msgs , const Args &...args );

// ./stan/math/prim/functor/map_rect.hpp:123
template <int call_id,
    typename F,
    typename T_shared_param,
    typename T_job_param,
    require_eigen_col_vector_t<T_shared_param>*>
Eigen ::Matrix <return_type_t <T_shared_param ,T_job_param >,Eigen ::Dynamic ,1>map_rect (const T_shared_param &shared_params , const std ::vector <Eigen ::Matrix <T_job_param , Eigen ::Dynamic , 1>>&job_params , const std ::vector <std ::vector <double >>&x_r , const std ::vector <std ::vector <int >>&x_i , std ::ostream *msgs =nullptr );

// ./stan/math/prim/functor/reduce_sum.hpp:198
template <typename ReduceFunction,
    typename Vec,
    typename >
inline auto reduce_sum (Vec &&vmapped , int grainsize , std ::ostream *msgs , Args &&...args );

// ./stan/math/prim/functor/reduce_sum_static.hpp:43
template <typename ReduceFunction,
    typename Vec,
    typename >
auto reduce_sum_static (Vec &&vmapped , int grainsize , std ::ostream *msgs , Args &&...args );

// ./stan/math/rev/functor/jacobian.hpp:13
template <typename F>
void jacobian (const F &f , const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&x , Eigen ::Matrix <double , Eigen ::Dynamic , 1>&fx , Eigen ::Matrix <double , Eigen ::Dynamic , Eigen ::Dynamic >&J );

// ./stan/math/rev/functor/algebra_system.hpp:100
template <typename T1,
    typename T2>
void algebra_solver_check (const Eigen ::Matrix <T1 , Eigen ::Dynamic , 1>&x , const Eigen ::Matrix <T2 , Eigen ::Dynamic , 1>y , const std ::vector <double >&dat , const std ::vector <int >&dat_int , double function_tolerance , long int max_num_steps );

// ./stan/math/rev/functor/algebra_solver_fp.hpp:364
template <typename F,
    typename T1,
    typename T2,
    typename T_u,
    typename T_f>
Eigen ::Matrix <T2 ,-1,1>algebra_solver_fp (const F &f , const Eigen ::Matrix <T1 , -1, 1>&x , const Eigen ::Matrix <T2 , -1, 1>&y , const std ::vector <double >&x_r , const std ::vector <int >&x_i , const std ::vector <T_u >&u_scale , const std ::vector <T_f >&f_scale , std ::ostream *msgs =nullptr , double f_tol =1e-8, int max_num_steps =200);

// ./stan/math/rev/functor/solve_powell.hpp:46
template <typename F,
    typename T,
    typename ...Args,
    require_eigen_vector_t<T>*>
T &solve_powell_call_solver (const F &f , T &x , std ::ostream *const msgs , const double relative_tolerance , const double function_tolerance , const int64_t max_num_steps , const Args &...args );

// ./stan/math/rev/functor/solve_powell.hpp:125
template <typename F,
    typename T,
    typename ...Args,
    require_eigen_vector_t<T>*>
Eigen ::VectorXd solve_powell_tol (const F &f , const T &x , const double relative_tolerance , const double function_tolerance , const int64_t max_num_steps , std ::ostream *const msgs , const Args &...args );

// ./stan/math/rev/functor/solve_powell.hpp:186
template <typename F,
    typename T,
    typename ...T_Args,
    require_eigen_vector_t<T>*>
Eigen ::Matrix <stan ::return_type_t <T_Args ...>,Eigen ::Dynamic ,1>solve_powell (const F &f , const T &x , std ::ostream *const msgs , const T_Args &...args );

// ./stan/math/rev/functor/solve_powell.hpp:249
template <typename F,
    typename T1,
    typename T2,
    require_all_eigen_vector_t<T1,
    T2>*>
Eigen ::Matrix <value_type_t <T2 >,Eigen ::Dynamic ,1>algebra_solver (const F &f , const T1 &x , const T2 &y , const std ::vector <double >&dat , const std ::vector <int >&dat_int , std ::ostream *msgs =nullptr , const double relative_tolerance =1e-10, const double function_tolerance =1e-6, const int64_t max_num_steps =1e+3);

// ./stan/math/rev/functor/solve_powell.hpp:322
template <typename F,
    typename T,
    typename ...T_Args,
    require_eigen_vector_t<T>*>
Eigen ::Matrix <var ,Eigen ::Dynamic ,1>solve_powell_tol (const F &f , const T &x , const double relative_tolerance , const double function_tolerance , const int64_t max_num_steps , std ::ostream *const msgs , const T_Args &...args );

// ./stan/math/rev/functor/kinsol_solve.hpp:60
template <typename F1,
    typename ...Args>
Eigen ::VectorXd kinsol_solve (const F1 &f , const Eigen ::VectorXd &x , const double scaling_step_tol , // = 1e-3const double function_tolerance , // = 1e-6const int64_t max_num_steps , // = 200const bool custom_jacobian , // = 1const int steps_eval_jacobian , // = 10const int global_line_search , // = KIN_LINESEARCHstd ::ostream *const msgs , const Args &...args );

// ./stan/math/rev/functor/solve_newton.hpp:56
template <typename F,
    typename T,
    typename ...Args,
    require_eigen_vector_t<T>*>
Eigen ::VectorXd solve_newton_tol (const F &f , const T &x , const double scaling_step_size , const double function_tolerance , const int64_t max_num_steps , std ::ostream *const msgs , const Args &...args );

// ./stan/math/rev/functor/solve_newton.hpp:135
template <typename F,
    typename T,
    typename ...T_Args,
    require_eigen_vector_t<T>*>
Eigen ::Matrix <var ,Eigen ::Dynamic ,1>solve_newton_tol (const F &f , const T &x , const double scaling_step_size , const double function_tolerance , const int64_t max_num_steps , std ::ostream *const msgs , const T_Args &...args );

// ./stan/math/rev/functor/solve_newton.hpp:237
template <typename F,
    typename T,
    typename ...T_Args,
    require_eigen_vector_t<T>*>
Eigen ::Matrix <stan ::return_type_t <T_Args ...>,Eigen ::Dynamic ,1>solve_newton (const F &f , const T &x , std ::ostream *const msgs , const T_Args &...args );

// ./stan/math/rev/functor/solve_newton.hpp:293
template <typename F,
    typename T1,
    typename T2,
    require_all_eigen_vector_t<T1,
    T2>*>
Eigen ::Matrix <scalar_type_t <T2 >,Eigen ::Dynamic ,1>algebra_solver_newton (const F &f , const T1 &x , const T2 &y , const std ::vector <double >&dat , const std ::vector <int >&dat_int , std ::ostream *const msgs =nullptr , const double scaling_step_size =1e-3, const double function_tolerance =1e-6, const long int max_num_steps =200);

// ./stan/math/rev/functor/apply_scalar_binary.hpp:28
template <typename T1,
    typename T2,
    typename F,
    require_any_var_matrix_t<T1,
    T2>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/rev/functor/apply_scalar_binary.hpp:48
template <typename T1,
    typename T2,
    typename F,
    require_any_var_matrix_t<T1,
    T2>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/rev/functor/apply_scalar_binary.hpp:70
template <typename T1,
    typename T2,
    typename F,
    require_any_std_vector_vt<is_std_vector,
    T1,
    T2>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/rev/functor/apply_scalar_binary.hpp:93
template <typename T1,
    typename T2,
    typename F,
    require_any_stan_scalar_t<T1,
    T2>*>
inline auto apply_scalar_binary (const T1 &x , const T2 &y , const F &f );

// ./stan/math/rev/functor/cvodes_utils.hpp:13

inline void cvodes_set_options (void *cvodes_mem , // NOLINTNEXTLINE(runtime/int)long int max_num_steps );

// ./stan/math/rev/functor/ode_store_sensitivities.hpp:32
template <typename F,
    typename T_y0_t0,
    typename T_t0,
    typename T_t,
    typename ...Args,
    require_any_autodiff_t<T_y0_t0,
    T_t0,
    T_t,
    scalar_type_t<Args>...>*>
Eigen ::Matrix <var ,Eigen ::Dynamic ,1>ode_store_sensitivities (const F &f , const std ::vector <double >&coupled_state , const Eigen ::Matrix <T_y0_t0 , Eigen ::Dynamic , 1>&y0 , const T_t0 &t0 , const T_t &t , std ::ostream *msgs , const Args &...args );

// ./stan/math/rev/functor/gradient.hpp:45
template <typename F>
void gradient (const F &f , const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&x , double &fx , Eigen ::Matrix <double , Eigen ::Dynamic , 1>&grad_fx );

// ./stan/math/rev/functor/gradient.hpp:100
template <typename F,
    typename EigVec,
    typename InputIt,
    require_eigen_vector_vt<std::is_arithmetic,
    EigVec>*>
void gradient (const F &f , const EigVec &x , double &fx , InputIt first_grad_fx , InputIt last_grad_fx );

// ./stan/math/rev/functor/integrate_1d.hpp:39
template <typename F,
    typename T_a,
    typename T_b,
    typename ...Args,
    require_any_st_var<T_a,
    T_b,
    Args...>*>
inline return_type_t <T_a ,T_b ,Args ...>integrate_1d_impl (const F &f , const T_a &a , const T_b &b , double relative_tolerance , std ::ostream *msgs , const Args &...args );

// ./stan/math/rev/functor/integrate_1d.hpp:214
template <typename F,
    typename T_a,
    typename T_b,
    typename T_theta,
    typename >
inline return_type_t <T_a ,T_b ,T_theta >integrate_1d (const F &f , const T_a &a , const T_b &b , const std ::vector <T_theta >&theta , const std ::vector <double >&x_r , const std ::vector <int >&x_i , std ::ostream *msgs , const double relative_tolerance =std ::sqrt (EPSILON ));

// ./stan/math/rev/functor/dae.hpp:50
template <typename F,
    typename T_yy,
    typename T_yp,
    typename ...T_Args,
    require_all_eigen_col_vector_t<T_yy,
    T_yp>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_yy ,T_yp ,T_Args ...>,-1,1>>dae_tol_impl (const char *func , const F &f , const T_yy &yy0 , const T_yp &yp0 , double t0 , const std ::vector <double >&ts , double rtol , double atol , int64_t max_num_steps , std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/dae.hpp:122
template <typename F,
    typename T_yy,
    typename T_yp,
    typename ...T_Args,
    require_all_eigen_col_vector_t<T_yy,
    T_yp>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_yy ,T_yp ,T_Args ...>,-1,1>>dae_tol (const F &f , const T_yy &yy0 , const T_yp &yp0 , double t0 , const std ::vector <double >&ts , double rtol , double atol , int64_t max_num_steps , std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/dae.hpp:164
template <typename F,
    typename T_yy,
    typename T_yp,
    typename ...T_Args,
    require_all_eigen_col_vector_t<T_yy,
    T_yp>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_yy ,T_yp ,T_Args ...>,-1,1>>dae (const F &f , const T_yy &yy0 , const T_yp &yp0 , double t0 , const std ::vector <double >&ts , std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/ode_adams.hpp:48
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...T_Args,
    require_eigen_col_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,T_Args ...>,Eigen ::Dynamic ,1>>ode_adams_tol_impl (const char *function_name , const F &f , const T_y0 &y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , double relative_tolerance , double absolute_tolerance , long int max_num_steps , // NOLINT(runtime/int)std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/ode_adams.hpp:102
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...T_Args,
    require_eigen_col_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,T_Args ...>,Eigen ::Dynamic ,1>>ode_adams_tol (const F &f , const T_y0 &y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , double relative_tolerance , double absolute_tolerance , long int max_num_steps , // NOLINT(runtime/int)std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/ode_adams.hpp:145
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...T_Args,
    require_eigen_col_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,T_Args ...>,Eigen ::Dynamic ,1>>ode_adams (const F &f , const T_y0 &y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/ode_bdf.hpp:49
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...T_Args,
    require_eigen_col_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,T_Args ...>,Eigen ::Dynamic ,1>>ode_bdf_tol_impl (const char *function_name , const F &f , const T_y0 &y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , double relative_tolerance , double absolute_tolerance , long int max_num_steps , // NOLINT(runtime/int)std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/ode_bdf.hpp:103
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...T_Args,
    require_eigen_col_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,T_Args ...>,Eigen ::Dynamic ,1>>ode_bdf_tol (const F &f , const T_y0 &y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , double relative_tolerance , double absolute_tolerance , long int max_num_steps , // NOLINT(runtime/int)std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/ode_bdf.hpp:146
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename ...T_Args,
    require_eigen_col_vector_t<T_y0>*>
std ::vector <Eigen ::Matrix <stan ::return_type_t <T_y0 ,T_t0 ,T_ts ,T_Args ...>,Eigen ::Dynamic ,1>>ode_bdf (const F &f , const T_y0 &y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/ode_adjoint.hpp:66
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename T_abs_tol_fwd,
    typename T_abs_tol_bwd,
    typename ...T_Args,
    require_all_eigen_col_vector_t<T_y0,
    T_abs_tol_fwd,
    T_abs_tol_bwd>*>
auto ode_adjoint_impl (const char *function_name , F &&f , const T_y0 &y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , double relative_tolerance_forward , const T_abs_tol_fwd &absolute_tolerance_forward , double relative_tolerance_backward , const T_abs_tol_bwd &absolute_tolerance_backward , double relative_tolerance_quadrature , double absolute_tolerance_quadrature , long int max_num_steps , // NOLINT(runtime/int)long int num_steps_between_checkpoints , // NOLINT(runtime/int)int interpolation_polynomial , int solver_forward , int solver_backward , std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/ode_adjoint.hpp:150
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename T_abs_tol_fwd,
    typename T_abs_tol_bwd,
    typename ...T_Args,
    require_all_eigen_col_vector_t<T_y0,
    T_abs_tol_fwd,
    T_abs_tol_bwd>*>
std ::vector <Eigen ::VectorXd >ode_adjoint_impl (const char *function_name , F &&f , const T_y0 &y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , double relative_tolerance_forward , const T_abs_tol_fwd &absolute_tolerance_forward , double relative_tolerance_backward , const T_abs_tol_bwd &absolute_tolerance_backward , double relative_tolerance_quadrature , double absolute_tolerance_quadrature , long int max_num_steps , // NOLINT(runtime/int)long int num_steps_between_checkpoints , // NOLINT(runtime/int)int interpolation_polynomial , int solver_forward , int solver_backward , std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/functor/ode_adjoint.hpp:239
template <typename F,
    typename T_y0,
    typename T_t0,
    typename T_ts,
    typename T_abs_tol_fwd,
    typename T_abs_tol_bwd,
    typename ...T_Args,
    require_all_eigen_col_vector_t<T_y0,
    T_abs_tol_fwd,
    T_abs_tol_bwd>*>
auto ode_adjoint_tol_ctl (F &&f , const T_y0 &y0 , const T_t0 &t0 , const std ::vector <T_ts >&ts , double relative_tolerance_forward , const T_abs_tol_fwd &absolute_tolerance_forward , double relative_tolerance_backward , const T_abs_tol_bwd &absolute_tolerance_backward , double relative_tolerance_quadrature , double absolute_tolerance_quadrature , long int max_num_steps , // NOLINT(runtime/int)long int num_steps_between_checkpoints , // NOLINT(runtime/int)int interpolation_polynomial , int solver_forward , int solver_backward , std ::ostream *msgs , const T_Args &...args );

// ./stan/math/rev/prob/std_normal_log_qf.hpp:20
template <typename T,
    require_stan_scalar_or_eigen_t<T>*>
inline auto std_normal_log_qf (const var_value <T >&log_p );

// ./stan/math/prim/prob/bernoulli_lccdf.hpp:29
template <typename T_n,
    typename T_prob,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_prob>*>
return_type_t <T_prob >bernoulli_lccdf (const T_n &n , const T_prob &theta );

// ./stan/math/prim/prob/bernoulli_ccdf_log.hpp:13
template <typename T_n,
    typename T_prob>
return_type_t <T_prob >bernoulli_ccdf_log (const T_n &n , const T_prob &theta );

// ./stan/math/prim/prob/bernoulli_cdf.hpp:27
template <typename T_n,
    typename T_prob,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_prob>*>
return_type_t <T_prob >bernoulli_cdf (const T_n &n , const T_prob &theta );

// ./stan/math/prim/prob/bernoulli_lcdf.hpp:27
template <typename T_n,
    typename T_prob,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_prob>*>
return_type_t <T_prob >bernoulli_lcdf (const T_n &n , const T_prob &theta );

// ./stan/math/prim/prob/bernoulli_cdf_log.hpp:13
template <typename T_n,
    typename T_prob>
return_type_t <T_prob >bernoulli_cdf_log (const T_n &n , const T_prob &theta );

// ./stan/math/prim/prob/bernoulli_logit_glm_lpmf.hpp:49
template <bool propto,
    typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta,
    require_matrix_t<T_x>*>
return_type_t <T_x ,T_alpha ,T_beta >bernoulli_logit_glm_lpmf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta );

// ./stan/math/prim/prob/bernoulli_logit_glm_lpmf.hpp:169
template <typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta>
inline return_type_t <T_x ,T_beta ,T_alpha >bernoulli_logit_glm_lpmf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta );

// ./stan/math/prim/prob/bernoulli_logit_glm_rng.hpp:42
template <typename T_x,
    typename T_alpha,
    typename T_beta,
    class RNG>
inline typename VectorBuilder <true ,int ,T_alpha >::type bernoulli_logit_glm_rng (const T_x &x , const T_alpha &alpha , const T_beta &beta , RNG &rng );

// ./stan/math/prim/prob/bernoulli_logit_lpmf.hpp:33
template <bool propto,
    typename T_n,
    typename T_prob,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_prob>*>
return_type_t <T_prob >bernoulli_logit_lpmf (const T_n &n , const T_prob &theta );

// ./stan/math/prim/prob/bernoulli_logit_lpmf.hpp:91
template <typename T_n,
    typename T_prob>
inline return_type_t <T_prob >bernoulli_logit_lpmf (const T_n &n , const T_prob &theta );

// ./stan/math/prim/prob/bernoulli_logit_rng.hpp:29
template <typename T_t,
    class RNG>
inline typename VectorBuilder <true ,int ,T_t >::type bernoulli_logit_rng (const T_t &t , RNG &rng );

// ./stan/math/prim/prob/bernoulli_lpmf.hpp:31
template <bool propto,
    typename T_n,
    typename T_prob,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_prob>*>
return_type_t <T_prob >bernoulli_lpmf (const T_n &n , const T_prob &theta );

// ./stan/math/prim/prob/bernoulli_lpmf.hpp:114
template <typename T_y,
    typename T_prob>
inline return_type_t <T_prob >bernoulli_lpmf (const T_y &n , const T_prob &theta );

// ./stan/math/prim/prob/bernoulli_rng.hpp:28
template <typename T_theta,
    class RNG>
inline typename VectorBuilder <true ,int ,T_theta >::type bernoulli_rng (const T_theta &theta , RNG &rng );

// ./stan/math/prim/prob/beta_binomial_lccdf.hpp:43
template <typename T_n,
    typename T_N,
    typename T_size1,
    typename T_size2>
return_type_t <T_size1 ,T_size2 >beta_binomial_lccdf (const T_n &n , const T_N &N , const T_size1 &alpha , const T_size2 &beta );

// ./stan/math/prim/prob/beta_binomial_ccdf_log.hpp:13
template <typename T_n,
    typename T_N,
    typename T_size1,
    typename T_size2>
return_type_t <T_size1 ,T_size2 >beta_binomial_ccdf_log (const T_n &n , const T_N &N , const T_size1 &alpha , const T_size2 &beta );

// ./stan/math/prim/prob/beta_binomial_cdf.hpp:42
template <typename T_n,
    typename T_N,
    typename T_size1,
    typename T_size2>
return_type_t <T_size1 ,T_size2 >beta_binomial_cdf (const T_n &n , const T_N &N , const T_size1 &alpha , const T_size2 &beta );

// ./stan/math/prim/prob/beta_binomial_lcdf.hpp:43
template <typename T_n,
    typename T_N,
    typename T_size1,
    typename T_size2>
return_type_t <T_size1 ,T_size2 >beta_binomial_lcdf (const T_n &n , const T_N &N , const T_size1 &alpha , const T_size2 &beta );

// ./stan/math/prim/prob/beta_binomial_cdf_log.hpp:13
template <typename T_n,
    typename T_N,
    typename T_size1,
    typename T_size2>
return_type_t <T_size1 ,T_size2 >beta_binomial_cdf_log (const T_n &n , const T_N &N , const T_size1 &alpha , const T_size2 &beta );

// ./stan/math/prim/prob/beta_binomial_lpmf.hpp:39
template <bool propto,
    typename T_n,
    typename T_N,
    typename T_size1,
    typename T_size2,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_N,
    T_size1,
    T_size2>*>
return_type_t <T_size1 ,T_size2 >beta_binomial_lpmf (const T_n &n , const T_N &N , const T_size1 &alpha , const T_size2 &beta );

// ./stan/math/prim/prob/beta_binomial_lpmf.hpp:168
template <typename T_n,
    typename T_N,
    typename T_size1,
    typename T_size2>
return_type_t <T_size1 ,T_size2 >beta_binomial_lpmf (const T_n &n , const T_N &N , const T_size1 &alpha , const T_size2 &beta );

// ./stan/math/prim/prob/binomial_rng.hpp:31
template <typename T_N,
    typename T_theta,
    class RNG>
inline typename VectorBuilder <true ,int ,T_N ,T_theta >::type binomial_rng (const T_N &N , const T_theta &theta , RNG &rng );

// ./stan/math/prim/prob/beta_rng.hpp:35
template <typename T_shape1,
    typename T_shape2,
    class RNG>
inline typename VectorBuilder <true ,double ,T_shape1 ,T_shape2 >::type beta_rng (const T_shape1 &alpha , const T_shape2 &beta , RNG &rng );

// ./stan/math/prim/prob/beta_binomial_rng.hpp:31
template <typename T_N,
    typename T_shape1,
    typename T_shape2,
    class RNG>
inline typename VectorBuilder <true ,int ,T_N ,T_shape1 ,T_shape2 >::type beta_binomial_rng (const T_N &N , const T_shape1 &alpha , const T_shape2 &beta , RNG &rng );

// ./stan/math/prim/prob/beta_lccdf.hpp:40
template <typename T_y,
    typename T_scale_succ,
    typename T_scale_fail>
return_type_t <T_y ,T_scale_succ ,T_scale_fail >beta_lccdf (const T_y &y , const T_scale_succ &alpha , const T_scale_fail &beta_param );

// ./stan/math/prim/prob/beta_ccdf_log.hpp:13
template <typename T_y,
    typename T_scale_succ,
    typename T_scale_fail>
return_type_t <T_y ,T_scale_succ ,T_scale_fail >beta_ccdf_log (const T_y &y , const T_scale_succ &alpha , const T_scale_fail &beta );

// ./stan/math/prim/fun/inc_beta_ddb.hpp:16
template <typename T>
T inc_beta_dda (T a , T b , T z , T digamma_a , T digamma_ab );

// ./stan/math/prim/fun/inc_beta_ddb.hpp:41
template <typename T>
T inc_beta_ddb (T a , T b , T z , T digamma_b , T digamma_ab );

// ./stan/math/prim/fun/inc_beta_dda.hpp:38
template <typename T>
T inc_beta_dda (T a , T b , T z , T digamma_a , T digamma_ab );

// ./stan/math/prim/fun/inc_beta_ddz.hpp:29
template <typename T>
T inc_beta_ddz (T a , T b , T z );

// ./stan/math/prim/fun/inc_beta_ddz.hpp:37
template <>
inline double inc_beta_ddz (double a , double b , double z );

// ./stan/math/prim/prob/beta_cdf.hpp:35
template <typename T_y,
    typename T_scale_succ,
    typename T_scale_fail>
return_type_t <T_y ,T_scale_succ ,T_scale_fail >beta_cdf (const T_y &y , const T_scale_succ &alpha , const T_scale_fail &beta );

// ./stan/math/prim/prob/beta_lcdf.hpp:40
template <typename T_y,
    typename T_scale_succ,
    typename T_scale_fail>
return_type_t <T_y ,T_scale_succ ,T_scale_fail >beta_lcdf (const T_y &y , const T_scale_succ &alpha , const T_scale_fail &beta_param );

// ./stan/math/prim/prob/beta_cdf_log.hpp:13
template <typename T_y,
    typename T_scale_succ,
    typename T_scale_fail>
return_type_t <T_y ,T_scale_succ ,T_scale_fail >beta_cdf_log (const T_y &y , const T_scale_succ &alpha , const T_scale_fail &beta );

// ./stan/math/prim/prob/beta_lpdf.hpp:43
template <bool propto,
    typename T_y,
    typename T_scale_succ,
    typename T_scale_fail,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale_succ,
    T_scale_fail>*>
return_type_t <T_y ,T_scale_succ ,T_scale_fail >beta_lpdf (const T_y &y , const T_scale_succ &alpha , const T_scale_fail &beta );

// ./stan/math/prim/prob/beta_lpdf.hpp:122
template <typename T_y,
    typename T_scale_succ,
    typename T_scale_fail>
inline return_type_t <T_y ,T_scale_succ ,T_scale_fail >beta_lpdf (const T_y &y , const T_scale_succ &alpha , const T_scale_fail &beta );

// ./stan/math/prim/prob/beta_proportion_lccdf.hpp:43
template <typename T_y,
    typename T_loc,
    typename T_prec>
return_type_t <T_y ,T_loc ,T_prec >beta_proportion_lccdf (const T_y &y , const T_loc &mu , const T_prec &kappa );

// ./stan/math/prim/prob/beta_proportion_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_prec>
return_type_t <T_y ,T_loc ,T_prec >beta_proportion_ccdf_log (const T_y &y , const T_loc &mu , const T_prec &kappa );

// ./stan/math/prim/prob/beta_proportion_lcdf.hpp:44
template <typename T_y,
    typename T_loc,
    typename T_prec>
return_type_t <T_y ,T_loc ,T_prec >beta_proportion_lcdf (const T_y &y , const T_loc &mu , const T_prec &kappa );

// ./stan/math/prim/prob/beta_proportion_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_prec>
return_type_t <T_y ,T_loc ,T_prec >beta_proportion_cdf_log (const T_y &y , const T_loc &mu , const T_prec &kappa );

// ./stan/math/prim/prob/beta_proportion_lpdf.hpp:48
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_prec,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_prec>*>
return_type_t <T_y ,T_loc ,T_prec >beta_proportion_lpdf (const T_y &y , const T_loc &mu , const T_prec &kappa );

// ./stan/math/prim/prob/beta_proportion_lpdf.hpp:126
template <typename T_y,
    typename T_loc,
    typename T_prec>
inline return_type_t <T_y ,T_loc ,T_prec >beta_proportion_lpdf (const T_y &y , const T_loc &mu , const T_prec &kappa );

// ./stan/math/prim/prob/beta_proportion_rng.hpp:34
template <typename T_loc,
    typename T_prec,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_prec >::type beta_proportion_rng (const T_loc &mu , const T_prec &kappa , RNG &rng );

// ./stan/math/prim/prob/binomial_lccdf.hpp:37
template <typename T_n,
    typename T_N,
    typename T_prob>
return_type_t <T_prob >binomial_lccdf (const T_n &n , const T_N &N , const T_prob &theta );

// ./stan/math/prim/prob/binomial_ccdf_log.hpp:13
template <typename T_n,
    typename T_N,
    typename T_prob>
return_type_t <T_prob >binomial_ccdf_log (const T_n &n , const T_N &N , const T_prob &theta );

// ./stan/math/prim/prob/binomial_cdf.hpp:35
template <typename T_n,
    typename T_N,
    typename T_prob>
return_type_t <T_prob >binomial_cdf (const T_n &n , const T_N &N , const T_prob &theta );

// ./stan/math/prim/prob/binomial_lcdf.hpp:37
template <typename T_n,
    typename T_N,
    typename T_prob>
return_type_t <T_prob >binomial_lcdf (const T_n &n , const T_N &N , const T_prob &theta );

// ./stan/math/prim/prob/binomial_cdf_log.hpp:13
template <typename T_n,
    typename T_N,
    typename T_prob>
return_type_t <T_prob >binomial_cdf_log (const T_n &n , const T_N &N , const T_prob &theta );

// ./stan/math/prim/prob/binomial_logit_lpmf.hpp:34
template <bool propto,
    typename T_n,
    typename T_N,
    typename T_prob,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_N,
    T_prob>*>
return_type_t <T_prob >binomial_logit_lpmf (const T_n &n , const T_N &N , const T_prob &alpha );

// ./stan/math/prim/prob/binomial_logit_lpmf.hpp:87
template <typename T_n,
    typename T_N,
    typename T_prob>
inline return_type_t <T_prob >binomial_logit_lpmf (const T_n &n , const T_N &N , const T_prob &alpha );

// ./stan/math/prim/prob/binomial_logit_glm_lpmf.hpp:54
template <bool propto,
    typename T_n,
    typename T_N,
    typename T_x,
    typename T_alpha,
    typename T_beta,
    require_matrix_t<T_x>*>
return_type_t <T_x ,T_alpha ,T_beta >binomial_logit_glm_lpmf (const T_n &n , const T_N &N , const T_x &x , const T_alpha &alpha , const T_beta &beta );

// ./stan/math/prim/prob/binomial_logit_glm_lpmf.hpp:162
template <typename T_n,
    typename T_N,
    typename T_x,
    typename T_alpha,
    typename T_beta>
inline return_type_t <T_x ,T_beta ,T_alpha >binomial_logit_glm_lpmf (const T_n &n , const T_N &N , const T_x &x , const T_alpha &alpha , const T_beta &beta );

// ./stan/math/prim/prob/binomial_lpmf.hpp:36
template <bool propto,
    typename T_n,
    typename T_N,
    typename T_prob,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_N,
    T_prob>*>
return_type_t <T_prob >binomial_lpmf (const T_n &n , const T_N &N , const T_prob &theta );

// ./stan/math/prim/prob/binomial_lpmf.hpp:139
template <typename T_n,
    typename T_N,
    typename T_prob>
inline return_type_t <T_prob >binomial_lpmf (const T_n &n , const T_N &N , const T_prob &theta );

// ./stan/math/prim/prob/categorical_logit_glm_lpmf.hpp:43
template <bool propto,
    typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta,
    require_matrix_t<T_x>*>
return_type_t <T_x ,T_alpha ,T_beta >categorical_logit_glm_lpmf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta );

// ./stan/math/prim/prob/categorical_logit_glm_lpmf.hpp:197
template <typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta>
return_type_t <T_x ,T_alpha ,T_beta >categorical_logit_glm_lpmf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta );

// ./stan/math/prim/prob/categorical_logit_lpmf.hpp:16
template <bool propto,
    typename T_prob,
    require_col_vector_t<T_prob>*>
return_type_t <T_prob >categorical_logit_lpmf (int n , const T_prob &beta );

// ./stan/math/prim/prob/categorical_logit_lpmf.hpp:33
template <bool propto,
    typename T_prob,
    require_col_vector_t<T_prob>*>
return_type_t <T_prob >categorical_logit_lpmf (const std ::vector <int >&ns , const T_prob &beta );

// ./stan/math/prim/prob/categorical_logit_lpmf.hpp:61
template <typename T_n,
    typename T_prob,
    require_st_integral<T_n>*>
inline return_type_t <T_prob >categorical_logit_lpmf (const T_n &ns , const T_prob &beta );

// ./stan/math/prim/prob/categorical_logit_rng.hpp:28
template <class RNG>
inline int categorical_logit_rng (const Eigen ::VectorXd &beta , RNG &rng );

// ./stan/math/prim/prob/categorical_lpmf.hpp:15
template <bool propto,
    typename T_prob,
    require_eigen_col_vector_t<T_prob>*>
return_type_t <T_prob >categorical_lpmf (int n , const T_prob &theta );

// ./stan/math/prim/prob/categorical_lpmf.hpp:31
template <bool propto,
    typename T_prob,
    require_eigen_col_vector_t<T_prob>*>
return_type_t <T_prob >categorical_lpmf (const std ::vector <int >&ns , const T_prob &theta );

// ./stan/math/prim/prob/categorical_lpmf.hpp:60
template <typename T_n,
    typename T_prob,
    require_st_integral<T_n>*>
inline return_type_t <T_prob >categorical_lpmf (const T_n &ns , const T_prob &theta );

// ./stan/math/prim/prob/categorical_rng.hpp:13
template <class RNG>
inline int categorical_rng (const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&theta , RNG &rng );

// ./stan/math/prim/prob/cauchy_lccdf.hpp:33
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >cauchy_lccdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/cauchy_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >cauchy_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/cauchy_cdf.hpp:32
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >cauchy_cdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/cauchy_lcdf.hpp:33
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >cauchy_lcdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/cauchy_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >cauchy_cdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/cauchy_lpdf.hpp:40
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >cauchy_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/cauchy_lpdf.hpp:118
template <typename T_y,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_loc ,T_scale >cauchy_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/cauchy_rng.hpp:32
template <typename T_loc,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale >::type cauchy_rng (const T_loc &mu , const T_scale &sigma , RNG &rng );

// ./stan/math/prim/prob/chi_square_lccdf.hpp:37
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >chi_square_lccdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/chi_square_ccdf_log.hpp:13
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >chi_square_ccdf_log (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/chi_square_cdf.hpp:36
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >chi_square_cdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/chi_square_lcdf.hpp:37
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >chi_square_lcdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/chi_square_cdf_log.hpp:13
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >chi_square_cdf_log (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/chi_square_lpdf.hpp:43
template <bool propto,
    typename T_y,
    typename T_dof,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_dof>*>
return_type_t <T_y ,T_dof >chi_square_lpdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/chi_square_lpdf.hpp:101
template <typename T_y,
    typename T_dof>
inline return_type_t <T_y ,T_dof >chi_square_lpdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/chi_square_rng.hpp:27
template <typename T_deg,
    class RNG>
inline typename VectorBuilder <true ,double ,T_deg >::type chi_square_rng (const T_deg &nu , RNG &rng );

// ./stan/math/prim/prob/dirichlet_lpdf.hpp:56
template <bool propto,
    typename T_prob,
    typename T_prior_size,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_prob,
    T_prior_size>*>
return_type_t <T_prob ,T_prior_size >dirichlet_lpdf (const T_prob &theta , const T_prior_size &alpha );

// ./stan/math/prim/prob/dirichlet_lpdf.hpp:128
template <typename T_prob,
    typename T_prior_size>
return_type_t <T_prob ,T_prior_size >dirichlet_lpdf (const T_prob &theta , const T_prior_size &alpha );

// ./stan/math/prim/prob/dirichlet_rng.hpp:38
template <class RNG>
inline Eigen ::VectorXd dirichlet_rng (const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&alpha , RNG &rng );

// ./stan/math/prim/prob/dirichlet_multinomial_lpmf.hpp:58
template <bool propto,
    typename T_prior_size,
    require_eigen_col_vector_t<T_prior_size>*>
return_type_t <T_prior_size >dirichlet_multinomial_lpmf (const std ::vector <int >&ns , const T_prior_size &alpha );

// ./stan/math/prim/prob/dirichlet_multinomial_lpmf.hpp:100
template <typename T_prior_size>
return_type_t <T_prior_size >dirichlet_multinomial_lpmf (const std ::vector <int >&ns , const T_prior_size &alpha );

// ./stan/math/prim/prob/multinomial_rng.hpp:26
template <class T_theta,
    class RNG,
    require_eigen_col_vector_t<T_theta>*>
inline std ::vector <int >multinomial_rng (const T_theta &theta , int N , RNG &rng );

// ./stan/math/prim/prob/dirichlet_multinomial_rng.hpp:36
template <class RNG>
inline std ::vector <int >dirichlet_multinomial_rng (const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&alpha , int N , RNG &rng );

// ./stan/math/prim/prob/discrete_range_lccdf.hpp:37
template <typename T_y,
    typename T_lower,
    typename T_upper>
double discrete_range_lccdf (const T_y &y , const T_lower &lower , const T_upper &upper );

// ./stan/math/prim/prob/discrete_range_ccdf_log.hpp:12
template <typename T_y,
    typename T_lower,
    typename T_upper>
double discrete_range_ccdf_log (const T_y &y , const T_lower &lower , const T_upper &upper );

// ./stan/math/prim/prob/discrete_range_cdf.hpp:37
template <typename T_y,
    typename T_lower,
    typename T_upper>
double discrete_range_cdf (const T_y &y , const T_lower &lower , const T_upper &upper );

// ./stan/math/prim/prob/discrete_range_lcdf.hpp:37
template <typename T_y,
    typename T_lower,
    typename T_upper>
double discrete_range_lcdf (const T_y &y , const T_lower &lower , const T_upper &upper );

// ./stan/math/prim/prob/discrete_range_cdf_log.hpp:12
template <typename T_y,
    typename T_lower,
    typename T_upper>
double discrete_range_cdf_log (const T_y &y , const T_lower &lower , const T_upper &upper );

// ./stan/math/prim/prob/discrete_range_lpmf.hpp:45
template <bool propto,
    typename T_y,
    typename T_lower,
    typename T_upper>
double discrete_range_lpmf (const T_y &y , const T_lower &lower , const T_upper &upper );

// ./stan/math/prim/prob/discrete_range_lpmf.hpp:92
template <typename T_y,
    typename T_lower,
    typename T_upper>
inline double discrete_range_lpmf (const T_y &y , const T_lower &lower , const T_upper &upper );

// ./stan/math/prim/prob/discrete_range_rng.hpp:34
template <typename T_lower,
    typename T_upper,
    class RNG>
inline typename VectorBuilder <true ,int ,T_lower ,T_upper >::type discrete_range_rng (const T_lower &lower , const T_upper &upper , RNG &rng );

// ./stan/math/prim/prob/double_exponential_lccdf.hpp:35
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >double_exponential_lccdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/double_exponential_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >double_exponential_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/double_exponential_cdf.hpp:36
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >double_exponential_cdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/double_exponential_lcdf.hpp:34
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >double_exponential_lcdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/double_exponential_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >double_exponential_cdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/double_exponential_lpdf.hpp:39
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >double_exponential_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/double_exponential_lpdf.hpp:108
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >double_exponential_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/double_exponential_rng.hpp:34
template <typename T_loc,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale >::type double_exponential_rng (const T_loc &mu , const T_scale &sigma , RNG &rng );

// ./stan/math/prim/prob/exp_mod_normal_lccdf.hpp:28
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_inv_scale>*>
return_type_t <T_y ,T_loc ,T_scale ,T_inv_scale >exp_mod_normal_lccdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_inv_scale &lambda );

// ./stan/math/prim/prob/exp_mod_normal_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_inv_scale>
return_type_t <T_y ,T_loc ,T_scale ,T_inv_scale >exp_mod_normal_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_inv_scale &lambda );

// ./stan/math/prim/prob/exp_mod_normal_cdf.hpp:26
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_inv_scale>*>
return_type_t <T_y ,T_loc ,T_scale ,T_inv_scale >exp_mod_normal_cdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_inv_scale &lambda );

// ./stan/math/prim/prob/exp_mod_normal_lcdf.hpp:28
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_inv_scale>*>
return_type_t <T_y ,T_loc ,T_scale ,T_inv_scale >exp_mod_normal_lcdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_inv_scale &lambda );

// ./stan/math/prim/prob/exp_mod_normal_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_inv_scale>
return_type_t <T_y ,T_loc ,T_scale ,T_inv_scale >exp_mod_normal_cdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_inv_scale &lambda );

// ./stan/math/prim/prob/exp_mod_normal_lpdf.hpp:25
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_inv_scale>*>
return_type_t <T_y ,T_loc ,T_scale ,T_inv_scale >exp_mod_normal_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_inv_scale &lambda );

// ./stan/math/prim/prob/exp_mod_normal_lpdf.hpp:121
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_inv_scale>
inline return_type_t <T_y ,T_loc ,T_scale ,T_inv_scale >exp_mod_normal_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_inv_scale &lambda );

// ./stan/math/prim/prob/exponential_rng.hpp:27
template <typename T_inv,
    class RNG>
inline typename VectorBuilder <true ,double ,T_inv >::type exponential_rng (const T_inv &beta , RNG &rng );

// ./stan/math/prim/prob/normal_rng.hpp:32
template <typename T_loc,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale >::type normal_rng (const T_loc &mu , const T_scale &sigma , RNG &rng );

// ./stan/math/prim/prob/exp_mod_normal_rng.hpp:38
template <typename T_loc,
    typename T_scale,
    typename T_inv_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale ,T_inv_scale >::type exp_mod_normal_rng (const T_loc &mu , const T_scale &sigma , const T_inv_scale &lambda , RNG &rng );

// ./stan/math/prim/prob/exponential_lccdf.hpp:17
template <typename T_y,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_inv_scale>*>
return_type_t <T_y ,T_inv_scale >exponential_lccdf (const T_y &y , const T_inv_scale &beta );

// ./stan/math/prim/prob/exponential_ccdf_log.hpp:13
template <typename T_y,
    typename T_inv_scale>
return_type_t <T_y ,T_inv_scale >exponential_ccdf_log (const T_y &y , const T_inv_scale &beta );

// ./stan/math/prim/prob/exponential_cdf.hpp:32
template <typename T_y,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_inv_scale>*>
return_type_t <T_y ,T_inv_scale >exponential_cdf (const T_y &y , const T_inv_scale &beta );

// ./stan/math/prim/prob/exponential_lcdf.hpp:21
template <typename T_y,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_inv_scale>*>
return_type_t <T_y ,T_inv_scale >exponential_lcdf (const T_y &y , const T_inv_scale &beta );

// ./stan/math/prim/prob/exponential_cdf_log.hpp:13
template <typename T_y,
    typename T_inv_scale>
return_type_t <T_y ,T_inv_scale >exponential_cdf_log (const T_y &y , const T_inv_scale &beta );

// ./stan/math/prim/prob/exponential_lpdf.hpp:49
template <bool propto,
    typename T_y,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_inv_scale>*>
return_type_t <T_y ,T_inv_scale >exponential_lpdf (const T_y &y , const T_inv_scale &beta );

// ./stan/math/prim/prob/exponential_lpdf.hpp:104
template <typename T_y,
    typename T_inv_scale>
inline return_type_t <T_y ,T_inv_scale >exponential_lpdf (const T_y &y , const T_inv_scale &beta );

// ./stan/math/prim/prob/frechet_lccdf.hpp:25
template <typename T_y,
    typename T_shape,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_scale>*>
return_type_t <T_y ,T_shape ,T_scale >frechet_lccdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/frechet_ccdf_log.hpp:13
template <typename T_y,
    typename T_shape,
    typename T_scale>
return_type_t <T_y ,T_shape ,T_scale >frechet_ccdf_log (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/frechet_cdf.hpp:25
template <typename T_y,
    typename T_shape,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_scale>*>
return_type_t <T_y ,T_shape ,T_scale >frechet_cdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/frechet_lcdf.hpp:23
template <typename T_y,
    typename T_shape,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_scale>*>
return_type_t <T_y ,T_shape ,T_scale >frechet_lcdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/frechet_cdf_log.hpp:13
template <typename T_y,
    typename T_shape,
    typename T_scale>
return_type_t <T_y ,T_shape ,T_scale >frechet_cdf_log (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/frechet_lpdf.hpp:28
template <bool propto,
    typename T_y,
    typename T_shape,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_scale>*>
return_type_t <T_y ,T_shape ,T_scale >frechet_lpdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/frechet_lpdf.hpp:99
template <typename T_y,
    typename T_shape,
    typename T_scale>
inline return_type_t <T_y ,T_shape ,T_scale >frechet_lpdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/frechet_rng.hpp:31
template <typename T_shape,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_shape ,T_scale >::type frechet_rng (const T_shape &alpha , const T_scale &sigma , RNG &rng );

// ./stan/math/prim/prob/gamma_lccdf.hpp:24
template <typename T_y,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_y ,T_shape ,T_inv_scale >gamma_lccdf (const T_y &y , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/gamma_ccdf_log.hpp:13
template <typename T_y,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_y ,T_shape ,T_inv_scale >gamma_ccdf_log (const T_y &y , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/gamma_cdf.hpp:40
template <typename T_y,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_y ,T_shape ,T_inv_scale >gamma_cdf (const T_y &y , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/gamma_lcdf.hpp:24
template <typename T_y,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_y ,T_shape ,T_inv_scale >gamma_lcdf (const T_y &y , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/gamma_cdf_log.hpp:13
template <typename T_y,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_y ,T_shape ,T_inv_scale >gamma_cdf_log (const T_y &y , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/gamma_lpdf.hpp:49
template <bool propto,
    typename T_y,
    typename T_shape,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_inv_scale>*>
return_type_t <T_y ,T_shape ,T_inv_scale >gamma_lpdf (const T_y &y , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/gamma_lpdf.hpp:120
template <typename T_y,
    typename T_shape,
    typename T_inv_scale>
inline return_type_t <T_y ,T_shape ,T_inv_scale >gamma_lpdf (const T_y &y , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/gamma_rng.hpp:32
template <typename T_shape,
    typename T_inv,
    class RNG>
inline typename VectorBuilder <true ,double ,T_shape ,T_inv >::type gamma_rng (const T_shape &alpha , const T_inv &beta , RNG &rng );

// ./stan/math/prim/prob/gaussian_dlm_obs_lpdf.hpp:65
template <bool propto,
    typename T_y,
    typename T_F,
    typename T_G,
    typename T_V,
    typename T_W,
    typename T_m0,
    typename T_C0,
    require_all_eigen_matrix_dynamic_t<T_y,
    T_F,
    T_G,
    T_V,
    T_W,
    T_C0>*>
inline return_type_t <T_y ,T_F ,T_G ,T_V ,T_W ,T_m0 ,T_C0 >gaussian_dlm_obs_lpdf (const T_y &y , const T_F &F , const T_G &G , const T_V &V , const T_W &W , const T_m0 &m0 , const T_C0 &C0 );

// ./stan/math/prim/prob/gaussian_dlm_obs_lpdf.hpp:190
template <bool propto,
    typename T_y,
    typename T_F,
    typename T_G,
    typename T_V,
    typename T_W,
    typename T_m0,
    typename T_C0,
    require_all_eigen_matrix_dynamic_t<T_y,
    T_F,
    T_G,
    T_W,
    T_C0>*>
inline return_type_t <T_y ,T_F ,T_G ,T_V ,T_W ,T_m0 ,T_C0 >gaussian_dlm_obs_lpdf (const T_y &y , const T_F &F , const T_G &G , const T_V &V , const T_W &W , const T_m0 &m0 , const T_C0 &C0 );

// ./stan/math/prim/prob/gaussian_dlm_obs_lpdf.hpp:286
template <typename T_y,
    typename T_F,
    typename T_G,
    typename T_V,
    typename T_W,
    typename T_m0,
    typename T_C0>
inline return_type_t <T_y ,T_F ,T_G ,T_V ,T_W ,T_m0 ,T_C0 >gaussian_dlm_obs_lpdf (const T_y &y , const T_F &F , const T_G &G , const T_V &V , const T_W &W , const T_m0 &m0 , const T_C0 &C0 );

// ./stan/math/prim/prob/gaussian_dlm_obs_rng.hpp:87
template <class RNG>
inline Eigen ::MatrixXd gaussian_dlm_obs_rng (const Eigen ::MatrixXd &F , const Eigen ::MatrixXd &G , const Eigen ::MatrixXd &V , const Eigen ::MatrixXd &W , const Eigen ::VectorXd &m0 , const Eigen ::MatrixXd &C0 , const int T , RNG &rng );

// ./stan/math/prim/prob/gaussian_dlm_obs_rng.hpp:140
template <class RNG>
inline Eigen ::Matrix <double ,Eigen ::Dynamic ,Eigen ::Dynamic >gaussian_dlm_obs_rng (const Eigen ::MatrixXd &F , const Eigen ::MatrixXd &G , const Eigen ::VectorXd &V , const Eigen ::MatrixXd &W , const Eigen ::VectorXd &m0 , const Eigen ::MatrixXd &C0 , const int T , RNG &rng );

// ./stan/math/prim/prob/gumbel_lccdf.hpp:36
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >gumbel_lccdf (const T_y &y , const T_loc &mu , const T_scale &beta );

// ./stan/math/prim/prob/gumbel_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >gumbel_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &beta );

// ./stan/math/prim/prob/gumbel_cdf.hpp:36
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >gumbel_cdf (const T_y &y , const T_loc &mu , const T_scale &beta );

// ./stan/math/prim/prob/gumbel_lcdf.hpp:35
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >gumbel_lcdf (const T_y &y , const T_loc &mu , const T_scale &beta );

// ./stan/math/prim/prob/gumbel_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >gumbel_cdf_log (const T_y &y , const T_loc &mu , const T_scale &beta );

// ./stan/math/prim/prob/gumbel_lpdf.hpp:37
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >gumbel_lpdf (const T_y &y , const T_loc &mu , const T_scale &beta );

// ./stan/math/prim/prob/gumbel_lpdf.hpp:103
template <typename T_y,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_loc ,T_scale >gumbel_lpdf (const T_y &y , const T_loc &mu , const T_scale &beta );

// ./stan/math/prim/prob/gumbel_rng.hpp:33
template <typename T_loc,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale >::type gumbel_rng (const T_loc &mu , const T_scale &beta , RNG &rng );

// ./stan/math/prim/prob/hmm_hidden_state_prob.hpp:45
template <typename T_omega,
    typename T_Gamma,
    typename T_rho,
    require_all_eigen_t<T_omega,
    T_Gamma>*>
inline Eigen ::MatrixXd hmm_hidden_state_prob (const T_omega &log_omegas , const T_Gamma &Gamma , const T_rho &rho );

// ./stan/math/prim/prob/hmm_latent_rng.hpp:41
template <typename T_omega,
    typename T_Gamma,
    typename T_rho,
    class RNG,
    require_all_eigen_t<T_omega,
    T_Gamma>*>
inline std ::vector <int >hmm_latent_rng (const T_omega &log_omegas , const T_Gamma &Gamma , const T_rho &rho , RNG &rng );

// ./stan/math/prim/prob/hmm_marginal.hpp:19
template <typename T_omega,
    typename T_Gamma,
    typename T_rho,
    typename T_alpha>
inline auto hmm_marginal_val (const Eigen ::Matrix <T_omega , Eigen ::Dynamic , Eigen ::Dynamic >&omegas , const T_Gamma &Gamma_val , const T_rho &rho_val , Eigen ::Matrix <T_alpha , Eigen ::Dynamic , Eigen ::Dynamic >&alphas , Eigen ::Matrix <T_alpha , Eigen ::Dynamic , 1>&alpha_log_norms , T_alpha &norm_norm );

// ./stan/math/prim/prob/hmm_marginal.hpp:73
template <typename T_omega,
    typename T_Gamma,
    typename T_rho,
    require_all_eigen_t<T_omega,
    T_Gamma>*>
inline auto hmm_marginal (const T_omega &log_omegas , const T_Gamma &Gamma , const T_rho &rho );

// ./stan/math/prim/prob/hypergeometric_lpmf.hpp:17
template <bool propto,
    typename T_n,
    typename T_N,
    typename T_a,
    typename T_b>
double hypergeometric_lpmf (const T_n &n , const T_N &N , const T_a &a , const T_b &b );

// ./stan/math/prim/prob/hypergeometric_lpmf.hpp:58
template <typename T_n,
    typename T_N,
    typename T_a,
    typename T_b>
inline double hypergeometric_lpmf (const T_n &n , const T_N &N , const T_a &a , const T_b &b );

// ./stan/math/prim/prob/uniform_rng.hpp:34
template <typename T_alpha,
    typename T_beta,
    class RNG>
inline typename VectorBuilder <true ,double ,T_alpha ,T_beta >::type uniform_rng (const T_alpha &alpha , const T_beta &beta , RNG &rng );

// ./stan/math/prim/prob/hypergeometric_rng.hpp:13
template <class RNG>
inline int hypergeometric_rng (int N , int a , int b , RNG &rng );

// ./stan/math/prim/prob/inv_chi_square_lccdf.hpp:37
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >inv_chi_square_lccdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/inv_chi_square_ccdf_log.hpp:13
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >inv_chi_square_ccdf_log (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/inv_chi_square_cdf.hpp:36
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >inv_chi_square_cdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/inv_chi_square_lcdf.hpp:37
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >inv_chi_square_lcdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/inv_chi_square_cdf_log.hpp:13
template <typename T_y,
    typename T_dof>
return_type_t <T_y ,T_dof >inv_chi_square_cdf_log (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/inv_chi_square_lpdf.hpp:46
template <bool propto,
    typename T_y,
    typename T_dof,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_dof>*>
return_type_t <T_y ,T_dof >inv_chi_square_lpdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/inv_chi_square_lpdf.hpp:103
template <typename T_y,
    typename T_dof>
inline return_type_t <T_y ,T_dof >inv_chi_square_lpdf (const T_y &y , const T_dof &nu );

// ./stan/math/prim/prob/inv_chi_square_rng.hpp:27
template <typename T_deg,
    class RNG>
inline typename VectorBuilder <true ,double ,T_deg >::type inv_chi_square_rng (const T_deg &nu , RNG &rng );

// ./stan/math/prim/prob/inv_gamma_lccdf.hpp:24
template <typename T_y,
    typename T_shape,
    typename T_scale>
return_type_t <T_y ,T_shape ,T_scale >inv_gamma_lccdf (const T_y &y , const T_shape &alpha , const T_scale &beta );

// ./stan/math/prim/prob/inv_gamma_ccdf_log.hpp:13
template <typename T_y,
    typename T_shape,
    typename T_scale>
return_type_t <T_y ,T_shape ,T_scale >inv_gamma_ccdf_log (const T_y &y , const T_shape &alpha , const T_scale &beta );

// ./stan/math/prim/prob/inv_gamma_cdf.hpp:39
template <typename T_y,
    typename T_shape,
    typename T_scale>
return_type_t <T_y ,T_shape ,T_scale >inv_gamma_cdf (const T_y &y , const T_shape &alpha , const T_scale &beta );

// ./stan/math/prim/prob/inv_gamma_lcdf.hpp:24
template <typename T_y,
    typename T_shape,
    typename T_scale>
return_type_t <T_y ,T_shape ,T_scale >inv_gamma_lcdf (const T_y &y , const T_shape &alpha , const T_scale &beta );

// ./stan/math/prim/prob/inv_gamma_cdf_log.hpp:13
template <typename T_y,
    typename T_shape,
    typename T_scale>
return_type_t <T_y ,T_shape ,T_scale >inv_gamma_cdf_log (const T_y &y , const T_shape &alpha , const T_scale &beta );

// ./stan/math/prim/prob/inv_gamma_lpdf.hpp:41
template <bool propto,
    typename T_y,
    typename T_shape,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_scale>*>
return_type_t <T_y ,T_shape ,T_scale >inv_gamma_lpdf (const T_y &y , const T_shape &alpha , const T_scale &beta );

// ./stan/math/prim/prob/inv_gamma_lpdf.hpp:114
template <typename T_y,
    typename T_shape,
    typename T_scale>
inline return_type_t <T_y ,T_shape ,T_scale >inv_gamma_lpdf (const T_y &y , const T_shape &alpha , const T_scale &beta );

// ./stan/math/prim/prob/inv_gamma_rng.hpp:32
template <typename T_shape,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_shape ,T_scale >::type inv_gamma_rng (const T_shape &alpha , const T_scale &beta , RNG &rng );

// ./stan/math/prim/prob/inv_wishart_cholesky_lpdf.hpp:40
template <bool propto,
    typename T_y,
    typename T_dof,
    typename T_scale,
    require_stan_scalar_t<T_dof>*>
return_type_t <T_y ,T_dof ,T_scale >inv_wishart_cholesky_lpdf (const T_y &L_Y , const T_dof &nu , const T_scale &L_S );

// ./stan/math/prim/prob/inv_wishart_cholesky_lpdf.hpp:92
template <typename T_y,
    typename T_dof,
    typename T_scale>
inline return_type_t <T_y ,T_dof ,T_scale >inv_wishart_cholesky_lpdf (const T_y &L_Y , const T_dof &nu , const T_scale &L_S );

// ./stan/math/prim/prob/inv_wishart_cholesky_rng.hpp:30
template <class RNG>
inline Eigen ::MatrixXd inv_wishart_cholesky_rng (double nu , const Eigen ::MatrixXd &L_S , RNG &rng );

// ./stan/math/prim/prob/inv_wishart_lpdf.hpp:44
template <bool propto,
    typename T_y,
    typename T_dof,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_scale >inv_wishart_lpdf (const T_y &W , const T_dof &nu , const T_scale &S );

// ./stan/math/prim/prob/inv_wishart_lpdf.hpp:91
template <typename T_y,
    typename T_dof,
    typename T_scale>
inline return_type_t <T_y ,T_dof ,T_scale >inv_wishart_lpdf (const T_y &W , const T_dof &nu , const T_scale &S );

// ./stan/math/prim/prob/wishart_rng.hpp:14
template <class RNG>
inline Eigen ::MatrixXd wishart_rng (double nu , const Eigen ::MatrixXd &S , RNG &rng );

// ./stan/math/prim/prob/inv_wishart_rng.hpp:12
template <class RNG>
inline Eigen ::MatrixXd inv_wishart_rng (double nu , const Eigen ::MatrixXd &S , RNG &rng );

// ./stan/math/prim/prob/lkj_corr_lpdf.hpp:15
template <typename T_shape>
return_type_t <double ,T_shape >do_lkj_constant (const T_shape &eta , const unsigned int &K );

// ./stan/math/prim/prob/lkj_corr_lpdf.hpp:48
template <bool propto,
    typename T_y,
    typename T_shape>
return_type_t <T_y ,T_shape >lkj_corr_lpdf (const T_y &y , const T_shape &eta );

// ./stan/math/prim/prob/lkj_corr_lpdf.hpp:79
template <typename T_y,
    typename T_shape>
inline return_type_t <T_y ,T_shape >lkj_corr_lpdf (const T_y &y , const T_shape &eta );

// ./stan/math/prim/prob/lkj_corr_cholesky_lpdf.hpp:17
template <bool propto,
    typename T_covar,
    typename T_shape>
return_type_t <T_covar ,T_shape >lkj_corr_cholesky_lpdf (const T_covar &L , const T_shape &eta );

// ./stan/math/prim/prob/lkj_corr_cholesky_lpdf.hpp:56
template <typename T_covar,
    typename T_shape>
inline return_type_t <T_covar ,T_shape >lkj_corr_cholesky_lpdf (const T_covar &L , const T_shape &eta );

// ./stan/math/prim/prob/lkj_corr_cholesky_rng.hpp:12
template <class RNG>
inline Eigen ::MatrixXd lkj_corr_cholesky_rng (size_t K , double eta , RNG &rng );

// ./stan/math/prim/prob/lkj_corr_rng.hpp:25
template <class RNG>
inline Eigen ::MatrixXd lkj_corr_rng (size_t K , double eta , RNG &rng );

// ./stan/math/prim/prob/lognormal_lpdf.hpp:25
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >lognormal_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/lognormal_lpdf.hpp:102
template <typename T_y,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_loc ,T_scale >lognormal_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/lkj_cov_lpdf.hpp:18
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape,
    require_eigen_matrix_dynamic_t<T_y>*>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >lkj_cov_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &eta );

// ./stan/math/prim/prob/lkj_cov_lpdf.hpp:63
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape,
    require_eigen_matrix_dynamic_t<T_y>*>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >lkj_cov_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &eta );

// ./stan/math/prim/prob/lkj_cov_lpdf.hpp:99
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
inline return_type_t <T_y ,T_loc ,T_scale ,T_shape >lkj_cov_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &eta );

// ./stan/math/prim/prob/logistic_lpdf.hpp:25
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >logistic_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/logistic_lpdf.hpp:94
template <typename T_y,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_loc ,T_scale >logistic_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/logistic_lccdf.hpp:22
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >logistic_lccdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/logistic_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >logistic_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/logistic_cdf.hpp:22
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >logistic_cdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/logistic_lcdf.hpp:22
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >logistic_lcdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/logistic_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >logistic_cdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/logistic_rng.hpp:33
template <typename T_loc,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale >::type logistic_rng (const T_loc &mu , const T_scale &sigma , RNG &rng );

// ./stan/math/prim/prob/loglogistic_cdf.hpp:43
template <typename T_y,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_scale ,T_shape >loglogistic_cdf (const T_y &y , const T_scale &alpha , const T_shape &beta );

// ./stan/math/prim/prob/loglogistic_lpdf.hpp:42
template <bool propto,
    typename T_y,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_scale ,T_shape >loglogistic_lpdf (const T_y &y , const T_scale &alpha , const T_shape &beta );

// ./stan/math/prim/prob/loglogistic_lpdf.hpp:143
template <typename T_y,
    typename T_scale,
    typename T_shape>
inline return_type_t <T_y ,T_scale ,T_shape >loglogistic_lpdf (const T_y &y , const T_scale &alpha , const T_shape &beta );

// ./stan/math/prim/prob/loglogistic_rng.hpp:32
template <typename T_scale,
    typename T_shape,
    class RNG>
inline typename VectorBuilder <true ,double ,T_scale ,T_shape >::type loglogistic_rng (const T_scale &alpha , const T_shape &beta , RNG &rng );

// ./stan/math/prim/prob/lognormal_lccdf.hpp:25
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >lognormal_lccdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/lognormal_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >lognormal_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/lognormal_cdf.hpp:25
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >lognormal_cdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/lognormal_lcdf.hpp:25
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_loc ,T_scale >lognormal_lcdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/lognormal_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >lognormal_cdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/lognormal_rng.hpp:32
template <typename T_loc,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale >::type lognormal_rng (const T_loc &mu , const T_scale &sigma , RNG &rng );

// ./stan/math/prim/prob/matrix_normal_prec_lpdf.hpp:33
template <bool propto,
    typename T_y,
    typename T_Mu,
    typename T_Sigma,
    typename T_D,
    require_all_matrix_t<T_y,
    T_Mu,
    T_Sigma,
    T_D>*>
return_type_t <T_y ,T_Mu ,T_Sigma ,T_D >matrix_normal_prec_lpdf (const T_y &y , const T_Mu &Mu , const T_Sigma &Sigma , const T_D &D );

// ./stan/math/prim/prob/matrix_normal_prec_lpdf.hpp:82
template <typename T_y,
    typename T_Mu,
    typename T_Sigma,
    typename T_D,
    require_all_matrix_t<T_y,
    T_Mu,
    T_Sigma,
    T_D>*>
return_type_t <T_y ,T_Mu ,T_Sigma ,T_D >matrix_normal_prec_lpdf (const T_y &y , const T_Mu &Mu , const T_Sigma &Sigma , const T_D &D );

// ./stan/math/prim/prob/matrix_normal_prec_rng.hpp:31
template <class RNG>
inline Eigen ::MatrixXd matrix_normal_prec_rng (const Eigen ::MatrixXd &Mu , const Eigen ::MatrixXd &Sigma , const Eigen ::MatrixXd &D , RNG &rng );

// ./stan/math/prim/prob/multi_gp_cholesky_lpdf.hpp:38
template <bool propto,
    typename T_y,
    typename T_covar,
    typename T_w,
    require_all_eigen_matrix_dynamic_t<T_y,
    T_covar>*>
return_type_t <T_y ,T_covar ,T_w >multi_gp_cholesky_lpdf (const T_y &y , const T_covar &L , const T_w &w );

// ./stan/math/prim/prob/multi_gp_cholesky_lpdf.hpp:89
template <typename T_y,
    typename T_covar,
    typename T_w>
inline return_type_t <T_y ,T_covar ,T_w >multi_gp_cholesky_lpdf (const T_y &y , const T_covar &L , const T_w &w );

// ./stan/math/prim/prob/multi_gp_lpdf.hpp:33
template <bool propto,
    typename T_y,
    typename T_covar,
    typename T_w,
    require_all_matrix_t<T_y,
    T_covar>*>
return_type_t <T_y ,T_covar ,T_w >multi_gp_lpdf (const T_y &y , const T_covar &Sigma , const T_w &w );

// ./stan/math/prim/prob/multi_gp_lpdf.hpp:84
template <typename T_y,
    typename T_covar,
    typename T_w>
inline return_type_t <T_y ,T_covar ,T_w >multi_gp_lpdf (const T_y &y , const T_covar &Sigma , const T_w &w );

// ./stan/math/prim/prob/multi_normal_cholesky_lpdf.hpp:44
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_covar,
    require_any_not_vector_vt<is_stan_scalar,
    T_y,
    T_loc>*>
return_type_t <T_y ,T_loc ,T_covar >multi_normal_cholesky_lpdf (const T_y &y , const T_loc &mu , const T_covar &L );

// ./stan/math/prim/prob/multi_normal_cholesky_lpdf.hpp:203
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_covar,
    require_all_vector_vt<is_stan_scalar,
    T_y,
    T_loc>*>
return_type_t <T_y ,T_loc ,T_covar >multi_normal_cholesky_lpdf (const T_y &y , const T_loc &mu , const T_covar &L );

// ./stan/math/prim/prob/multi_normal_cholesky_lpdf.hpp:299
template <typename T_y,
    typename T_loc,
    typename T_covar>
inline return_type_t <T_y ,T_loc ,T_covar >multi_normal_cholesky_lpdf (const T_y &y , const T_loc &mu , const T_covar &L );

// ./stan/math/prim/prob/multi_normal_cholesky_rng.hpp:32
template <typename T_loc,
    class RNG>
inline typename StdVectorBuilder <true ,Eigen ::VectorXd ,T_loc >::type multi_normal_cholesky_rng (const T_loc &mu , const Eigen ::Matrix <double , Eigen ::Dynamic , Eigen ::Dynamic >&L , RNG &rng );

// ./stan/math/prim/prob/multi_normal_lpdf.hpp:24
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_covar,
    require_any_not_vector_vt<is_stan_scalar,
    T_y,
    T_loc>*>
return_type_t <T_y ,T_loc ,T_covar >multi_normal_lpdf (const T_y &y , const T_loc &mu , const T_covar &Sigma );

// ./stan/math/prim/prob/multi_normal_lpdf.hpp:158
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_covar,
    require_all_vector_vt<is_stan_scalar,
    T_y,
    T_loc>*>
return_type_t <T_y ,T_loc ,T_covar >multi_normal_lpdf (const T_y &y , const T_loc &mu , const T_covar &Sigma );

// ./stan/math/prim/prob/multi_normal_lpdf.hpp:254
template <typename T_y,
    typename T_loc,
    typename T_covar>
inline return_type_t <T_y ,T_loc ,T_covar >multi_normal_lpdf (const T_y &y , const T_loc &mu , const T_covar &Sigma );

// ./stan/math/prim/prob/multi_normal_prec_lpdf.hpp:19
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_covar>
return_type_t <T_y ,T_loc ,T_covar >multi_normal_prec_lpdf (const T_y &y , const T_loc &mu , const T_covar &Sigma );

// ./stan/math/prim/prob/multi_normal_prec_lpdf.hpp:105
template <typename T_y,
    typename T_loc,
    typename T_covar>
inline return_type_t <T_y ,T_loc ,T_covar >multi_normal_prec_lpdf (const T_y &y , const T_loc &mu , const T_covar &Sigma );

// ./stan/math/prim/prob/multi_normal_prec_rng.hpp:32
template <typename T_loc,
    class RNG>
inline typename StdVectorBuilder <true ,Eigen ::VectorXd ,T_loc >::type multi_normal_prec_rng (const T_loc &mu , const Eigen ::MatrixXd &S , RNG &rng );

// ./stan/math/prim/prob/multi_normal_rng.hpp:32
template <typename T_loc,
    class RNG>
inline typename StdVectorBuilder <true ,Eigen ::VectorXd ,T_loc >::type multi_normal_rng (const T_loc &mu , const Eigen ::Matrix <double , Eigen ::Dynamic , Eigen ::Dynamic >&S , RNG &rng );

// ./stan/math/prim/prob/multi_student_t_cholesky_lpdf.hpp:52
template <bool propto,
    typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_covar,
    require_any_not_vector_vt<is_stan_scalar,
    T_y,
    T_dof,
    T_loc>*>
return_type_t <T_y ,T_dof ,T_loc ,T_covar >multi_student_t_cholesky_lpdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_covar &L );

// ./stan/math/prim/prob/multi_student_t_cholesky_lpdf.hpp:239
template <bool propto,
    typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_covar,
    require_all_vector_vt<is_stan_scalar,
    T_y,
    T_dof,
    T_loc>*>
return_type_t <T_y ,T_dof ,T_loc ,T_covar >multi_student_t_cholesky_lpdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_covar &L );

// ./stan/math/prim/prob/multi_student_t_cholesky_lpdf.hpp:358
template <typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_covar>
inline return_type_t <T_y ,T_dof ,T_loc ,T_covar >multi_student_t_cholesky_lpdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_covar &L );

// ./stan/math/prim/prob/multi_student_t_cholesky_rng.hpp:40
template <typename T_loc,
    class RNG>
inline typename StdVectorBuilder <true ,Eigen ::VectorXd ,T_loc >::type multi_student_t_cholesky_rng (double nu , const T_loc &mu , const Eigen ::MatrixXd &L , RNG &rng );

// ./stan/math/prim/prob/multi_student_t_lpdf.hpp:41
template <bool propto,
    typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_loc ,T_scale >multi_student_t_lpdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_scale &Sigma );

// ./stan/math/prim/prob/multi_student_t_lpdf.hpp:134
template <typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_dof ,T_loc ,T_scale >multi_student_t_lpdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_scale &Sigma );

// ./stan/math/prim/prob/multi_student_t_rng.hpp:40
template <typename T_loc,
    class RNG>
inline typename StdVectorBuilder <true ,Eigen ::VectorXd ,T_loc >::type multi_student_t_rng (double nu , const T_loc &mu , const Eigen ::Matrix <double , Eigen ::Dynamic , Eigen ::Dynamic >&S , RNG &rng );

// ./stan/math/prim/prob/multinomial_logit_lpmf.hpp:23
template <bool propto,
    typename T_beta,
    typename T_prob>
return_type_t <T_prob >multinomial_logit_lpmf (const std ::vector <int >&ns , const T_beta &beta );

// ./stan/math/prim/prob/multinomial_logit_lpmf.hpp:54
template <typename T_beta,
    require_eigen_col_vector_t<T_beta>*>
return_type_t <T_beta >multinomial_logit_lpmf (const std ::vector <int >&ns , const T_beta &beta );

// ./stan/math/prim/prob/multinomial_logit_rng.hpp:27
template <class RNG,
    typename T_beta,
    require_eigen_col_vector_t<T_beta>*>
inline std ::vector <int >multinomial_logit_rng (const T_beta &beta , int N , RNG &rng );

// ./stan/math/prim/prob/multinomial_lpmf.hpp:15
template <bool propto,
    typename T_prob,
    require_eigen_col_vector_t<T_prob>*>
return_type_t <T_prob >multinomial_lpmf (const std ::vector <int >&ns , const T_prob &theta );

// ./stan/math/prim/prob/multinomial_lpmf.hpp:44
template <typename T_prob>
return_type_t <T_prob >multinomial_lpmf (const std ::vector <int >&ns , const T_prob &theta );

// ./stan/math/prim/prob/neg_binomial_lccdf.hpp:26
template <typename T_n,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_shape ,T_inv_scale >neg_binomial_lccdf (const T_n &n , const T_shape &alpha , const T_inv_scale &beta_param );

// ./stan/math/prim/prob/neg_binomial_ccdf_log.hpp:13
template <typename T_n,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_shape ,T_inv_scale >neg_binomial_ccdf_log (const T_n &n , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/neg_binomial_2_lccdf.hpp:16
template <typename T_n,
    typename T_location,
    typename T_precision>
return_type_t <T_location ,T_precision >neg_binomial_2_lccdf (const T_n &n , const T_location &mu , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_ccdf_log.hpp:13
template <typename T_n,
    typename T_location,
    typename T_precision>
return_type_t <T_location ,T_precision >neg_binomial_2_ccdf_log (const T_n &n , const T_location &mu , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_cdf.hpp:23
template <typename T_n,
    typename T_location,
    typename T_precision>
return_type_t <T_location ,T_precision >neg_binomial_2_cdf (const T_n &n , const T_location &mu , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_lcdf.hpp:17
template <typename T_n,
    typename T_location,
    typename T_precision>
return_type_t <T_location ,T_precision >neg_binomial_2_lcdf (const T_n &n , const T_location &mu , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_cdf_log.hpp:13
template <typename T_n,
    typename T_location,
    typename T_precision>
return_type_t <T_location ,T_precision >neg_binomial_2_cdf_log (const T_n &n , const T_location &mu , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_log_glm_lpmf.hpp:64
template <bool propto,
    typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta,
    typename T_precision,
    require_matrix_t<T_x>*>
return_type_t <T_x ,T_alpha ,T_beta ,T_precision >neg_binomial_2_log_glm_lpmf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_log_glm_lpmf.hpp:250
template <typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta,
    typename T_precision>
inline return_type_t <T_x ,T_alpha ,T_beta ,T_precision >neg_binomial_2_log_glm_lpmf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_log_lpmf.hpp:24
template <bool propto,
    typename T_n,
    typename T_log_location,
    typename T_precision,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_log_location,
    T_precision>*>
return_type_t <T_log_location ,T_precision >neg_binomial_2_log_lpmf (const T_n &n , const T_log_location &eta , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_log_lpmf.hpp:134
template <typename T_n,
    typename T_log_location,
    typename T_precision>
inline return_type_t <T_log_location ,T_precision >neg_binomial_2_log_lpmf (const T_n &n , const T_log_location &eta , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_log_rng.hpp:35
template <typename T_loc,
    typename T_inv,
    class RNG>
inline typename VectorBuilder <true ,int ,T_loc ,T_inv >::type neg_binomial_2_log_rng (const T_loc &eta , const T_inv &phi , RNG &rng );

// ./stan/math/prim/prob/neg_binomial_2_lpmf.hpp:22
template <bool propto,
    typename T_n,
    typename T_location,
    typename T_precision,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_location,
    T_precision>*>
return_type_t <T_location ,T_precision >neg_binomial_2_lpmf (const T_n &n , const T_location &mu , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_lpmf.hpp:119
template <typename T_n,
    typename T_location,
    typename T_precision>
inline return_type_t <T_location ,T_precision >neg_binomial_2_lpmf (const T_n &n , const T_location &mu , const T_precision &phi );

// ./stan/math/prim/prob/neg_binomial_2_rng.hpp:34
template <typename T_loc,
    typename T_prec,
    class RNG>
inline typename VectorBuilder <true ,int ,T_loc ,T_prec >::type neg_binomial_2_rng (const T_loc &mu , const T_prec &phi , RNG &rng );

// ./stan/math/prim/prob/neg_binomial_cdf.hpp:24
template <typename T_n,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_shape ,T_inv_scale >neg_binomial_cdf (const T_n &n , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/neg_binomial_lcdf.hpp:26
template <typename T_n,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_shape ,T_inv_scale >neg_binomial_lcdf (const T_n &n , const T_shape &alpha , const T_inv_scale &beta_param );

// ./stan/math/prim/prob/neg_binomial_cdf_log.hpp:13
template <typename T_n,
    typename T_shape,
    typename T_inv_scale>
return_type_t <T_shape ,T_inv_scale >neg_binomial_cdf_log (const T_n &n , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/neg_binomial_lpmf.hpp:30
template <bool propto,
    typename T_n,
    typename T_shape,
    typename T_inv_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_shape,
    T_inv_scale>*>
return_type_t <T_shape ,T_inv_scale >neg_binomial_lpmf (const T_n &n , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/neg_binomial_lpmf.hpp:122
template <typename T_n,
    typename T_shape,
    typename T_inv_scale>
inline return_type_t <T_shape ,T_inv_scale >neg_binomial_lpmf (const T_n &n , const T_shape &alpha , const T_inv_scale &beta );

// ./stan/math/prim/prob/neg_binomial_rng.hpp:34
template <typename T_shape,
    typename T_inv,
    class RNG>
inline typename VectorBuilder <true ,int ,T_shape ,T_inv >::type neg_binomial_rng (const T_shape &alpha , const T_inv &beta , RNG &rng );

// ./stan/math/prim/prob/normal_lccdf.hpp:21
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
inline return_type_t <T_y ,T_loc ,T_scale >normal_lccdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/normal_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_loc ,T_scale >normal_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/normal_cdf.hpp:35
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
inline return_type_t <T_y ,T_loc ,T_scale >normal_cdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/normal_lcdf.hpp:24
template <typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
inline return_type_t <T_y ,T_loc ,T_scale >normal_lcdf (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/normal_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_loc ,T_scale >normal_cdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/normal_id_glm_lpdf.hpp:54
template <bool propto,
    typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta,
    typename T_scale,
    require_matrix_t<T_x>*>
return_type_t <T_y ,T_x ,T_alpha ,T_beta ,T_scale >normal_id_glm_lpdf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta , const T_scale &sigma );

// ./stan/math/prim/prob/normal_id_glm_lpdf.hpp:218
template <typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta,
    typename T_scale>
inline return_type_t <T_y ,T_x ,T_alpha ,T_beta ,T_scale >normal_id_glm_lpdf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta , const T_scale &sigma );

// ./stan/math/prim/prob/normal_lpdf.hpp:41
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale>*>
inline return_type_t <T_y ,T_loc ,T_scale >normal_lpdf (T_y &&y , T_loc &&mu , T_scale &&sigma );

// ./stan/math/prim/prob/normal_lpdf.hpp:107
template <typename T_y,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_loc ,T_scale >normal_lpdf (T_y &&y , T_loc &&mu , T_scale &&sigma );

// ./stan/math/prim/prob/normal_sufficient_lpdf.hpp:52
template <bool propto,
    typename T_y,
    typename T_s,
    typename T_n,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_s ,T_loc ,T_scale >normal_sufficient_lpdf (const T_y &y_bar , const T_s &s_squared , const T_n &n_obs , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/normal_sufficient_lpdf.hpp:155
template <typename T_y,
    typename T_s,
    typename T_n,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_s ,T_loc ,T_scale >normal_sufficient_lpdf (const T_y &y_bar , const T_s &s_squared , const T_n &n_obs , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/ordered_logistic_glm_lpmf.hpp:46
template <bool propto,
    typename T_y,
    typename T_x,
    typename T_beta,
    typename T_cuts,
    require_matrix_t<T_x>*>
return_type_t <T_x ,T_beta ,T_cuts >ordered_logistic_glm_lpmf (const T_y &y , const T_x &x , const T_beta &beta , const T_cuts &cuts );

// ./stan/math/prim/prob/ordered_logistic_glm_lpmf.hpp:212
template <typename T_y,
    typename T_x,
    typename T_beta,
    typename T_cuts>
return_type_t <T_x ,T_beta ,T_cuts >ordered_logistic_glm_lpmf (const T_y &y , const T_x &x , const T_beta &beta , const T_cuts &cuts );

// ./stan/math/prim/fun/is_integer.hpp:17
template <typename T>
inline bool is_integer (T x );

// ./stan/math/prim/prob/ordered_logistic_lpmf.hpp:72
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_cut,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_cut>*>
return_type_t <T_loc ,T_cut >ordered_logistic_lpmf (const T_y &y , const T_loc &lambda , const T_cut &c );

// ./stan/math/prim/prob/ordered_logistic_lpmf.hpp:206
template <typename T_y,
    typename T_loc,
    typename T_cut>
return_type_t <T_loc ,T_cut >ordered_logistic_lpmf (const T_y &y , const T_loc &lambda , const T_cut &c );

// ./stan/math/prim/prob/ordered_logistic_rng.hpp:13
template <class RNG>
inline int ordered_logistic_rng (double eta , const Eigen ::Matrix <double , Eigen ::Dynamic , 1>&c , RNG &rng );

// ./stan/math/prim/prob/std_normal_lcdf.hpp:24
template <typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>*>
inline return_type_t <T_y >std_normal_lcdf (const T_y &y );

// ./stan/math/prim/prob/ordered_probit_lpmf.hpp:47
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_cut>
return_type_t <T_loc ,T_cut >ordered_probit_lpmf (const T_y &y , const T_loc &lambda , const T_cut &c );

// ./stan/math/prim/prob/ordered_probit_lpmf.hpp:103
template <typename T_y,
    typename T_loc,
    typename T_cut>
return_type_t <T_loc ,T_cut >ordered_probit_lpmf (const T_y &y , const T_loc &lambda , const T_cut &c );

// ./stan/math/prim/prob/ordered_probit_rng.hpp:12
template <class RNG>
inline int ordered_probit_rng (double eta , const Eigen ::VectorXd &c , RNG &rng );

// ./stan/math/prim/prob/pareto_lccdf.hpp:24
template <typename T_y,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_scale ,T_shape >pareto_lccdf (const T_y &y , const T_scale &y_min , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_ccdf_log.hpp:13
template <typename T_y,
    typename T_scale,
    typename T_shape>
return_type_t <T_y ,T_scale ,T_shape >pareto_ccdf_log (const T_y &y , const T_scale &y_min , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_cdf.hpp:20
template <typename T_y,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_scale ,T_shape >pareto_cdf (const T_y &y , const T_scale &y_min , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_lcdf.hpp:24
template <typename T_y,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_scale ,T_shape >pareto_lcdf (const T_y &y , const T_scale &y_min , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_cdf_log.hpp:13
template <typename T_y,
    typename T_scale,
    typename T_shape>
return_type_t <T_y ,T_scale ,T_shape >pareto_cdf_log (const T_y &y , const T_scale &y_min , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_lpdf.hpp:24
template <bool propto,
    typename T_y,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_scale ,T_shape >pareto_lpdf (const T_y &y , const T_scale &y_min , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_lpdf.hpp:94
template <typename T_y,
    typename T_scale,
    typename T_shape>
inline return_type_t <T_y ,T_scale ,T_shape >pareto_lpdf (const T_y &y , const T_scale &y_min , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_rng.hpp:34
template <typename T_shape,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_shape ,T_scale >::type pareto_rng (const T_scale &y_min , const T_shape &alpha , RNG &rng );

// ./stan/math/prim/prob/pareto_type_2_lccdf.hpp:21
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >pareto_type_2_lccdf (const T_y &y , const T_loc &mu , const T_scale &lambda , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_type_2_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >pareto_type_2_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &lambda , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_type_2_cdf.hpp:21
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >pareto_type_2_cdf (const T_y &y , const T_loc &mu , const T_scale &lambda , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_type_2_lcdf.hpp:21
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >pareto_type_2_lcdf (const T_y &y , const T_loc &mu , const T_scale &lambda , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_type_2_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >pareto_type_2_cdf_log (const T_y &y , const T_loc &mu , const T_scale &lambda , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_type_2_lpdf.hpp:23
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >pareto_type_2_lpdf (const T_y &y , const T_loc &mu , const T_scale &lambda , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_type_2_lpdf.hpp:111
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
inline return_type_t <T_y ,T_loc ,T_scale ,T_shape >pareto_type_2_lpdf (const T_y &y , const T_loc &mu , const T_scale &lambda , const T_shape &alpha );

// ./stan/math/prim/prob/pareto_type_2_rng.hpp:40
template <typename T_loc,
    typename T_scale,
    typename T_shape,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale ,T_shape >::type pareto_type_2_rng (const T_loc &mu , const T_scale &lambda , const T_shape &alpha , RNG &rng );

// ./stan/math/prim/prob/poisson_binomial_lccdf.hpp:36
template <bool propto,
    typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_lccdf (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_lccdf.hpp:76
template <typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_lccdf (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_ccdf_log.hpp:13
template <typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_ccdf_log (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_cdf.hpp:35
template <bool propto,
    typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_cdf (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_cdf.hpp:65
template <typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_cdf (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_lcdf.hpp:35
template <bool propto,
    typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_lcdf (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_lcdf.hpp:66
template <typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_lcdf (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_cdf_log.hpp:13
template <typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_cdf_log (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_lpmf.hpp:27
template <bool propto,
    typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_lpmf (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_lpmf.hpp:59
template <typename T_y,
    typename T_theta>
return_type_t <T_theta >poisson_binomial_lpmf (const T_y &y , const T_theta &theta );

// ./stan/math/prim/prob/poisson_binomial_rng.hpp:24
template <typename T_theta,
    typename RNG,
    require_eigen_vt<std::is_arithmetic,
    T_theta>*>
inline int poisson_binomial_rng (const T_theta &theta , RNG &rng );

// ./stan/math/prim/prob/poisson_lccdf.hpp:27
template <typename T_n,
    typename T_rate>
return_type_t <T_rate >poisson_lccdf (const T_n &n , const T_rate &lambda );

// ./stan/math/prim/prob/poisson_ccdf_log.hpp:13
template <typename T_n,
    typename T_rate>
return_type_t <T_rate >poisson_ccdf_log (const T_n &n , const T_rate &lambda );

// ./stan/math/prim/prob/poisson_cdf.hpp:26
template <typename T_n,
    typename T_rate>
return_type_t <T_rate >poisson_cdf (const T_n &n , const T_rate &lambda );

// ./stan/math/prim/prob/poisson_lcdf.hpp:27
template <typename T_n,
    typename T_rate>
return_type_t <T_rate >poisson_lcdf (const T_n &n , const T_rate &lambda );

// ./stan/math/prim/prob/poisson_cdf_log.hpp:13
template <typename T_n,
    typename T_rate>
return_type_t <T_rate >poisson_cdf_log (const T_n &n , const T_rate &lambda );

// ./stan/math/prim/prob/poisson_log_glm_lpmf.hpp:51
template <bool propto,
    typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta,
    require_matrix_t<T_x>*>
return_type_t <T_x ,T_alpha ,T_beta >poisson_log_glm_lpmf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta );

// ./stan/math/prim/prob/poisson_log_glm_lpmf.hpp:165
template <typename T_y,
    typename T_x,
    typename T_alpha,
    typename T_beta>
inline return_type_t <T_x ,T_alpha ,T_beta >poisson_log_glm_lpmf (const T_y &y , const T_x &x , const T_alpha &alpha , const T_beta &beta );

// ./stan/math/prim/prob/poisson_log_lpmf.hpp:26
template <bool propto,
    typename T_n,
    typename T_log_rate,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_log_rate>*>
return_type_t <T_log_rate >poisson_log_lpmf (const T_n &n , const T_log_rate &alpha );

// ./stan/math/prim/prob/poisson_log_lpmf.hpp:88
template <typename T_n,
    typename T_log_rate>
inline return_type_t <T_log_rate >poisson_log_lpmf (const T_n &n , const T_log_rate &alpha );

// ./stan/math/prim/prob/poisson_log_rng.hpp:30
template <typename T_rate,
    class RNG>
inline typename VectorBuilder <true ,int ,T_rate >::type poisson_log_rng (const T_rate &alpha , RNG &rng );

// ./stan/math/prim/prob/poisson_lpmf.hpp:26
template <bool propto,
    typename T_n,
    typename T_rate,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_n,
    T_rate>*>
return_type_t <T_rate >poisson_lpmf (const T_n &n , const T_rate &lambda );

// ./stan/math/prim/prob/poisson_lpmf.hpp:84
template <typename T_n,
    typename T_rate>
inline return_type_t <T_rate >poisson_lpmf (const T_n &n , const T_rate &lambda );

// ./stan/math/prim/prob/poisson_rng.hpp:29
template <typename T_rate,
    class RNG>
inline typename VectorBuilder <true ,int ,T_rate >::type poisson_rng (const T_rate &lambda , RNG &rng );

// ./stan/math/prim/prob/rayleigh_lccdf.hpp:20
template <typename T_y,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale>*>
return_type_t <T_y ,T_scale >rayleigh_lccdf (const T_y &y , const T_scale &sigma );

// ./stan/math/prim/prob/rayleigh_ccdf_log.hpp:13
template <typename T_y,
    typename T_scale>
return_type_t <T_y ,T_scale >rayleigh_ccdf_log (const T_y &y , const T_scale &sigma );

// ./stan/math/prim/prob/rayleigh_cdf.hpp:22
template <typename T_y,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale>*>
return_type_t <T_y ,T_scale >rayleigh_cdf (const T_y &y , const T_scale &sigma );

// ./stan/math/prim/prob/rayleigh_lcdf.hpp:22
template <typename T_y,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale>*>
return_type_t <T_y ,T_scale >rayleigh_lcdf (const T_y &y , const T_scale &sigma );

// ./stan/math/prim/prob/rayleigh_cdf_log.hpp:13
template <typename T_y,
    typename T_scale>
return_type_t <T_y ,T_scale >rayleigh_cdf_log (const T_y &y , const T_scale &sigma );

// ./stan/math/prim/prob/rayleigh_lpdf.hpp:22
template <bool propto,
    typename T_y,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_scale>*>
return_type_t <T_y ,T_scale >rayleigh_lpdf (const T_y &y , const T_scale &sigma );

// ./stan/math/prim/prob/rayleigh_lpdf.hpp:81
template <typename T_y,
    typename T_scale>
inline return_type_t <T_y ,T_scale >rayleigh_lpdf (const T_y &y , const T_scale &sigma );

// ./stan/math/prim/prob/rayleigh_rng.hpp:29
template <typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_scale >::type rayleigh_rng (const T_scale &sigma , RNG &rng );

// ./stan/math/prim/prob/scaled_inv_chi_square_lccdf.hpp:24
template <typename T_y,
    typename T_dof,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_scale >scaled_inv_chi_square_lccdf (const T_y &y , const T_dof &nu , const T_scale &s );

// ./stan/math/prim/prob/scaled_inv_chi_square_ccdf_log.hpp:13
template <typename T_y,
    typename T_dof,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_scale >scaled_inv_chi_square_ccdf_log (const T_y &y , const T_dof &nu , const T_scale &s );

// ./stan/math/prim/prob/scaled_inv_chi_square_cdf.hpp:37
template <typename T_y,
    typename T_dof,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_scale >scaled_inv_chi_square_cdf (const T_y &y , const T_dof &nu , const T_scale &s );

// ./stan/math/prim/prob/scaled_inv_chi_square_lcdf.hpp:24
template <typename T_y,
    typename T_dof,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_scale >scaled_inv_chi_square_lcdf (const T_y &y , const T_dof &nu , const T_scale &s );

// ./stan/math/prim/prob/scaled_inv_chi_square_cdf_log.hpp:13
template <typename T_y,
    typename T_dof,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_scale >scaled_inv_chi_square_cdf_log (const T_y &y , const T_dof &nu , const T_scale &s );

// ./stan/math/prim/prob/scaled_inv_chi_square_lpdf.hpp:44
template <bool propto,
    typename T_y,
    typename T_dof,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_dof,
    T_scale>*>
return_type_t <T_y ,T_dof ,T_scale >scaled_inv_chi_square_lpdf (const T_y &y , const T_dof &nu , const T_scale &s );

// ./stan/math/prim/prob/scaled_inv_chi_square_lpdf.hpp:174
template <typename T_y,
    typename T_dof,
    typename T_scale>
inline return_type_t <T_y ,T_dof ,T_scale >scaled_inv_chi_square_lpdf (const T_y &y , const T_dof &nu , const T_scale &s );

// ./stan/math/prim/prob/scaled_inv_chi_square_rng.hpp:34
template <typename T_deg,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_deg ,T_scale >::type scaled_inv_chi_square_rng (const T_deg &nu , const T_scale &s , RNG &rng );

// ./stan/math/prim/prob/skew_double_exponential_lccdf.hpp:37
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_skewness,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_skewness>*>
return_type_t <T_y ,T_loc ,T_scale ,T_skewness >skew_double_exponential_lccdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_skewness &tau );

// ./stan/math/prim/prob/skew_double_exponential_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_skewness>
inline return_type_t <T_y ,T_loc ,T_scale ,T_skewness >skew_double_exponential_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_skewness &tau );

// ./stan/math/prim/prob/skew_double_exponential_cdf.hpp:38
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_skewness,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_skewness>*>
return_type_t <T_y ,T_loc ,T_scale ,T_skewness >skew_double_exponential_cdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_skewness &tau );

// ./stan/math/prim/prob/skew_double_exponential_lcdf.hpp:38
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_skewness,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_skewness>*>
return_type_t <T_y ,T_loc ,T_scale ,T_skewness >skew_double_exponential_lcdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_skewness &tau );

// ./stan/math/prim/prob/skew_double_exponential_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_skewness>
inline return_type_t <T_y ,T_loc ,T_scale ,T_skewness >skew_double_exponential_cdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_skewness &tau );

// ./stan/math/prim/prob/skew_double_exponential_lpdf.hpp:39
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_skewness,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_skewness>*>
return_type_t <T_y ,T_loc ,T_scale ,T_skewness >skew_double_exponential_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_skewness &tau );

// ./stan/math/prim/prob/skew_double_exponential_lpdf.hpp:130
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_skewness>
return_type_t <T_y ,T_loc ,T_scale ,T_skewness >skew_double_exponential_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_skewness &tau );

// ./stan/math/prim/prob/skew_double_exponential_rng.hpp:36
template <typename T_loc,
    typename T_scale,
    typename T_skewness,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale ,T_skewness >::type skew_double_exponential_rng (const T_loc &mu , const T_scale &sigma , const T_skewness &tau , RNG &rng );

// ./stan/math/prim/prob/skew_normal_lccdf.hpp:25
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >skew_normal_lccdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &alpha );

// ./stan/math/prim/prob/skew_normal_ccdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >skew_normal_ccdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &alpha );

// ./stan/math/prim/prob/skew_normal_cdf.hpp:27
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >skew_normal_cdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &alpha );

// ./stan/math/prim/prob/skew_normal_lcdf.hpp:27
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >skew_normal_lcdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &alpha );

// ./stan/math/prim/prob/skew_normal_cdf_log.hpp:13
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >skew_normal_cdf_log (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &alpha );

// ./stan/math/prim/prob/skew_normal_lpdf.hpp:25
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_loc,
    T_scale,
    T_shape>*>
return_type_t <T_y ,T_loc ,T_scale ,T_shape >skew_normal_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &alpha );

// ./stan/math/prim/prob/skew_normal_lpdf.hpp:118
template <typename T_y,
    typename T_loc,
    typename T_scale,
    typename T_shape>
inline return_type_t <T_y ,T_loc ,T_scale ,T_shape >skew_normal_lpdf (const T_y &y , const T_loc &mu , const T_scale &sigma , const T_shape &alpha );

// ./stan/math/prim/prob/skew_normal_rng.hpp:38
template <typename T_loc,
    typename T_scale,
    typename T_shape,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_scale ,T_shape >::type skew_normal_rng (const T_loc &mu , const T_scale &sigma , const T_shape &alpha , RNG &rng );

// ./stan/math/prim/prob/std_normal_lccdf.hpp:21
template <typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>*>
inline return_type_t <T_y >std_normal_lccdf (const T_y &y );

// ./stan/math/prim/prob/std_normal_ccdf_log.hpp:13
template <typename T_y>
inline return_type_t <T_y >std_normal_ccdf_log (const T_y &y );

// ./stan/math/prim/prob/std_normal_cdf.hpp:30
template <typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>*>
inline return_type_t <T_y >std_normal_cdf (const T_y &y );

// ./stan/math/prim/prob/std_normal_cdf_log.hpp:13
template <typename T_y>
inline return_type_t <T_y >std_normal_cdf_log (const T_y &y );

// ./stan/math/prim/prob/std_normal_lpdf.hpp:29
template <bool propto,
    typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>*>
return_type_t <T_y >std_normal_lpdf (const T_y &y );

// ./stan/math/prim/prob/std_normal_lpdf.hpp:67
template <typename T_y>
inline return_type_t <T_y >std_normal_lpdf (const T_y &y );

// ./stan/math/prim/prob/std_normal_rng.hpp:19
template <class RNG>
inline double std_normal_rng (RNG &rng );

// ./stan/math/prim/prob/student_t_lccdf.hpp:23
template <typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_loc ,T_scale >student_t_lccdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/student_t_ccdf_log.hpp:13
template <typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_loc ,T_scale >student_t_ccdf_log (const T_y &y , const T_dof &nu , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/student_t_cdf.hpp:22
template <typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_loc ,T_scale >student_t_cdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/student_t_lcdf.hpp:23
template <typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_loc ,T_scale >student_t_lcdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/student_t_cdf_log.hpp:13
template <typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_dof ,T_loc ,T_scale >student_t_cdf_log (const T_y &y , const T_dof &nu , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/student_t_lpdf.hpp:55
template <bool propto,
    typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_dof,
    T_loc,
    T_scale>*>
return_type_t <T_y ,T_dof ,T_loc ,T_scale >student_t_lpdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/student_t_lpdf.hpp:153
template <typename T_y,
    typename T_dof,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_dof ,T_loc ,T_scale >student_t_lpdf (const T_y &y , const T_dof &nu , const T_loc &mu , const T_scale &sigma );

// ./stan/math/prim/prob/student_t_rng.hpp:36
template <typename T_deg,
    typename T_loc,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_deg ,T_loc ,T_scale >::type student_t_rng (const T_deg &nu , const T_loc &mu , const T_scale &sigma , RNG &rng );

// ./stan/math/prim/prob/uniform_lccdf.hpp:22
template <typename T_y,
    typename T_low,
    typename T_high,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_low,
    T_high>*>
return_type_t <T_y ,T_low ,T_high >uniform_lccdf (const T_y &y , const T_low &alpha , const T_high &beta );

// ./stan/math/prim/prob/uniform_ccdf_log.hpp:13
template <typename T_y,
    typename T_low,
    typename T_high>
return_type_t <T_y ,T_low ,T_high >uniform_ccdf_log (const T_y &y , const T_low &alpha , const T_high &beta );

// ./stan/math/prim/prob/uniform_cdf.hpp:20
template <typename T_y,
    typename T_low,
    typename T_high,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_low,
    T_high>*>
return_type_t <T_y ,T_low ,T_high >uniform_cdf (const T_y &y , const T_low &alpha , const T_high &beta );

// ./stan/math/prim/prob/uniform_lcdf.hpp:22
template <typename T_y,
    typename T_low,
    typename T_high,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_low,
    T_high>*>
return_type_t <T_y ,T_low ,T_high >uniform_lcdf (const T_y &y , const T_low &alpha , const T_high &beta );

// ./stan/math/prim/prob/uniform_cdf_log.hpp:13
template <typename T_y,
    typename T_low,
    typename T_high>
return_type_t <T_y ,T_low ,T_high >uniform_cdf_log (const T_y &y , const T_low &alpha , const T_high &beta );

// ./stan/math/prim/prob/uniform_lpdf.hpp:44
template <bool propto,
    typename T_y,
    typename T_low,
    typename T_high,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_low,
    T_high>*>
return_type_t <T_y ,T_low ,T_high >uniform_lpdf (const T_y &y , const T_low &alpha , const T_high &beta );

// ./stan/math/prim/prob/uniform_lpdf.hpp:115
template <typename T_y,
    typename T_low,
    typename T_high>
inline return_type_t <T_y ,T_low ,T_high >uniform_lpdf (const T_y &y , const T_low &alpha , const T_high &beta );

// ./stan/math/prim/prob/von_mises_lpdf.hpp:25
template <bool propto,
    typename T_y,
    typename T_loc,
    typename T_scale>
return_type_t <T_y ,T_loc ,T_scale >von_mises_lpdf (T_y const &y , T_loc const &mu , T_scale const &kappa );

// ./stan/math/prim/prob/von_mises_lpdf.hpp:91
template <typename T_y,
    typename T_loc,
    typename T_scale>
inline return_type_t <T_y ,T_loc ,T_scale >von_mises_lpdf (T_y const &y , T_loc const &mu , T_scale const &kappa );

// ./stan/math/prim/prob/von_mises_rng.hpp:45
template <typename T_loc,
    typename T_conc,
    class RNG>
inline typename VectorBuilder <true ,double ,T_loc ,T_conc >::type von_mises_rng (const T_loc &mu , const T_conc &kappa , RNG &rng );

// ./stan/math/prim/prob/von_mises_cdf.hpp:117
template <typename T_x,
    typename T_mu,
    typename T_k>
inline return_type_t <T_x ,T_mu ,T_k >von_mises_cdf (const T_x &x , const T_mu &mu , const T_k &k );

// ./stan/math/prim/prob/von_mises_lccdf.hpp:31
template <typename T_x,
    typename T_mu,
    typename T_k>
inline return_type_t <T_x ,T_mu ,T_k >von_mises_lccdf (const T_x &x , const T_mu &mu , const T_k &k );

// ./stan/math/prim/prob/von_mises_ccdf_log.hpp:13
template <typename T_x,
    typename T_mu,
    typename T_k>
inline return_type_t <T_x ,T_mu ,T_k >von_mises_ccdf_log (const T_x &x , const T_mu &mu , const T_k &k );

// ./stan/math/prim/prob/von_mises_lcdf.hpp:31
template <typename T_x,
    typename T_mu,
    typename T_k>
inline return_type_t <T_x ,T_mu ,T_k >von_mises_lcdf (const T_x &x , const T_mu &mu , const T_k &k );

// ./stan/math/prim/prob/von_mises_cdf_log.hpp:13
template <typename T_x,
    typename T_mu,
    typename T_k>
inline return_type_t <T_x ,T_mu ,T_k >von_mises_cdf_log (const T_x &x , const T_mu &mu , const T_k &k );

// ./stan/math/prim/prob/weibull_lccdf.hpp:34
template <typename T_y,
    typename T_shape,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_scale>*>
return_type_t <T_y ,T_shape ,T_scale >weibull_lccdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/weibull_ccdf_log.hpp:13
template <typename T_y,
    typename T_shape,
    typename T_scale>
return_type_t <T_y ,T_shape ,T_scale >weibull_ccdf_log (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/weibull_cdf.hpp:36
template <typename T_y,
    typename T_shape,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_scale>*>
return_type_t <T_y ,T_shape ,T_scale >weibull_cdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/weibull_lcdf.hpp:36
template <typename T_y,
    typename T_shape,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_scale>*>
return_type_t <T_y ,T_shape ,T_scale >weibull_lcdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/weibull_cdf_log.hpp:13
template <typename T_y,
    typename T_shape,
    typename T_scale>
return_type_t <T_y ,T_shape ,T_scale >weibull_cdf_log (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/weibull_lpdf.hpp:37
template <bool propto,
    typename T_y,
    typename T_shape,
    typename T_scale,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y,
    T_shape,
    T_scale>*>
return_type_t <T_y ,T_shape ,T_scale >weibull_lpdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/weibull_lpdf.hpp:115
template <typename T_y,
    typename T_shape,
    typename T_scale>
inline return_type_t <T_y ,T_shape ,T_scale >weibull_lpdf (const T_y &y , const T_shape &alpha , const T_scale &sigma );

// ./stan/math/prim/prob/weibull_rng.hpp:32
template <typename T_shape,
    typename T_scale,
    class RNG>
inline typename VectorBuilder <true ,double ,T_shape ,T_scale >::type weibull_rng (const T_shape &alpha , const T_scale &sigma , RNG &rng );

// ./stan/math/prim/prob/wiener5_lpdf.hpp:672
template <bool propto>
inline auto wiener_lpdf (const T_y &y , const T_a &a , const T_t0 &t0 , const T_w &w , const T_v &v , const T_sv &sv , const double &precision_derivatives =1e-4);

// ./stan/math/prim/prob/wiener_lpdf.hpp:75
template <bool propto,
    typename T_y,
    typename T_alpha,
    typename T_tau,
    typename T_beta,
    typename T_delta>
return_type_t <T_y ,T_alpha ,T_tau ,T_beta ,T_delta >wiener_lpdf (const T_y &y , const T_alpha &alpha , const T_tau &tau , const T_beta &beta , const T_delta &delta );

// ./stan/math/prim/prob/wiener_lpdf.hpp:209
template <typename T_y,
    typename T_alpha,
    typename T_tau,
    typename T_beta,
    typename T_delta>
inline return_type_t <T_y ,T_alpha ,T_tau ,T_beta ,T_delta >wiener_lpdf (const T_y &y , const T_alpha &alpha , const T_tau &tau , const T_beta &beta , const T_delta &delta );

// ./stan/math/prim/prob/wiener_full_lpdf.hpp:317
template <bool propto>
inline auto wiener_lpdf (const T_y &y , const T_a &a , const T_t0 &t0 , const T_w &w , const T_v &v , const T_sv &sv , const T_sw &sw , const T_st0 &st0 , const double &precision_derivatives =1e-4);

// ./stan/math/prim/prob/wishart_cholesky_lpdf.hpp:39
template <bool propto,
    typename T_y,
    typename T_dof,
    typename T_scale,
    require_stan_scalar_t<T_dof>*>
return_type_t <T_y ,T_dof ,T_scale >wishart_cholesky_lpdf (const T_y &L_Y , const T_dof &nu , const T_scale &L_S );

// ./stan/math/prim/prob/wishart_cholesky_lpdf.hpp:90
template <typename T_y,
    typename T_dof,
    typename T_scale>
inline return_type_t <T_y ,T_dof ,T_scale >wishart_cholesky_lpdf (const T_y &LW , const T_dof &nu , const T_scale &L_S );

// ./stan/math/prim/prob/wishart_cholesky_rng.hpp:30
template <class RNG>
inline Eigen ::MatrixXd wishart_cholesky_rng (double nu , const Eigen ::MatrixXd &L_S , RNG &rng );

// ./stan/math/prim/prob/wishart_lpdf.hpp:46
template <bool propto,
    typename T_y,
    typename T_dof,
    typename T_scale,
    require_stan_scalar_t<T_dof>*>
return_type_t <T_y ,T_dof ,T_scale >wishart_lpdf (const T_y &W , const T_dof &nu , const T_scale &S );

// ./stan/math/prim/prob/wishart_lpdf.hpp:101
template <typename T_y,
    typename T_dof,
    typename T_scale>
inline return_type_t <T_y ,T_dof ,T_scale >wishart_lpdf (const T_y &W , const T_dof &nu , const T_scale &S );

// Functions Parsed: 2357

} // namespace math
} // namespace stan
#endif // STAN_MATH_FORWARD_DECL_HPP

