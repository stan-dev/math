
#ifndef STAN_MATH_FORWARD_DECL_HPP
#define STAN_MATH_FORWARD_DECL_HPP
#include <stan/math/manual_forward_decl.hpp>
namespace stan {
namespace math {
namespace internal {

template <typename F>
struct reverse_pass_callback_vari;

template <typename T, typename F>
struct callback_vari;

class multiply_vv_vari;

class multiply_vd_vari;

struct hash_profile_key;

struct equal_profile_key;

struct nonexisting_adjoint;

template <typename Result_, typename WMat_, typename B_>
struct csr_adjoint;

class fdim_vv_vari;

class fdim_vd_vari;

class fdim_dv_vari;

class fmod_vv_vari;

class fmod_vd_vari;

class fmod_dv_vari;

class gamma_p_vv_vari;

class gamma_p_vd_vari;

class gamma_p_dv_vari;

class gamma_q_vv_vari;

class gamma_q_vd_vari;

class gamma_q_dv_vari;

class inc_beta_vvv_vari;

class lbeta_vv_vari;

class lbeta_vd_vari;

class lbeta_dv_vari;

class lmgamma_dv_vari;

class lmultiply_vv_vari;

class lmultiply_vd_vari;

class lmultiply_dv_vari;

class log_diff_exp_vv_vari;

class log_diff_exp_vd_vari;

class log_diff_exp_dv_vari;

class log_falling_factorial_vv_vari;

class log_falling_factorial_vd_vari;

class log_falling_factorial_dv_vari;

class log_inv_logit_diff_vv_vari;

class log_inv_logit_diff_vd_vari;

class log_inv_logit_diff_dv_vari;







class log_rising_factorial_vv_vari;

class log_rising_factorial_vd_vari;

class log_rising_factorial_dv_vari;

class log_softmax_elt_vari;

class log_sum_exp_vv_vari;

class log_sum_exp_vd_vari;

template <int R1, int C1, int R2, int C2>
class mdivide_left_spd_alloc;

template <int R1, int C1, int R2, int C2>
class mdivide_left_spd_vv_vari;

template <int R1, int C1, int R2, int C2>
class mdivide_left_spd_dv_vari;

template <int R1, int C1, int R2, int C2>
class mdivide_left_spd_vd_vari;

template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
class mdivide_left_tri_vv_vari;

template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
class mdivide_left_tri_dv_vari;

template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
class mdivide_left_tri_vd_vari;

class modified_bessel_first_kind_dv_vari;

class modified_bessel_second_kind_dv_vari;

class multiply_log_vv_vari;

class multiply_log_vd_vari;

class multiply_log_dv_vari;

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class quad_form_vari_alloc;

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class quad_form_vari;

class rising_factorial_vd_vari;

class squared_distance_vv_vari;

class squared_distance_vd_vari;

template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class trace_gen_quad_form_vari_alloc;

template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class trace_gen_quad_form_vari;

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class trace_quad_form_vari_alloc;

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class trace_quad_form_vari;

template <typename ReduceFunction, typename Enable, typename ReturnType, typename Vec, typename... Args>
struct reduce_sum_impl;

template <typename ReduceFunction, typename Enable, typename ReturnType, typename Vec, typename... Args>
struct reduce_sum_impl;

 inline  char * eight_byte_aligned_malloc(size_t size);

template <typename VecVar, require_std_vector_vt<is_var, VecVar> * , typename... Pargs>
 inline  size_t count_vars_impl(size_t count, VecVar && x, Pargs &&... args);

template <typename VecContainer, require_std_vector_st<is_var, VecContainer> * , require_std_vector_vt<is_container, VecContainer> * , typename... Pargs>
 inline  size_t count_vars_impl(size_t count, VecContainer && x, Pargs &&... args);

template <typename EigT, require_eigen_vt<is_var, EigT> * , typename... Pargs>
 inline  size_t count_vars_impl(size_t count, EigT && x, Pargs &&... args);

template <typename... Pargs>
 inline  size_t count_vars_impl(size_t count, const var & x, Pargs &&... args);

template <typename Arith, require_arithmetic_t<scalar_type_t<Arith>> * , typename... Pargs>
 inline  size_t count_vars_impl(size_t count, Arith & x, Pargs &&... args);

 inline  size_t count_vars_impl(size_t count);

template <typename VecVar, require_std_vector_vt<is_var, VecVar> * , typename... Pargs>
 inline  size_t count_vars_impl(size_t count, VecVar && x, Pargs &&... args);

template <typename VecContainer, require_std_vector_st<is_var, VecContainer> * , require_std_vector_vt<is_container, VecContainer> * , typename... Pargs>
 inline  size_t count_vars_impl(size_t count, VecContainer && x, Pargs &&... args);

template <typename EigT, require_eigen_vt<is_var, EigT> * , typename... Pargs>
 inline  size_t count_vars_impl(size_t count, EigT && x, Pargs &&... args);

template <typename... Pargs>
 inline  size_t count_vars_impl(size_t count, const var & x, Pargs &&... args);

template <typename Arith, require_arithmetic_t<scalar_type_t<Arith>> * , typename... Pargs>
 inline  size_t count_vars_impl(size_t count, Arith & x, Pargs &&... args);

 inline  size_t count_vars_impl(size_t count);

template <typename LMat, typename LAMat>
 inline  void initialize_return(LMat & L, const LAMat & L_A, vari *& dummy);

template <typename T1, typename T2, typename T3>
 inline  auto unblocked_cholesky_lambda(T1 & L_A, T2 & L, T3 & A);

template <typename T1, typename T2, typename T3>
 inline  auto cholesky_lambda(T1 & L_A, T2 & L, T3 & A);

template <typename Result_, typename WMat_, typename B_>
 inline  void make_csr_adjoint(Result_ && res, WMat_ && w_mat, B_ && b);

template <typename T1, typename T2, typename T3, typename T4, require_all_matrix_t<T1, T2, T3> * >
 inline  auto fma_reverse_pass(T1 & arena_x, T2 & arena_y, T3 & arena_z, T4 & ret);

template <typename T1, typename T2, typename T3, typename T4, require_all_matrix_t<T2, T3> * , require_stan_scalar_t<T1> * >
 inline  auto fma_reverse_pass(T1 & arena_x, T2 & arena_y, T3 & arena_z, T4 & ret);

template <typename T1, typename T2, typename T3, typename T4, require_all_matrix_t<T1, T3> * , require_stan_scalar_t<T2> * >
 inline  auto fma_reverse_pass(T1 & arena_x, T2 & arena_y, T3 & arena_z, T4 & ret);

template <typename T1, typename T2, typename T3, typename T4, require_matrix_t<T3> * , require_all_stan_scalar_t<T1, T2> * >
 inline  auto fma_reverse_pass(T1 & arena_x, T2 & arena_y, T3 & arena_z, T4 & ret);

template <typename T1, typename T2, typename T3, typename T4, require_all_matrix_t<T1, T2> * , require_stan_scalar_t<T3> * >
 inline  auto fma_reverse_pass(T1 & arena_x, T2 & arena_y, T3 & arena_z, T4 & ret);

template <typename T1, typename T2, typename T3, typename T4, require_matrix_t<T2> * , require_all_stan_scalar_t<T1, T3> * >
 inline  auto fma_reverse_pass(T1 & arena_x, T2 & arena_y, T3 & arena_z, T4 & ret);

template <typename T1, typename T2, typename T3, typename T4, require_matrix_t<T1> * , require_all_stan_scalar_t<T2, T3> * >
 inline  auto fma_reverse_pass(T1 & arena_x, T2 & arena_y, T3 & arena_z, T4 & ret);

template <typename T1, typename T2>
 inline  auto generalized_inverse_lambda(T1 & G_arena, T2 & inv_G);

template <typename T1, typename T2, require_all_kernel_expressions_and_none_scalar_t<T1, T2> * >
 inline  void update_adjoints(var_value<T1> & x, const T2 & y, const vari & z);

template <typename T1, typename T2, require_all_kernel_expressions_and_none_scalar_t<T1, T2> * >
 inline  void update_adjoints(var_value<T1> & x, const T2 & y, const var & z);

template <typename Scalar1, typename Scalar2, require_var_t<Scalar1> * , require_not_var_matrix_t<Scalar1> * , require_arithmetic_t<Scalar2> * >
 inline  void update_adjoints(Scalar1 x, Scalar2 y, const vari & z) noexcept;

template <typename Scalar1, typename Scalar2, require_var_t<Scalar1> * , require_not_var_matrix_t<Scalar1> * , require_arithmetic_t<Scalar2> * >
 inline  void update_adjoints(Scalar1 x, Scalar2 y, const var & z) noexcept;

template <typename Matrix1, typename Matrix2, require_rev_matrix_t<Matrix1> * , require_st_arithmetic<Matrix2> * >
 inline  void update_adjoints(Matrix1 & x, const Matrix2 & y, const vari & z);

template <typename Matrix1, typename Matrix2, require_rev_matrix_t<Matrix1> * , require_st_arithmetic<Matrix2> * >
 inline  void update_adjoints(Matrix1 & x, const Matrix2 & y, const var & z);

template <typename Arith, typename Alt, require_st_arithmetic<Arith> * >
 inline constexpr  void update_adjoints(Arith && , Alt && , const vari & ) noexcept;

template <typename Arith, typename Alt, require_st_arithmetic<Arith> * >
 inline constexpr  void update_adjoints(Arith && , Alt && , const var & ) noexcept;

template <typename StdVec1, typename Vec2, require_std_vector_t<StdVec1> * , require_st_arithmetic<Vec2> * >
 inline  void update_adjoints(StdVec1 & x, const Vec2 & y, const vari & z);

template <typename StdVec1, typename Vec2, require_std_vector_t<StdVec1> * , require_st_arithmetic<Vec2> * >
 inline  void update_adjoints(StdVec1 & x, const Vec2 & y, const var & z);

template <typename Mat1, typename Mat2, require_all_matrix_t<Mat1, Mat2> * , require_any_var_matrix_t<Mat1, Mat2> * >
 inline  auto quad_form_impl(const Mat1 & A, const Mat2 & B, bool symmetric);

 inline  var calc_variance(size_t size, const var * dtrs);

template <int call_id, typename F, typename T_shared_param, typename T_job_param, require_eigen_col_vector_t<T_shared_param> * >
  Eigen::Matrix<return_type_t<T_shared_param, T_job_param>, Eigen::Dynamic, 1> map_rect_concurrent(const T_shared_param & shared_params, const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>> & job_params, const std::vector<std::vector<double>> & x_r, const std::vector<std::vector<int>> & x_i, std::ostream * msgs);

template <typename F>
  void finite_diff_hessian_auto(const F & f, const Eigen::VectorXd & x, double & fx, Eigen::VectorXd & grad_fx, Eigen::MatrixXd & hess_fx);

template <typename F>
  void finite_diff_hessian_times_vector_auto(const F & f, const Eigen::VectorXd & x, const Eigen::VectorXd & v, double & fx, Eigen::VectorXd & hvp);

template <typename FuncTangent, typename InputArg, require_not_st_fvar<InputArg> * >
 inline constexpr  double aggregate_tangent(const FuncTangent & tangent, const InputArg & arg);

template <typename FuncTangent, typename InputArg, require_st_fvar<InputArg> * >
 inline  auto aggregate_tangent(const FuncTangent & tangent, const InputArg & arg);

static  constexpr  auto sum_dx();

template <typename T1>
 inline  auto sum_dx(T1 & a);

template <typename T1, typename T2, typename... Types>
 inline  auto sum_dx(T1 & a, T2 & b, Types &... args);


} // namespace internal

template <typename T, typename >
class vari_value;

template <typename T, typename >
class var_value;

class chainable_alloc;

template <typename T>
struct arena_allocator;

template <typename MatrixType, typename >
class arena_matrix;

class stack_alloc;

template <typename ChainableT, typename ChainableAllocT>
struct AutodiffStackSingleton;

class chainable_alloc;

class vari_base;

class chainable_alloc;

template <typename T>
class chainable_object;

template <typename T>
class unsafe_chainable_object;

class vari_base;

template <typename , typename >
class var_value;

template <typename T, typename >
class vari_view;

template <typename Derived>
class vari_view_eigen;

template <typename , typename >
class var_value;

template <typename , typename >
class var_value;

template <typename T>
struct arena_allocator;

class ad_tape_observer;

class op_ddv_vari;

class op_dv_vari;

class op_dvd_vari;

class op_dvv_vari;

class gevv_vvv_vari;

template <typename EigRev, typename EigVari, typename EigDbl>
class vi_val_adj_functor;

template <typename EigRev, typename EigDbl>
class val_adj_functor;

template <typename EigVar, typename EigVari>
class vi_val_functor;

template <typename EigVar, typename EigVari>
class vi_adj_functor;

class nested_rev_autodiff;

class op_matrix_vari;

class op_vv_vari;

class op_vd_vari;

class precomp_vv_vari;

class op_vvv_vari;

class precomp_vvv_vari;

template <typename ContainerOperands, typename ContainerGradients>
class precomputed_gradients_vari_template;

class profile_info;

template <typename T>
class profile;

class ScopedChainableStack;

class stored_gradient_vari;

class op_vdd_vari;

class op_vdv_vari;

class op_vector_vari;

class op_vvd_vari;

template <typename T>
struct fvar;

template <typename EigFvar, typename EigOut>
class read_fvar_functor;

template <typename T_x, typename T_sigma, typename T_l>
class cov_exp_quad_vari;

template <typename T, int NX, int NY>
struct nlo_functor;

template <typename S>
struct hybrj_functor_solver;

template <typename F1, typename... Args>
class kinsol_system_data;

template <typename F>
struct KinsolFixedPointEnv;

struct FixedPointADJac;

template <typename fp_env_type, typename fp_jac_type>
struct FixedPointSolver;

template <int Lmm, typename F, typename T_y0, typename T_t0, typename T_ts, typename... T_Args>
class cvodes_integrator;

template <typename dae_type>
struct idas_service;

class idas_integrator;

template <typename F, typename Tyy, typename Typ, typename... T_par>
class dae_system;

template <typename F, typename T_y0, typename T_t0, typename T_ts, typename... T_Args>
class cvodes_integrator_adjoint_vari;

template <typename T, typename >
class accumulator;

template <typename Vari>
  void grad(Vari * vi);

template <typename Ta1, typename Ta2, typename Tb, typename Tz, require_all_stan_scalar_t<Ta1, Ta2, Tb, Tz> * , require_any_var_t<Ta1, Ta2, Tb, Tz> * >
 inline  return_type_t<Ta1, Ta1, Tb, Tz> hypergeometric_2F1(const Ta1 & a1, const Ta2 & a2, const Tb & b, const Tz & z);

template <typename T>
  bool is_aligned(T * ptr, unsigned int bytes_aligned);

template <typename T>
  auto make_chainable_ptr(T && obj);

template <typename T>
  auto make_unsafe_chainable_ptr(T && obj);







static  inline  bool empty_nested();

static  inline  size_t nested_size();

static  inline  void grad();

template <typename Vari>
  void grad(Vari * vi);

template <typename F>
 inline  void reverse_pass_callback(F && functor);





template <typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, const var & x, Pargs &&... args);

template <typename VarVec, require_std_vector_vt<is_var, VarVec> * , typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, VarVec && x, Pargs &&... args);

template <typename VecContainer, require_std_vector_st<is_var, VecContainer> * , require_std_vector_vt<is_container, VecContainer> * , typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, VecContainer && x, Pargs &&... args);

template <typename EigT, require_eigen_vt<is_var, EigT> * , typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, EigT && x, Pargs &&... args);

template <typename Arith, require_st_arithmetic<Arith> * , typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, Arith && x, Pargs &&... args);

 inline  double * accumulate_adjoints(double * dest);

template <typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, const var & x, Pargs &&... args);

template <typename VarVec, require_std_vector_vt<is_var, VarVec> * , typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, VarVec && x, Pargs &&... args);

template <typename VecContainer, require_std_vector_st<is_var, VecContainer> * , require_std_vector_vt<is_container, VecContainer> * , typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, VecContainer && x, Pargs &&... args);

template <typename EigT, require_eigen_vt<is_var, EigT> * , typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, EigT && x, Pargs &&... args);

template <typename Arith, require_st_arithmetic<Arith> * , typename... Pargs>
 inline  double * accumulate_adjoints(double * dest, Arith && x, Pargs &&... args);

 inline  double * accumulate_adjoints(double * dest);

template <int R, int C>
  vari ** build_vari_array(const Eigen::Matrix<var, R, C> & x);

template <typename... Pargs>
 inline  size_t count_vars(Pargs &&... args);

template <typename T, typename F>
 inline  internal::callback_vari<plain_type_t<T>, F> * make_callback_vari(T && value, F && functor);

template <typename T, typename F>
 inline  var_value<plain_type_t<T>> make_callback_var(T && value, F && functor);

template <typename Arith, typename >
 inline  Arith deep_copy_vars(Arith && arg);

 inline  auto deep_copy_vars(const var & arg);

template <typename VarVec, require_std_vector_vt<is_var, VarVec> * >
 inline  auto deep_copy_vars(VarVec && arg);

template <typename VecContainer, require_std_vector_st<is_var, VecContainer> * , require_std_vector_vt<is_container, VecContainer> * >
 inline  auto deep_copy_vars(VecContainer && arg);

template <typename EigT, require_eigen_vt<is_var, EigT> * >
 inline  auto deep_copy_vars(EigT && arg);

template <typename EigVar, typename EigVari, typename EigDbl>
 inline  void read_vi_val_adj(const EigVar & VarMat, EigVari & VariMat, EigDbl & ValMat, EigDbl & AdjMat);

template <typename EigRev, typename EigDbl>
 inline  void read_val_adj(const EigRev & VarMat, EigDbl & ValMat, EigDbl & AdjMat);

template <typename EigVar, typename EigVari, typename EigDbl>
 inline  void read_vi_val(const EigVar & VarMat, EigVari & VariMat, EigDbl & ValMat);

template <typename EigVar, typename EigVari, typename EigDbl>
 inline  void read_vi_adj(const EigVar & VarMat, EigVari & VariMat, EigDbl & AdjMat);

static  inline  void recover_memory_nested();

static  inline  void set_zero_all_adjoints_nested();

static  inline  void start_nested();

 inline  var operator+(const var & a, const var & b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  var operator+(const var & a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  var operator+(Arith a, const var & b);

template <typename VarMat1, typename VarMat2, require_all_rev_matrix_t<VarMat1, VarMat2> * >
 inline  auto add(VarMat1 && a, VarMat2 && b);

template <typename Arith, typename VarMat, require_st_arithmetic<Arith> * , require_rev_matrix_t<VarMat> * >
 inline  auto add(VarMat && a, const Arith & b);

template <typename Arith, typename VarMat, require_st_arithmetic<Arith> * , require_rev_matrix_t<VarMat> * >
 inline  auto add(const Arith & a, VarMat && b);

template <typename Var, typename EigMat, require_var_vt<std::is_arithmetic, Var> * , require_eigen_vt<std::is_arithmetic, EigMat> * >
 inline  auto add(const Var & a, const EigMat & b);

template <typename EigMat, typename Var, require_eigen_vt<std::is_arithmetic, EigMat> * , require_var_vt<std::is_arithmetic, Var> * >
 inline  auto add(const EigMat & a, const Var & b);

template <typename Var, typename VarMat, require_var_vt<std::is_arithmetic, Var> * , require_rev_matrix_t<VarMat> * >
 inline  auto add(const Var & a, VarMat && b);

template <typename Var, typename VarMat, require_var_vt<std::is_arithmetic, Var> * , require_rev_matrix_t<VarMat> * >
 inline  auto add(VarMat && a, const Var & b);

template <typename T1, typename T2, require_any_var_vt<std::is_arithmetic, T1, T2> * , require_any_arithmetic_t<T1, T2> * >
 inline  auto add(const T1 & a, const T2 & b);

template <typename T1, typename T2, require_all_var_vt<std::is_arithmetic, T1, T2> * >
 inline  auto add(const T1 & a, const T2 & b);

template <typename VarMat1, typename VarMat2, require_any_var_matrix_t<VarMat1, VarMat2> * >
 inline  auto operator+(VarMat1 && a, VarMat2 && b);

 inline  var operator-(const var & a, const var & b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  var operator-(const var & a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  var operator-(Arith a, const var & b);

template <typename VarMat1, typename VarMat2, require_all_rev_matrix_t<VarMat1, VarMat2> * >
 inline  auto subtract(const VarMat1 & a, const VarMat2 & b);

template <typename Arith, typename VarMat, require_st_arithmetic<Arith> * , require_rev_matrix_t<VarMat> * >
 inline  auto subtract(const VarMat & a, const Arith & b);

template <typename Arith, typename VarMat, require_st_arithmetic<Arith> * , require_rev_matrix_t<VarMat> * >
 inline  auto subtract(const Arith & a, const VarMat & b);

template <typename Var, typename EigMat, require_var_vt<std::is_arithmetic, Var> * , require_eigen_vt<std::is_arithmetic, EigMat> * >
 inline  auto subtract(const Var & a, const EigMat & b);

template <typename EigMat, typename Var, require_eigen_vt<std::is_arithmetic, EigMat> * , require_var_vt<std::is_arithmetic, Var> * >
 inline  auto subtract(const EigMat & a, const Var & b);

template <typename Var, typename VarMat, require_var_vt<std::is_arithmetic, Var> * , require_rev_matrix_t<VarMat> * >
 inline  auto subtract(const Var & a, const VarMat & b);

template <typename Var, typename VarMat, require_rev_matrix_t<VarMat> * , require_var_vt<std::is_arithmetic, Var> * >
 inline  auto subtract(const VarMat & a, const Var & b);

template <typename T1, typename T2, require_any_var_vt<std::is_arithmetic, T1, T2> * , require_any_arithmetic_t<T1, T2> * >
 inline  auto subtract(const T1 & a, const T2 & b);

template <typename T1, typename T2, require_all_var_vt<std::is_arithmetic, T1, T2> * >
 inline  auto subtract(const T1 & a, const T2 & b);

template <typename VarMat1, typename VarMat2, require_any_var_matrix_t<VarMat1, VarMat2> * >
 inline  auto operator-(const VarMat1 & a, const VarMat2 & b);

 inline  var operator*(const var & a, const var & b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  var operator*(const var & a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  var operator*(Arith a, const var & b);

 inline  std::complex<stan::math::var> operator*(const std::complex<stan::math::var> & x, const std::complex<stan::math::var> & y);

 inline  std::complex<stan::math::var> operator*(const std::complex<double> & x, const std::complex<stan::math::var> & y);

 inline  std::complex<stan::math::var> operator*(const std::complex<stan::math::var> & x, const std::complex<double> & y);

template <typename T, require_not_same_t<T, arena_t<T>> * , require_not_container_t<T> * , require_not_matrix_cl_t<T> * >
 inline  arena_t<T> to_arena(T && a);

template <typename T, require_same_t<T, arena_t<T>> * , require_not_matrix_cl_t<T> * , require_not_std_vector_t<T> * >
 inline  std::remove_reference_t<T> to_arena(T && a);

template <typename T, require_eigen_t<T> * , require_not_same_t<T, arena_t<T>> * >
 inline  arena_t<T> to_arena(const T & a);

template <typename T>
 inline  std::vector<T, arena_allocator<T>> to_arena(const std::vector<T, arena_allocator<T>> & a);

template <typename T, require_same_t<T, arena_t<T>> * >
 inline  arena_t<std::vector<T>> to_arena(const std::vector<T> & a);

template <typename T, require_not_same_t<T, arena_t<T>> * >
 inline  arena_t<std::vector<T>> to_arena(const std::vector<T> & a);

template <bool Condition, typename T, require_not_t<std::bool_constant<Condition>> * >
 inline  T to_arena_if(T && a);

template <bool Condition, typename T, require_t<std::bool_constant<Condition>> * >
 inline  arena_t<T> to_arena_if(const T & a);

template <typename T>
 inline  auto & value_of(const var_value<T> & v);

 inline  var operator/(const var & dividend, const var & divisor);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  var operator/(const var & dividend, Arith divisor);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  var operator/(Arith dividend, const var & divisor);

template <typename Scalar, typename Mat, require_matrix_t<Mat> * , require_stan_scalar_t<Scalar> * , require_all_st_var_or_arithmetic<Scalar, Mat> * , require_any_st_var<Scalar, Mat> * >
 inline  auto divide(const Mat & m, Scalar c);

template <typename Scalar, typename Mat, require_matrix_t<Mat> * , require_stan_scalar_t<Scalar> * , require_all_st_var_or_arithmetic<Scalar, Mat> * , require_any_st_var<Scalar, Mat> * >
 inline  auto divide(Scalar c, const Mat & m);

template <typename Mat1, typename Mat2, require_all_matrix_st<is_var_or_arithmetic, Mat1, Mat2> * , require_any_matrix_st<is_var, Mat1, Mat2> * >
 inline  auto divide(const Mat1 & m1, const Mat2 & m2);

template <typename T1, typename T2, require_any_var_matrix_t<T1, T2> * >
 inline  auto operator/(const T1 & dividend, const T2 & divisor);

 inline  std::complex<var> operator/(const std::complex<var> & x1, const std::complex<var> & x2);

template <typename T, require_var_t<T> * >
 inline  bool operator==(const T & a, const T & b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator==(const var & a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator==(Arith a, const var & b);

 inline  bool operator==(const var & x, const std::complex<var> & z);

 inline  bool operator==(const std::complex<var> & z, const var & y);

 inline  bool operator>(const var & a, const var & b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator>(const var & a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator>(Arith a, const var & b);

 inline  bool operator>=(const var & a, const var & b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator>=(const var & a, Arith b);

template <typename Arith, typename Var, require_arithmetic_t<Arith> * >
 inline  bool operator>=(Arith a, const var & b);

 inline  bool operator<(const var & a, const var & b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator<(const var & a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator<(Arith a, const var & b);

 inline  bool operator<=(const var & a, const var & b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator<=(const var & a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator<=(Arith a, const var & b);

 inline  bool operator&&(const var & x, const var & y);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator&&(const var & x, Arith y);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator&&(Arith x, const var & y);

 inline  bool operator||(const var & x, const var & y);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator||(const var & x, Arith y);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator||(Arith x, const var & y);

 inline  bool operator!=(const var & a, const var & b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator!=(const var & a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * >
 inline  bool operator!=(Arith a, const var & b);

 inline  bool operator!=(const var & x, const std::complex<var> & z);

 inline  bool operator!=(const std::complex<var> & z, const var & y);

 inline  var & operator--(var & a);

 inline  var operator--(var & a, int );

 inline  var & operator++(var & a);

 inline  var operator++(var & a, int );

 inline  var operator-(const var & a);

template <typename T, require_var_matrix_t<T> * >
 inline  auto operator-(const T & a);

 inline  bool operator!(const var & x);

 inline  var operator+(const var & a);

template <typename T>
 inline  void dims(const var_value<T> & x, std::vector<int> & result);

template <typename T>
 inline  void dims(const vari_value<T> & x, std::vector<int> & result);



 inline  void print_stack(std::ostream & o);

static  inline  void recover_memory();

static  inline  void set_zero_all_adjoints();

template <typename... Pargs>
 inline  vari ** save_varis(vari ** dest, const var & x, Pargs &&... args);

template <typename VarVec, require_std_vector_vt<is_var, VarVec> * , typename... Pargs>
 inline  vari ** save_varis(vari ** dest, VarVec && x, Pargs &&... args);

template <typename VecContainer, require_std_vector_st<is_var, VecContainer> * , require_std_vector_vt<is_container, VecContainer> * , typename... Pargs>
 inline  vari ** save_varis(vari ** dest, VecContainer && x, Pargs &&... args);

template <typename EigT, require_eigen_vt<is_var, EigT> * , typename... Pargs>
 inline  vari ** save_varis(vari ** dest, EigT && x, Pargs &&... args);

template <typename Arith, require_st_arithmetic<Arith> * , typename... Pargs>
 inline  vari ** save_varis(vari ** dest, Arith && x, Pargs &&... args);

 inline  vari ** save_varis(vari ** dest);

template <typename... Pargs>
 inline  vari ** save_varis(vari ** dest, const var & x, Pargs &&... args);

template <typename VarVec, require_std_vector_vt<is_var, VarVec> * , typename... Pargs>
 inline  vari ** save_varis(vari ** dest, VarVec && x, Pargs &&... args);

template <typename VecContainer, require_std_vector_st<is_var, VecContainer> * , require_std_vector_vt<is_container, VecContainer> * , typename... Pargs>
 inline  vari ** save_varis(vari ** dest, VecContainer && x, Pargs &&... args);

template <typename EigT, require_eigen_vt<is_var, EigT> * , typename... Pargs>
 inline  vari ** save_varis(vari ** dest, EigT && x, Pargs &&... args);

template <typename Arith, require_st_arithmetic<Arith> * , typename... Pargs>
 inline  vari ** save_varis(vari ** dest, Arith && x, Pargs &&... args);

 inline  vari ** save_varis(vari ** dest);

 inline  void zero_adjoints() noexcept;

template <typename T, require_st_arithmetic<T> * >
 inline  void zero_adjoints(T & x) noexcept;

 inline  void zero_adjoints(var & x);

template <typename EigMat, require_eigen_vt<is_autodiff, EigMat> * >
 inline  void zero_adjoints(EigMat & x);

template <typename StdVec, require_std_vector_st<is_autodiff, StdVec> * >
 inline  void zero_adjoints(StdVec & x);

template <typename T>
 inline  fvar<T> operator+(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> operator+(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> operator+(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> operator/(const fvar<T> & x1, const fvar<T> & x2);

template <typename T, typename U, require_arithmetic_t<U> * >
 inline  fvar<T> operator/(const fvar<T> & x1, U x2);

template <typename T, typename U, require_arithmetic_t<U> * >
 inline  fvar<T> operator/(U x1, const fvar<T> & x2);

template <typename T>
 inline  std::complex<fvar<T>> operator/(const std::complex<fvar<T>> & x1, const std::complex<fvar<T>> & x2);

template <typename T, typename U, require_arithmetic_t<U> * >
 inline  std::complex<fvar<T>> operator/(const std::complex<fvar<T>> & x1, const std::complex<U> & x2);

template <typename T>
 inline  std::complex<fvar<T>> operator/(const std::complex<fvar<T>> & x1, const fvar<T> & x2);

template <typename T, typename U, require_arithmetic_t<U> * >
 inline  std::complex<fvar<T>> operator/(const std::complex<fvar<T>> & x1, U x2);

template <typename T, typename U, require_arithmetic_t<U> * >
 inline  std::complex<fvar<T>> operator/(const std::complex<U> & x1, const std::complex<fvar<T>> & x2);

template <typename T, typename U, require_arithmetic_t<U> * >
 inline  std::complex<fvar<T>> operator/(const std::complex<U> & x1, const fvar<T> & x2);

template <typename T>
 inline  std::complex<fvar<T>> operator/(const fvar<T> & x1, const std::complex<fvar<T>> & x2);

template <typename T, typename U, typename >
 inline  std::complex<fvar<T>> operator/(const fvar<T> & x1, const std::complex<U> & x2);

template <typename T, typename U, require_arithmetic_t<U> * >
 inline  std::complex<fvar<T>> operator/(U x1, const std::complex<fvar<T>> & x2);

template <typename T>
 inline  bool operator==(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  bool operator==(const fvar<T> & x, double y);

template <typename T>
 inline  bool operator==(double x, const fvar<T> & y);

template <typename T>
 inline  bool operator>(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  bool operator>(const fvar<T> & x, double y);

template <typename T>
 inline  bool operator>(double x, const fvar<T> & y);

template <typename T>
 inline  bool operator>=(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  bool operator>=(const fvar<T> & x, double y);

template <typename T>
 inline  bool operator>=(double x, const fvar<T> & y);

template <typename T>
 inline  bool operator<(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  bool operator<(double x, const fvar<T> & y);

template <typename T>
 inline  bool operator<(const fvar<T> & x, double y);

template <typename T>
 inline  bool operator<=(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  bool operator<=(const fvar<T> & x, double y);

template <typename T>
 inline  bool operator<=(double x, const fvar<T> & y);

template <typename T>
 inline  bool operator&&(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  bool operator&&(const fvar<T> & x, double y);

template <typename T>
 inline  bool operator&&(double x, const fvar<T> & y);

template <typename T>
 inline  bool operator||(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  bool operator||(const fvar<T> & x, double y);

template <typename T>
 inline  bool operator||(double x, const fvar<T> & y);

template <typename T>
 inline  fvar<T> operator*(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  fvar<T> operator*(double x, const fvar<T> & y);

template <typename T>
 inline  fvar<T> operator*(const fvar<T> & x, double y);

template <typename T>
 inline  std::complex<stan::math::fvar<T>> operator*(const std::complex<stan::math::fvar<T>> & x, const std::complex<stan::math::fvar<T>> & y);

template <typename T>
 inline  std::complex<stan::math::fvar<T>> operator*(const std::complex<double> & x, const std::complex<stan::math::fvar<T>> & y);

template <typename T>
 inline  std::complex<stan::math::fvar<T>> operator*(const std::complex<stan::math::fvar<T>> & x, const std::complex<double> & y);

template <typename T>
 inline  bool operator!=(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  bool operator!=(const fvar<T> & x, double y);

template <typename T>
 inline  bool operator!=(double x, const fvar<T> & y);

template <typename T>
 inline  fvar<T> operator-(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> operator-(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> operator-(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> operator-(const fvar<T> & x);

template <typename T>
 inline  bool operator!(const fvar<T> & x);

template <typename T>
 inline  fvar<T> operator+(const fvar<T> & x);

template <typename EigFvar, typename EigOut>
 inline  void read_fvar(const EigFvar & FvarMat, EigOut & ValMat, EigOut & DMat);

template <typename T, typename F>
  void derivative(const F & f, const T & x, T & fx, T & dfx_dx);

template <typename F>
  void hessian(const F & f, const Eigen::Matrix<double, Eigen::Dynamic, 1> & x, double & fx, Eigen::Matrix<double, Eigen::Dynamic, 1> & grad, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & H);



template <typename F>
  void finite_diff_grad_hessian_auto(const F & f, const Eigen::VectorXd & x, double & fx, Eigen::MatrixXd & hess, std::vector<Eigen::MatrixXd> & grad_hess_fx);

template <typename F>
  void grad_hessian(const F & f, const Eigen::Matrix<double, Eigen::Dynamic, 1> & x, double & fx, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & H, std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> & grad_H);

template <typename T1, typename T2, typename F>
  void gradient_dot_vector(const F & f, const Eigen::Matrix<T1, Eigen::Dynamic, 1> & x, const Eigen::Matrix<T2, Eigen::Dynamic, 1> & v, T1 & fx, T1 & grad_fx_dot_v);

template <typename F>
  void grad_tr_mat_times_hessian(const F & f, const Eigen::Matrix<double, Eigen::Dynamic, 1> & x, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & M, Eigen::Matrix<double, Eigen::Dynamic, 1> & grad_tr_MH);

template <typename F>
  void hessian_times_vector(const F & f, const Eigen::Matrix<double, Eigen::Dynamic, 1> & x, const Eigen::Matrix<double, Eigen::Dynamic, 1> & v, double & fx, Eigen::Matrix<double, Eigen::Dynamic, 1> & Hv);

template <typename T, typename F>
  void hessian_times_vector(const F & f, const Eigen::Matrix<T, Eigen::Dynamic, 1> & x, const Eigen::Matrix<T, Eigen::Dynamic, 1> & v, T & fx, Eigen::Matrix<T, Eigen::Dynamic, 1> & Hv);

template <typename T, typename F>
  void partial_derivative(const F & f, const Eigen::Matrix<T, Eigen::Dynamic, 1> & x, int n, T & fx, T & dfx_dxn);

 inline  var Phi(const var & a);

template <typename T, require_var_matrix_t<T> * >
 inline  auto Phi(const T & a);

 inline  var Phi_approx(const var & a);

template <typename T, require_var_matrix_t<T> * >
 inline  auto Phi_approx(const T & a);

 inline  var fabs(const var & a);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto fabs(const VarMat & a);

 inline  var hypot(const var & a, const var & b);

 inline  var hypot(const var & a, double b);

 inline  var hypot(double a, const var & b);

template <typename T>
 inline  auto abs(const var_value<T> & a);

 inline  var abs(const std::complex<var> & z);

template <typename Alloc>
 inline  var sum(const std::vector<var, Alloc> & m);

template <typename T, require_rev_matrix_t<T> * >
 inline  var sum(T && x);

 inline  var atan2(const var & a, const var & b);

 inline  var atan2(const var & a, double b);

 inline  var atan2(double a, const var & b);

template <typename Mat1, typename Mat2, require_any_var_matrix_t<Mat1, Mat2> * , require_all_matrix_t<Mat1, Mat2> * >
 inline  auto atan2(const Mat1 & a, const Mat2 & b);

template <typename Scalar, typename VarMat, require_var_matrix_t<VarMat> * , require_stan_scalar_t<Scalar> * >
 inline  auto atan2(const Scalar & a, const VarMat & b);

template <typename VarMat, typename Scalar, require_var_matrix_t<VarMat> * , require_stan_scalar_t<Scalar> * >
 inline  auto atan2(const VarMat & a, const Scalar & b);

 inline  var arg(const std::complex<var> & z);

template <typename T>
 inline  auto & value_of_rec(const var_value<T> & v);

 inline  int is_inf(const var & v);

 inline  bool is_nan(const var & v);

 inline  var sinh(const var & a);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto sinh(const VarMat & a);

 inline  std::complex<var> sinh(const std::complex<var> & z);

 inline  var cos(var a);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto cos(const VarMat & a);

 inline  std::complex<var> cos(const std::complex<var> & z);

 inline  var sin(const var & a);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto sin(const VarMat & a);

 inline  std::complex<var> sin(const std::complex<var> & z);

 inline  var exp(const var & a);

 inline  std::complex<var> exp(const std::complex<var> & z);

template <typename T, require_var_matrix_t<T> * >
 inline  auto exp(const T & x);

 inline  var cosh(const var & a);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto cosh(const VarMat & a);

 inline  std::complex<var> cosh(const std::complex<var> & z);

 inline  var square(const var & x);

template <typename T, require_var_matrix_t<T> * >
 inline  auto square(const T & x);

 inline  var norm(const std::complex<var> & z);

 inline  var sqrt(const var & a);

template <typename T, require_var_matrix_t<T> * >
 inline  auto sqrt(const T & a);

 inline  std::complex<var> sqrt(const std::complex<var> & z);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto log(const var_value<T> & a);

 inline  std::complex<var> log(const std::complex<var> & z);

 inline  std::complex<var> polar(const var & r, const var & theta);

template <typename T>
 inline  std::complex<var> polar(T r, const var & theta);

template <typename T>
 inline  std::complex<var> polar(const var & r, T theta);

 inline  var asinh(const var & x);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto asinh(const VarMat & x);

 inline  std::complex<var> asinh(const std::complex<var> & z);

 inline  var asin(const var & x);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto asin(const VarMat & x);

 inline  std::complex<var> asin(const std::complex<var> & z);

 inline  var acos(const var & x);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto acos(const VarMat & x);

 inline  std::complex<var> acos(const std::complex<var> & x);

 inline  var acosh(const var & x);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto acosh(const VarMat & x);

 inline  std::complex<var> acosh(const std::complex<var> & z);

template <typename T1, typename T2, require_any_var_matrix_t<T1, T2> * >
 inline  auto append_col(const T1 & A, const T2 & B);

template <typename Scal, typename RowVec, require_stan_scalar_t<Scal> * , require_t<is_eigen_row_vector<RowVec>> * >
 inline  auto append_col(const Scal & A, const var_value<RowVec> & B);

template <typename RowVec, typename Scal, require_t<is_eigen_row_vector<RowVec>> * , require_stan_scalar_t<Scal> * >
 inline  auto append_col(const var_value<RowVec> & A, const Scal & B);

template <typename T1, typename T2, require_any_var_matrix_t<T1, T2> * >
 inline  auto append_row(const T1 & A, const T2 & B);

template <typename Scal, typename ColVec, require_stan_scalar_t<Scal> * , require_t<is_eigen_col_vector<ColVec>> * >
 inline  auto append_row(const Scal & A, const var_value<ColVec> & B);

template <typename ColVec, typename Scal, require_t<is_eigen_col_vector<ColVec>> * , require_stan_scalar_t<Scal> * >
 inline  auto append_row(const var_value<ColVec> & A, const Scal & B);

 inline  int as_bool(const var & v);

template <typename T, require_var_matrix_t<T> * >
 inline  auto as_array_or_scalar(T && v);

template <typename T, require_var_col_vector_t<T> * >
 inline  auto as_column_vector_or_scalar(T && a);

template <typename T, require_var_row_vector_t<T> * , require_not_var_col_vector_t<T> * >
 inline  auto as_column_vector_or_scalar(T && a);

 inline  var atanh(const var & x);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto atanh(const VarMat & x);

 inline  std::complex<var> atanh(const std::complex<var> & z);

 inline  var atan(const var & x);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto atan(const VarMat & x);

 inline  std::complex<var> atan(const std::complex<var> & z);

 inline  var bessel_first_kind(int v, const var & a);

template <typename T1, typename T2, require_st_integral<T1> * , require_eigen_t<T2> * >
 inline  auto bessel_first_kind(const T1 & v, const var_value<T2> & a);

 inline  var bessel_second_kind(int v, const var & a);

template <typename T1, typename T2, require_st_integral<T1> * , require_eigen_t<T2> * >
 inline  auto bessel_second_kind(const T1 & v, const var_value<T2> & a);

 inline  var beta(const var & a, const var & b);

 inline  var beta(const var & a, double b);

 inline  var beta(double a, const var & b);

template <typename Mat1, typename Mat2, require_any_var_matrix_t<Mat1, Mat2> * , require_all_matrix_t<Mat1, Mat2> * >
 inline  auto beta(const Mat1 & a, const Mat2 & b);

template <typename Scalar, typename VarMat, require_var_matrix_t<VarMat> * , require_stan_scalar_t<Scalar> * >
 inline  auto beta(const Scalar & a, const VarMat & b);

template <typename VarMat, typename Scalar, require_var_matrix_t<VarMat> * , require_stan_scalar_t<Scalar> * >
 inline  auto beta(const VarMat & a, const Scalar & b);

 inline  var binary_log_loss(int y, const var & y_hat);

template <typename Mat, require_eigen_t<Mat> * >
 inline  auto binary_log_loss(int y, const var_value<Mat> & y_hat);

template <typename StdVec, typename Mat, require_eigen_t<Mat> * , require_st_integral<StdVec> * >
 inline  auto binary_log_loss(const StdVec & y, const var_value<Mat> & y_hat);

 inline  var cbrt(const var & a);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto cbrt(const VarMat & a);

 inline  var ceil(const var & a);

template <typename T, require_matrix_t<T> * >
 inline  auto ceil(const var_value<T> & a);

template <typename EigMat, require_eigen_vt<is_var, EigMat> * >
 inline  auto cholesky_decompose(const EigMat & A);

template <typename T, require_var_matrix_t<T> * >
 inline  auto cholesky_decompose(const T & A);

template <typename EigVec, require_rev_vector_t<EigVec> * >
 inline  auto cumulative_sum(const EigVec & x);

template <typename T1, typename T2, require_all_vector_t<T1, T2> * , require_not_complex_t<return_type_t<T1, T2>> * , require_all_not_std_vector_t<T1, T2> * , require_any_st_var<T1, T2> * >
 inline  var dot_product(const T1 & v1, const T2 & v2);

template <typename Mat1, typename Mat2, require_all_eigen_t<Mat1, Mat2> * , require_any_eigen_vt<is_var, Mat1, Mat2> * >
 inline  Eigen::Matrix<return_type_t<Mat1, Mat2>, 1, Mat1::ColsAtCompileTime> columns_dot_product(const Mat1 & v1, const Mat2 & v2);

template <typename Mat1, typename Mat2, require_all_matrix_t<Mat1, Mat2> * , require_any_var_matrix_t<Mat1, Mat2> * >
 inline  auto columns_dot_product(const Mat1 & v1, const Mat2 & v2);

template <typename T, require_eigen_vector_vt<is_var, T> * >
 inline  var dot_self(const T & v);

template <typename T, require_var_matrix_t<T> * >
 inline  var dot_self(const T & v);

template <typename Mat, require_eigen_vt<is_var, Mat> * >
 inline  Eigen::Matrix<var, 1, Mat::ColsAtCompileTime> columns_dot_self(const Mat & x);

template <typename Mat, require_var_matrix_t<Mat> * >
 inline  auto columns_dot_self(const Mat & x);

 inline  std::complex<var> conj(const std::complex<var> & z);

template <typename T, require_var_t<T> * >
  auto & adjoint_of(const T & x);

template <typename T, require_not_var_t<T> * >
  internal::nonexisting_adjoint adjoint_of(const T & x);

template <typename T_x, typename T_sigma, require_st_arithmetic<T_x> * , require_stan_scalar_t<T_sigma> * >
 inline  Eigen::Matrix<var, -1, -1> gp_exp_quad_cov(const std::vector<T_x> & x, const T_sigma sigma, const var length_scale);

template <typename T_x, typename >
 inline  Eigen::Matrix<var, -1, -1> cov_exp_quad(const std::vector<T_x> & x, const var & sigma, const var & l);

template <typename T_x, typename >
 inline  Eigen::Matrix<var, -1, -1> cov_exp_quad(const std::vector<T_x> & x, double sigma, const var & l);

template <int Options, typename VarMatrix, typename Vec1, typename Vec2, require_var_t<VarMatrix> * , require_eigen_dense_base_t<value_type_t<VarMatrix>> * , require_all_std_vector_vt<std::is_integral, Vec1, Vec2> * >
 inline  auto to_soa_sparse_matrix(int m, int n, VarMatrix && w, Vec1 && u, Vec2 && v);

template <int Options, typename MatrixVar, typename Vec1, typename Vec2, require_eigen_dense_base_vt<is_var, MatrixVar> * , require_all_std_vector_vt<std::is_integral, Vec1, Vec2> * >
 inline  auto to_soa_sparse_matrix(int m, int n, MatrixVar && w, Vec1 && u, Vec2 && v);

template <int Options, typename Mat, typename Vec1, typename Vec2, require_eigen_dense_base_vt<std::is_arithmetic, Mat> * , require_all_std_vector_vt<std::is_integral, Vec1, Vec2> * >
 inline  auto to_soa_sparse_matrix(int m, int n, Mat && w, Vec1 && u, Vec2 && v);

template <typename T1, typename T2, require_any_rev_matrix_t<T1, T2> * >
 inline  auto csr_matrix_times_vector(int m, int n, const T1 & w, const std::vector<int> & v, const std::vector<int> & u, const T2 & b);

template <typename T, require_rev_matrix_t<T> * >
 inline  var determinant(const T & m);

template <typename T1, typename T2, require_vector_t<T1> * , require_matrix_t<T2> * , require_any_st_var<T1, T2> * >
  auto diag_pre_multiply(const T1 & m1, const T2 & m2);

template <typename T1, typename T2, require_matrix_t<T1> * , require_vector_t<T2> * , require_any_st_var<T1, T2> * >
  auto diag_post_multiply(const T1 & m1, const T2 & m2);

 inline  var digamma(const var & a);

template <typename T, require_var_matrix_t<T> * >
 inline  auto digamma(const T & a);

template <typename T, require_rev_matrix_t<T> * >
 inline  auto eigendecompose_sym(const T & m);

template <typename T, require_rev_matrix_t<T> * >
 inline  auto eigenvalues_sym(const T & m);

template <typename T, require_rev_matrix_t<T> * >
 inline  auto eigenvectors_sym(const T & m);

template <typename Mat1, typename Mat2, require_all_matrix_t<Mat1, Mat2> * , require_any_rev_matrix_t<Mat1, Mat2> * >
  auto elt_divide(const Mat1 & m1, const Mat2 & m2);

template <typename Scal, typename Mat, require_stan_scalar_t<Scal> * , require_var_matrix_t<Mat> * >
  auto elt_divide(Scal s, const Mat & m);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * , require_return_type_t<is_var, T1, T2> * , require_not_row_and_col_vector_t<T1, T2> * >
 inline  auto multiply(T1 && A, T2 && B);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * , require_return_type_t<is_var, T1, T2> * , require_row_and_col_vector_t<T1, T2> * >
 inline  var multiply(const T1 & A, const T2 & B);

template <typename T1, typename T2, require_not_matrix_t<T1> * , require_matrix_t<T2> * , require_return_type_t<is_var, T1, T2> * , require_not_row_and_col_vector_t<T1, T2> * >
 inline  auto multiply(const T1 & a, T2 && B);

template <typename T1, typename T2, require_matrix_t<T1> * , require_not_matrix_t<T2> * , require_any_st_var<T1, T2> * , require_not_complex_t<value_type_t<T1>> * , require_not_complex_t<value_type_t<T2>> * , require_not_row_and_col_vector_t<T1, T2> * >
 inline  auto multiply(T1 && A, T2 && B);

template <typename T1, typename T2, require_any_var_matrix_t<T1, T2> * >
 inline  auto operator*(T1 && a, T2 && b);

template <typename Mat1, typename Mat2, require_all_matrix_t<Mat1, Mat2> * , require_any_rev_matrix_t<Mat1, Mat2> * >
  auto elt_multiply(const Mat1 & m1, const Mat2 & m2);

 inline  var erf(const var & a);

template <typename T, require_matrix_t<T> * >
 inline  auto erf(const var_value<T> & a);

 inline  var erfc(const var & a);

template <typename T, require_matrix_t<T> * >
 inline  auto erfc(const var_value<T> & a);

 inline  var exp2(const var & a);

template <typename T, require_eigen_t<T> * >
 inline  auto exp2(const var_value<T> & a);

 inline  var expm1(const var & a);

template <typename T, require_eigen_t<T> * >
 inline  auto expm1(const var_value<T> & a);

 inline  var falling_factorial(const var & a, int b);

template <typename T1, typename T2, require_eigen_t<T1> * , require_st_integral<T2> * >
 inline  auto falling_factorial(const var_value<T1> & a, const T2 & b);

 inline  var fdim(const var & a, const var & b);

 inline  var fdim(double a, const var & b);

 inline  var fdim(const var & a, double b);

template <typename V, require_eigen_vector_vt<is_complex, V> * , require_var_t<base_type_t<value_type_t<V>>> * >
 inline  plain_type_t<V> fft(const V & x);

template <typename V, require_eigen_vector_vt<is_complex, V> * , require_var_t<base_type_t<value_type_t<V>>> * >
 inline  plain_type_t<V> inv_fft(const V & y);

template <typename M, require_eigen_dense_dynamic_vt<is_complex, M> * , require_var_t<base_type_t<value_type_t<M>>> * >
 inline  plain_type_t<M> fft2(const M & x);

template <typename M, require_eigen_dense_dynamic_vt<is_complex, M> * , require_var_t<base_type_t<value_type_t<M>>> * >
 inline  plain_type_t<M> inv_fft2(const M & y);

template <typename VarMat, typename S, require_var_matrix_t<VarMat> * , require_var_t<S> * >
 inline  void fill(VarMat & x, const S & y);

template <typename VarMat, typename S, require_var_matrix_t<VarMat> * , require_arithmetic_t<S> * >
 inline  void fill(VarMat & x, const S & y);

 inline  var floor(const var & a);

template <typename T, require_eigen_t<T> * >
 inline  auto floor(const var_value<T> & a);

 inline  var fma(const var & x, const var & y, const var & z);

template <typename Tc, require_arithmetic_t<Tc> * >
 inline  var fma(const var & x, const var & y, Tc && z);

template <typename Tb, require_arithmetic_t<Tb> * >
 inline  var fma(const var & x, Tb && y, const var & z);

template <typename Tb, typename Tc, require_all_arithmetic_t<Tb, Tc> * >
 inline  var fma(const var & x, Tb && y, Tc && z);

template <typename Ta, typename Tc, require_all_arithmetic_t<Ta, Tc> * >
 inline  var fma(Ta && x, const var & y, Tc && z);

template <typename Ta, typename Tb, require_all_arithmetic_t<Ta, Tb> * >
 inline  var fma(Ta && x, Tb && y, const var & z);

template <typename Ta, require_arithmetic_t<Ta> * >
 inline  var fma(Ta && x, const var & y, const var & z);

template <typename T1, typename T2, typename T3, require_any_matrix_t<T1, T2, T3> * , require_var_t<return_type_t<T1, T2, T3>> * >
 inline  auto fma(const T1 & x, const T2 & y, const T3 & z);

 inline  var fmax(const var & a, const var & b);

 inline  var fmax(const var & a, double b);

 inline  var fmax(double a, const var & b);

 inline  var fmin(const var & a, const var & b);

 inline  var fmin(const var & a, double b);

 inline  var fmin(double a, const var & b);

 inline  var fmod(const var & a, const var & b);

 inline  var fmod(const var & a, double b);

 inline  var fmod(double a, const var & b);

template <typename T, require_var_matrix_t<T> * >
  Eigen::Matrix<var, T::RowsAtCompileTime, T::ColsAtCompileTime> from_var_value(const T & a);

template <typename T, require_any_t<conjunction<is_eigen<T>, is_var<scalar_type_t<T>>>, std::is_same<std::decay_t<T>, var>, bool_constant<! std::is_same<scalar_type_t<T>, var>::value>> * >
  T from_var_value(T && a);

template <typename T>
  auto from_var_value(const std::vector<T> & a);

 inline  var gamma_p(const var & a, const var & b);

 inline  var gamma_p(const var & a, double b);

 inline  var gamma_p(double a, const var & b);

 inline  var gamma_q(const var & a, const var & b);

 inline  var gamma_q(const var & a, double b);

 inline  var gamma_q(double a, const var & b);

template <typename T, require_rev_matrix_t<T> * >
 inline  auto inverse(const T & m);

template <typename VarMat, require_rev_matrix_t<VarMat> * >
 inline  auto generalized_inverse(const VarMat & G);

template <typename T_x, typename T_sigma, require_st_arithmetic<T_x> * , require_stan_scalar_t<T_sigma> * >
 inline  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> gp_periodic_cov(const std::vector<T_x> & x, const T_sigma sigma, const var l, const var p);

 inline  void grad(var & v, Eigen::Matrix<var, Eigen::Dynamic, 1> & x, Eigen::VectorXd & g);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto inv(const var_value<T> & a);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto inv_sqrt(const var_value<T> & a);

 inline  var inv_square(const var & a);

template <typename Scal1, typename Scal2, require_any_var_t<base_type_t<Scal1>, base_type_t<Scal2>> * , require_all_stan_scalar_t<Scal1, Scal2> * >
 inline  var pow(const Scal1 & base, const Scal2 & exponent);

template <typename Mat1, typename Mat2, require_all_st_var_or_arithmetic<Mat1, Mat2> * , require_any_matrix_st<is_var, Mat1, Mat2> * , require_all_not_stan_scalar_t<Mat1, Mat2> * >
 inline  auto pow(const Mat1 & base, const Mat2 & exponent);

template <typename Mat1, typename Scal1, require_all_st_var_or_arithmetic<Mat1, Scal1> * , require_all_matrix_st<is_var, Mat1> * , require_stan_scalar_t<Scal1> * >
 inline  auto pow(const Mat1 & base, const Scal1 & exponent);

template <typename Scal1, typename Mat1, require_all_st_var_or_arithmetic<Scal1, Mat1> * , require_stan_scalar_t<Scal1> * , require_all_matrix_st<is_var, Mat1> * >
 inline  auto pow(Scal1 base, const Mat1 & exponent);

template <typename T1, typename T2, require_any_container_t<T1, T2> * , require_all_not_matrix_st<is_var, T1, T2> * , require_any_var_t<base_type_t<T1>, base_type_t<T2>> * >
 inline  auto pow(const T1 & a, const T2 & b);

 inline  var inc_beta(const var & a, const var & b, const var & c);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto lgamma(const var_value<T> & a);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto log1m(const var_value<T> & a);

template <typename Ta1, typename Ta2, typename Tb, typename Tz, require_all_stan_scalar_t<Ta1, Ta2, Tb, Tz> * , require_any_var_t<Ta1, Ta2, Tb, Tz> * >
 inline  return_type_t<Ta1, Ta1, Tb, Tz> hypergeometric_2F1(const Ta1 & a1, const Ta2 & a2, const Tb & b, const Tz & z);

 inline  void grad_inc_beta(var & g1, var & g2, const var & a, const var & b, const var & z);

template <typename Ta, typename Tz, require_all_stan_scalar_t<Ta, Tz> * , require_any_var_t<Ta, Tz> * >
  var hypergeometric_1f0(const Ta & a, const Tz & z);

template <typename Ta, typename Tb, typename Tz, bool grad_a, bool grad_b, bool grad_z, require_all_matrix_t<Ta, Tb> * , require_return_type_t<is_var, Ta, Tb, Tz> * >
 inline  var hypergeometric_pFq(const Ta & a, const Tb & b, const Tz & z);

 inline  var if_else(bool c, const var & y_true, const var & y_false);

 inline  var if_else(bool c, double y_true, const var & y_false);

 inline  var if_else(bool c, const var & y_true, double y_false);

template <typename T1, typename T2, typename T3, require_all_stan_scalar_t<T1, T2, T3> * , require_any_var_t<T1, T2, T3> * >
 inline  var inv_inc_beta(const T1 & a, const T2 & b, const T3 & p);

template <typename VarMat, typename S, require_var_matrix_t<VarMat> * , require_stan_scalar_t<S> * >
 inline  void initialize_fill(VarMat & x, const S & y);

 inline  void initialize_variable(var & variable, const var & value);

template <int R, int C>
 inline  void initialize_variable(Eigen::Matrix<var, R, C> & matrix, const var & value);

template <typename T>
 inline  void initialize_variable(std::vector<T> & variables, const var & value);

 inline  var inv_Phi(const var & p);

template <typename T, require_var_matrix_t<T> * >
 inline  auto inv_Phi(const T & p);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto inv_cloglog(const var_value<T> & a);

 inline  var inv_erfc(const var & a);

template <typename T, require_matrix_t<T> * >
 inline  auto inv_erfc(const var_value<T> & a);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto inv_logit(const var_value<T> & a);

 inline  bool is_uninitialized(var x);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto lambert_w0(const var_value<T> & a);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto lambert_wm1(const var_value<T> & a);

 inline  var lbeta(const var & a, const var & b);

 inline  var lbeta(const var & a, double b);

 inline  var lbeta(double a, const var & b);

 inline  var ldexp(const var & a, int b);

 inline  var lmgamma(int a, const var & b);

 inline  var lmultiply(const var & a, const var & b);

 inline  var lmultiply(const var & a, double b);

 inline  var lmultiply(double a, const var & b);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * , require_any_var_matrix_t<T1, T2> * >
 inline  auto lmultiply(const T1 & a, const T2 & b);

template <typename T1, typename T2, require_var_matrix_t<T1> * , require_stan_scalar_t<T2> * >
 inline  auto lmultiply(const T1 & a, const T2 & b);

template <typename T1, typename T2, require_stan_scalar_t<T1> * , require_var_matrix_t<T2> * >
 inline  auto lmultiply(const T1 & a, const T2 & b);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto log10(const var_value<T> & a);

 inline  std::complex<var> log10(const std::complex<var> & z);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto log1m_exp(const var_value<T> & x);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto log1m_inv_logit(const var_value<T> & u);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto log1p(const var_value<T> & a);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto log1p_exp(const var_value<T> & a);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto log2(const var_value<T> & a);

template <typename T, require_rev_matrix_t<T> * >
 inline  var log_determinant(const T & m);

template <typename T, require_rev_matrix_t<T> * >
  var log_determinant_ldlt(LDLT_factor<T> & A);

template <typename T, require_rev_matrix_t<T> * >
 inline  var log_determinant_spd(const T & M);

 inline  var log_diff_exp(const var & a, const var & b);

 inline  var log_diff_exp(const var & a, double b);

 inline  var log_diff_exp(double a, const var & b);

 inline  var log_falling_factorial(const var & a, double b);

 inline  var log_falling_factorial(const var & a, const var & b);

 inline  var log_falling_factorial(double a, const var & b);

template <typename T, require_arithmetic_t<T> * >
 inline  auto log_inv_logit(const var_value<T> & u);

template <typename T, require_eigen_t<T> * >
 inline  auto log_inv_logit(const var_value<T> & u);

 inline  var log_inv_logit_diff(const var & a, double b);

 inline  var log_inv_logit_diff(const var & a, const var & b);

 inline  var log_inv_logit_diff(double a, const var & b);

 inline  void log_mix_partial_helper(double theta_val, double lambda1_val, double lambda2_val, double & one_m_exp_lam2_m_lam1, double & one_m_t_prod_exp_lam2_m_lam1, double & one_d_t_plus_one_m_t_prod_exp_lam2_m_lam1);

template <typename T_theta, typename T_lambda1, typename T_lambda2, require_any_var_t<T_theta, T_lambda1, T_lambda2> * >
 inline  return_type_t<T_theta, T_lambda1, T_lambda2> log_mix(const T_theta & theta, const T_lambda1 & lambda1, const T_lambda2 & lambda2);

 inline  var log_rising_factorial(const var & a, double b);

 inline  var log_rising_factorial(const var & a, const var & b);

 inline  var log_rising_factorial(double a, const var & b);

template <typename T, require_eigen_st<is_var, T> * >
  auto log_softmax(const T & x);

template <typename T, require_var_matrix_t<T> * >
 inline  auto log_softmax(const T & x);

template <typename T, require_std_vector_st<is_var, T> * >
 inline  auto log_softmax(const T & x);

 inline  var log_sum_exp(const var & a, const var & b);

 inline  var log_sum_exp(const var & a, double b);

 inline  var log_sum_exp(double a, const var & b);

template <typename T, require_eigen_st<is_var, T> * , require_not_var_matrix_t<T> * >
 inline  var log_sum_exp(const T & v);

template <typename T, require_var_matrix_t<T> * >
 inline  var log_sum_exp(const T & x);

template <typename T, require_std_vector_st<is_var, T> * >
 inline  auto log_sum_exp(const T & x);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto logit(const var_value<T> & u);

template <typename Ta, typename Tb, require_all_eigen_t<Ta, Tb> * , require_any_st_autodiff<Ta, Tb> * >
 inline  Eigen::Matrix<return_type_t<Ta, Tb>, -1, Tb::ColsAtCompileTime> matrix_exp_multiply(const Ta & A, const Tb & B);

template <typename T, require_rev_matrix_t<T> * >
 inline  plain_type_t<T> matrix_power(const T & M, const int n);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * , require_any_st_var<T1, T2> * >
 inline  auto mdivide_left(const T1 & A, const T2 & B);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * , require_any_st_var<T1, T2> * >
 inline  auto mdivide_left_ldlt(LDLT_factor<T1> & A, const T2 & B);

template <typename EigMat1, typename EigMat2, require_all_eigen_matrix_base_vt<is_var, EigMat1, EigMat2> * >
 inline  Eigen::Matrix<var, EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime> mdivide_left_spd(const EigMat1 & A, const EigMat2 & b);

template <typename EigMat1, typename EigMat2, require_eigen_matrix_base_vt<is_var, EigMat1> * , require_eigen_matrix_base_vt<std::is_arithmetic, EigMat2> * >
 inline  Eigen::Matrix<var, EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime> mdivide_left_spd(const EigMat1 & A, const EigMat2 & b);

template <typename EigMat1, typename EigMat2, require_eigen_matrix_base_vt<std::is_arithmetic, EigMat1> * , require_eigen_matrix_base_vt<is_var, EigMat2> * >
 inline  Eigen::Matrix<var, EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime> mdivide_left_spd(const EigMat1 & A, const EigMat2 & b);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * , require_any_var_matrix_t<T1, T2> * >
 inline  auto mdivide_left_spd(const T1 & A, const T2 & B);

template <Eigen::UpLoType TriView, typename T1, typename T2, require_all_eigen_vt<is_var, T1, T2> * >
 inline  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> mdivide_left_tri(const T1 & A, const T2 & b);

template <Eigen::UpLoType TriView, typename T1, typename T2, require_eigen_vt<std::is_arithmetic, T1> * , require_eigen_vt<is_var, T2> * >
 inline  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> mdivide_left_tri(const T1 & A, const T2 & b);

template <Eigen::UpLoType TriView, typename T1, typename T2, require_eigen_vt<is_var, T1> * , require_eigen_vt<std::is_arithmetic, T2> * >
 inline  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> mdivide_left_tri(const T1 & A, const T2 & b);

template <Eigen::UpLoType TriView, typename T1, typename T2, require_all_matrix_t<T1, T2> * , require_any_var_matrix_t<T1, T2> * >
 inline  auto mdivide_left_tri(const T1 & A, const T2 & B);

 inline  var modified_bessel_first_kind(int v, const var & a);

 inline  var modified_bessel_second_kind(int v, const var & a);

 inline  var multiply_log(const var & a, const var & b);

 inline  var multiply_log(const var & a, double b);

 inline  var multiply_log(double a, const var & b);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * , require_any_var_matrix_t<T1, T2> * >
 inline  auto multiply_log(const T1 & a, const T2 & b);

template <typename T1, typename T2, require_var_matrix_t<T1> * , require_stan_scalar_t<T2> * >
 inline  auto multiply_log(const T1 & a, const T2 & b);

template <typename T1, typename T2, require_stan_scalar_t<T1> * , require_var_matrix_t<T2> * >
 inline  auto multiply_log(const T1 & a, const T2 & b);

template <typename T, require_rev_matrix_t<T> * >
 inline  auto multiply_lower_tri_self_transpose(const T & L);

template <typename T, require_eigen_vector_vt<is_var, T> * >
 inline  var norm1(const T & v);

template <typename T, require_var_matrix_t<T> * >
 inline  var norm1(const T & v);

template <typename T, require_eigen_vector_vt<is_var, T> * >
 inline  var norm2(const T & v);

template <typename T, require_var_matrix_t<T> * >
 inline  var norm2(const T & v);

template <typename Var1, typename Var2, require_all_st_var<Var1, Var2> * , require_all_not_std_vector_t<Var1, Var2> * >
 inline  auto owens_t(const Var1 & h, const Var2 & a);

template <typename Var, typename Arith, require_st_arithmetic<Arith> * , require_all_not_std_vector_t<Var, Arith> * , require_st_var<Var> * >
 inline  auto owens_t(const Var & h, const Arith & a);

template <typename Arith, typename Var, require_st_arithmetic<Arith> * , require_all_not_std_vector_t<Var, Arith> * , require_st_var<Var> * >
 inline  auto owens_t(const Arith & h, const Var & a);

 inline  double primitive_value(const var & v);

 inline  std::complex<var> proj(const std::complex<var> & z);

template <typename T, require_eigen_vt<is_var, T> * >
 inline  var_value<Eigen::Matrix<double, T::RowsAtCompileTime, T::ColsAtCompileTime>> to_var_value(const T & a);

template <typename T, require_var_t<T> * >
 inline  T to_var_value(T && a);

template <typename T>
 inline  auto to_var_value(const std::vector<T> & a);









template <typename EigMat1, typename EigMat2, require_all_eigen_t<EigMat1, EigMat2> * , require_any_vt_var<EigMat1, EigMat2> * >
 inline  auto quad_form_sym(const EigMat1 & A, const EigMat2 & B);

template <typename T, require_var_vector_t<T> * >
  auto read_corr_L(const T & CPCs, size_t K);

template <typename T1, typename T2, require_var_vector_t<T1> * , require_stan_scalar_t<T2> * >
  auto read_corr_L(const T1 & CPCs, size_t K, T2 & log_prob);

template <typename T_CPCs, require_var_vector_t<T_CPCs> * >
 inline  var_value<Eigen::MatrixXd> read_corr_matrix(const T_CPCs & CPCs, size_t K);

template <typename T_CPCs, require_var_vector_t<T_CPCs> * >
 inline  var_value<Eigen::MatrixXd> read_corr_matrix(const T_CPCs & CPCs, size_t K, scalar_type_t<T_CPCs> & log_prob);

template <typename T_CPCs, typename T_sds, require_any_var_vector_t<T_CPCs, T_sds> * , require_vt_same<T_CPCs, T_sds> * >
 inline  auto read_cov_L(const T_CPCs & CPCs, const T_sds & sds, scalar_type_t<T_CPCs> & log_prob);

template <typename Mat1, typename Mat2, require_all_eigen_t<Mat1, Mat2> * , require_any_eigen_vt<is_var, Mat1, Mat2> * >
 inline  Eigen::Matrix<var, Mat1::RowsAtCompileTime, 1> rows_dot_product(const Mat1 & v1, const Mat2 & v2);

template <typename Mat1, typename Mat2, require_all_matrix_t<Mat1, Mat2> * , require_any_var_matrix_t<Mat1, Mat2> * >
 inline  auto rows_dot_product(const Mat1 & v1, const Mat2 & v2);

template <typename T_CPCs, typename T_sds, require_all_var_vector_t<T_CPCs, T_sds> * >
  var_value<Eigen::MatrixXd> read_cov_matrix(const T_CPCs & CPCs, const T_sds & sds, scalar_type_t<T_CPCs> & log_prob);

template <typename T_CPCs, typename T_sds, require_all_var_vector_t<T_CPCs, T_sds> * >
 inline  var_value<Eigen::MatrixXd> read_cov_matrix(const T_CPCs & CPCs, const T_sds & sds);

template <typename Ret, typename T, require_var_matrix_t<Ret> * , require_var_t<T> * >
 inline  auto rep_matrix(const T & x, int m, int n);

template <typename Ret, typename Vec, require_var_matrix_t<Ret> * , require_var_matrix_t<Vec> * >
 inline  auto rep_matrix(const Vec & x, int n);

template <typename T_ret, require_var_matrix_t<T_ret> * , require_eigen_row_vector_t<value_type_t<T_ret>> * >
 inline  auto rep_row_vector(var x, int n);

template <typename T_ret, require_var_matrix_t<T_ret> * , require_eigen_col_vector_t<value_type_t<T_ret>> * >
 inline  auto rep_vector(var x, int n);

 inline  var rising_factorial(const var & a, int b);

 inline  var round(const var & a);

template <typename Mat, require_eigen_vt<is_var, Mat> * >
 inline  Eigen::Matrix<var, Mat::RowsAtCompileTime, 1> rows_dot_self(const Mat & x);

template <typename Mat, require_var_matrix_t<Mat> * >
 inline  auto rows_dot_self(const Mat & x);

template <typename T, require_eigen_st<is_var, T> * >
  var sd(const T & x);

template <typename T, require_var_matrix_t<T> * >
  var sd(const T & x);

template <typename T, require_std_vector_st<is_var, T> * >
  auto sd(const T & m);

template <typename EigMat, require_rev_matrix_t<EigMat> * >
 inline  auto singular_values(const EigMat & m);

template <typename EigMat, require_rev_matrix_t<EigMat> * >
 inline  auto svd(const EigMat & m);

template <typename EigMat, require_rev_matrix_t<EigMat> * >
 inline  auto svd_U(const EigMat & m);

template <typename EigMat, require_rev_matrix_t<EigMat> * >
 inline  auto svd_V(const EigMat & m);

template <typename Mat, require_rev_matrix_t<Mat> * >
 inline  auto softmax(const Mat & alpha);

 inline  var squared_distance(const var & a, const var & b);

 inline  var squared_distance(const var & a, double b);

 inline  var squared_distance(double a, const var & b);

template <typename EigVecVar1, typename EigVecVar2, require_all_eigen_vector_vt<is_var, EigVecVar1, EigVecVar2> * >
 inline  var squared_distance(const EigVecVar1 & v1, const EigVecVar2 & v2);

template <typename EigVecVar, typename EigVecArith, require_eigen_vector_vt<is_var, EigVecVar> * , require_eigen_vector_vt<std::is_arithmetic, EigVecArith> * >
 inline  var squared_distance(const EigVecVar & v1, const EigVecArith & v2);

template <typename EigVecArith, typename EigVecVar, require_eigen_vector_vt<std::is_arithmetic, EigVecArith> * , require_eigen_vector_vt<is_var, EigVecVar> * >
 inline  var squared_distance(const EigVecArith & v1, const EigVecVar & v2);

template <typename T1, typename T2, require_all_vector_t<T1, T2> * , require_any_var_vector_t<T1, T2> * >
 inline  var squared_distance(const T1 & A, const T2 & B);

 inline  void stan_print(std::ostream * o, const var & x);

 inline  var step(const var & a);

 inline  var tanh(const var & a);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto tanh(const VarMat & a);

 inline  std::complex<var> tanh(const std::complex<var> & z);

 inline  var tan(const var & a);

template <typename VarMat, require_var_matrix_t<VarMat> * >
 inline  auto tan(const VarMat & a);

 inline  std::complex<var> tan(const std::complex<var> & z);

template <typename T, require_rev_matrix_t<T> * >
 inline  auto tcrossprod(const T & M);

 inline  var tgamma(const var & a);

template <typename T, require_var_matrix_t<T> * >
 inline  auto tgamma(const T & a);

 inline  var to_var(double x);

 inline  var & to_var(var & x);

 inline  const var & to_var(const var & x);

 inline  std::vector<var> to_var(const std::vector<double> & v);

 inline  const std::vector<var> & to_var(const std::vector<var> & v);

 inline  std::vector<var> & to_var(std::vector<var> & v);

 inline  matrix_v to_var(const matrix_d & m);

 inline  matrix_v & to_var(matrix_v & m);

 inline  const matrix_v & to_var(const matrix_v & m);

 inline  vector_v to_var(const vector_d & v);

 inline  const vector_v & to_var(const vector_v & v);

 inline  vector_v & to_var(vector_v & v);

 inline  row_vector_v to_var(const row_vector_d & rv);

 inline  const row_vector_v & to_var(const row_vector_v & rv);

 inline  row_vector_v & to_var(row_vector_v & rv);

template <typename EigMat, require_eigen_t<EigMat> * >
 inline  auto to_vector(const var_value<EigMat> & x);

template <typename T, require_rev_matrix_t<T> * >
 inline  auto trace(const T & m);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * , require_any_st_var<T1, T2> * >
 inline  var trace_inv_quad_form_ldlt(LDLT_factor<T1> & A, const T2 & B);

template <typename Td, typename Ta, typename Tb, require_not_col_vector_t<Td> * , require_all_matrix_t<Td, Ta, Tb> * , require_any_st_var<Td, Ta, Tb> * >
 inline  var trace_gen_inv_quad_form_ldlt(const Td & D, LDLT_factor<Ta> & A, const Tb & B);

template <typename Td, typename Ta, typename Tb, require_col_vector_t<Td> * , require_all_matrix_t<Ta, Tb> * , require_any_st_var<Td, Ta, Tb> * >
 inline  var trace_gen_inv_quad_form_ldlt(const Td & D, const LDLT_factor<Ta> & A, const Tb & B);

template <typename Td, typename Ta, typename Tb, typename , typename >
 inline  var trace_gen_quad_form(const Td & D, const Ta & A, const Tb & B);

template <typename Td, typename Ta, typename Tb, require_any_var_matrix_t<Td, Ta, Tb> * , require_all_matrix_t<Td, Ta, Tb> * >
 inline  var trace_gen_quad_form(const Td & D, const Ta & A, const Tb & B);

template <typename EigMat1, typename EigMat2, require_all_eigen_t<EigMat1, EigMat2> * , require_any_st_var<EigMat1, EigMat2> * >
 inline  return_type_t<EigMat1, EigMat2> trace_quad_form(const EigMat1 & A, const EigMat2 & B);

template <typename Mat1, typename Mat2, require_all_matrix_t<Mat1, Mat2> * , require_any_var_matrix_t<Mat1, Mat2> * >
 inline  var trace_quad_form(const Mat1 & A, const Mat2 & B);

 inline  var trigamma(const var & u);

 inline  var trunc(const var & a);

 inline  var variance(const std::vector<var> & v);

template <typename EigMat, require_eigen_vt<is_var, EigMat> * >
  var variance(const EigMat & m);

template <typename Mat, require_var_matrix_t<Mat> * >
 inline  var variance(const Mat & x);

template <typename F>
  void jacobian(const F & f, const Eigen::Matrix<double, Eigen::Dynamic, 1> & x, Eigen::Matrix<double, Eigen::Dynamic, 1> & fx, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & J);

template <typename T1, typename T2>
  void algebra_solver_check(const Eigen::Matrix<T1, Eigen::Dynamic, 1> & x, const Eigen::Matrix<T2, Eigen::Dynamic, 1> y, const std::vector<double> & dat, const std::vector<int> & dat_int, double function_tolerance, long max_num_steps);



template <typename F, typename T, typename... Args, require_eigen_vector_t<T> * >
  T & solve_powell_call_solver(const F & f, T & x, std::ostream *const msgs, const double relative_tolerance, const double function_tolerance, const int64_t max_num_steps, const Args &... args);

template <typename F, typename T, typename... Args, require_eigen_vector_t<T> * , require_all_st_arithmetic<Args...> * >
  Eigen::VectorXd solve_powell_tol(const F & f, const T & x, const double relative_tolerance, const double function_tolerance, const int64_t max_num_steps, std::ostream *const msgs, const Args &... args);

template <typename F, typename T, typename... T_Args, require_eigen_vector_t<T> * >
  Eigen::Matrix<stan::return_type_t<T_Args...>, Eigen::Dynamic, 1> solve_powell(const F & f, const T & x, std::ostream *const msgs, const T_Args &... args);



template <typename F, typename T, typename... T_Args, require_eigen_vector_t<T> * , require_any_st_var<T_Args...> * >
  Eigen::Matrix<var, Eigen::Dynamic, 1> solve_powell_tol(const F & f, const T & x, const double relative_tolerance, const double function_tolerance, const int64_t max_num_steps, std::ostream *const msgs, const T_Args &... args);

template <typename F1, typename... Args>
  Eigen::VectorXd kinsol_solve(const F1 & f, const Eigen::VectorXd & x, const double scaling_step_tol, const double function_tolerance, const int64_t max_num_steps, const bool custom_jacobian, const int steps_eval_jacobian, const int global_line_search, std::ostream *const msgs, const Args &... args);

template <typename F, typename T, typename... Args, require_eigen_vector_t<T> * , require_all_st_arithmetic<Args...> * >
  Eigen::VectorXd solve_newton_tol(const F & f, const T & x, const double scaling_step_size, const double function_tolerance, const int64_t max_num_steps, std::ostream *const msgs, const Args &... args);

template <typename F, typename T, typename... T_Args, require_eigen_vector_t<T> * , require_any_st_var<T_Args...> * >
  Eigen::Matrix<var, Eigen::Dynamic, 1> solve_newton_tol(const F & f, const T & x, const double scaling_step_size, const double function_tolerance, const int64_t max_num_steps, std::ostream *const msgs, const T_Args &... args);

template <typename F, typename T, typename... T_Args, require_eigen_vector_t<T> * >
  Eigen::Matrix<stan::return_type_t<T_Args...>, Eigen::Dynamic, 1> solve_newton(const F & f, const T & x, std::ostream *const msgs, const T_Args &... args);



template <typename T1, typename T2, typename F, require_any_var_matrix_t<T1, T2> * , require_all_matrix_t<T1, T2> * >
 inline  auto apply_scalar_binary(const T1 & x, const T2 & y, const F & f);

template <typename T1, typename T2, typename F, require_any_var_matrix_t<T1, T2> * , require_any_std_vector_vt<std::is_integral, T1, T2> * >
 inline  auto apply_scalar_binary(const T1 & x, const T2 & y, const F & f);

template <typename T1, typename T2, typename F, require_any_std_vector_vt<is_std_vector, T1, T2> * , require_any_std_vector_st<std::is_integral, T1, T2> * , require_any_var_matrix_t<T1, T2> * >
 inline  auto apply_scalar_binary(const T1 & x, const T2 & y, const F & f);

template <typename T1, typename T2, typename F, require_any_stan_scalar_t<T1, T2> * , require_any_var_matrix_t<T1, T2> * >
 inline  auto apply_scalar_binary(const T1 & x, const T2 & y, const F & f);

 inline  void cvodes_set_options(void * cvodes_mem, long max_num_steps);

template <typename F, typename T_y0_t0, typename T_t0, typename T_t, typename... Args, require_any_autodiff_t<T_y0_t0, T_t0, T_t, scalar_type_t<Args>...> * >
  Eigen::Matrix<var, Eigen::Dynamic, 1> ode_store_sensitivities(const F & f, const std::vector<double> & coupled_state, const Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1> & y0, const T_t0 & t0, const T_t & t, std::ostream * msgs, const Args &... args);

template <typename F>
  void gradient(const F & f, const Eigen::Matrix<double, Eigen::Dynamic, 1> & x, double & fx, Eigen::Matrix<double, Eigen::Dynamic, 1> & grad_fx);

template <typename F, typename EigVec, typename InputIt, require_eigen_vector_vt<std::is_arithmetic, EigVec> * >
  void gradient(const F & f, const EigVec & x, double & fx, InputIt first_grad_fx, InputIt last_grad_fx);

template <typename F, typename T_a, typename T_b, typename... Args, require_any_st_var<T_a, T_b, Args...> * >
 inline  return_type_t<T_a, T_b, Args...> integrate_1d_impl(const F & f, const T_a & a, const T_b & b, double relative_tolerance, std::ostream * msgs, const Args &... args);



template <typename F, typename T_yy, typename T_yp, typename... T_Args, require_all_eigen_col_vector_t<T_yy, T_yp> * >
  std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>> dae_tol_impl(const char * func, const F & f, const T_yy & yy0, const T_yp & yp0, double t0, const std::vector<double> & ts, double rtol, double atol, int64_t max_num_steps, std::ostream * msgs, const T_Args &... args);

template <typename F, typename T_yy, typename T_yp, typename... T_Args, require_all_eigen_col_vector_t<T_yy, T_yp> * >
  std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>> dae_tol(const F & f, const T_yy & yy0, const T_yp & yp0, double t0, const std::vector<double> & ts, double rtol, double atol, int64_t max_num_steps, std::ostream * msgs, const T_Args &... args);

template <typename F, typename T_yy, typename T_yp, typename... T_Args, require_all_eigen_col_vector_t<T_yy, T_yp> * >
  std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>> dae(const F & f, const T_yy & yy0, const T_yp & yp0, double t0, const std::vector<double> & ts, std::ostream * msgs, const T_Args &... args);

template <typename F, typename T_y0, typename T_t0, typename T_ts, typename... T_Args, require_eigen_col_vector_t<T_y0> * >
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>, Eigen::Dynamic, 1>> ode_adams_tol_impl(const char * function_name, const F & f, const T_y0 & y0, const T_t0 & t0, const std::vector<T_ts> & ts, double relative_tolerance, double absolute_tolerance, long max_num_steps, std::ostream * msgs, const T_Args &... args);

template <typename F, typename T_y0, typename T_t0, typename T_ts, typename... T_Args, require_eigen_col_vector_t<T_y0> * >
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>, Eigen::Dynamic, 1>> ode_adams_tol(const F & f, const T_y0 & y0, const T_t0 & t0, const std::vector<T_ts> & ts, double relative_tolerance, double absolute_tolerance, long max_num_steps, std::ostream * msgs, const T_Args &... args);

template <typename F, typename T_y0, typename T_t0, typename T_ts, typename... T_Args, require_eigen_col_vector_t<T_y0> * >
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>, Eigen::Dynamic, 1>> ode_adams(const F & f, const T_y0 & y0, const T_t0 & t0, const std::vector<T_ts> & ts, std::ostream * msgs, const T_Args &... args);



template <typename F, typename T_y0, typename T_t0, typename T_ts, typename... T_Args, require_eigen_col_vector_t<T_y0> * >
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>, Eigen::Dynamic, 1>> ode_bdf_tol_impl(const char * function_name, const F & f, const T_y0 & y0, const T_t0 & t0, const std::vector<T_ts> & ts, double relative_tolerance, double absolute_tolerance, long max_num_steps, std::ostream * msgs, const T_Args &... args);

template <typename F, typename T_y0, typename T_t0, typename T_ts, typename... T_Args, require_eigen_col_vector_t<T_y0> * >
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>, Eigen::Dynamic, 1>> ode_bdf_tol(const F & f, const T_y0 & y0, const T_t0 & t0, const std::vector<T_ts> & ts, double relative_tolerance, double absolute_tolerance, long max_num_steps, std::ostream * msgs, const T_Args &... args);

template <typename F, typename T_y0, typename T_t0, typename T_ts, typename... T_Args, require_eigen_col_vector_t<T_y0> * >
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>, Eigen::Dynamic, 1>> ode_bdf(const F & f, const T_y0 & y0, const T_t0 & t0, const std::vector<T_ts> & ts, std::ostream * msgs, const T_Args &... args);



template <typename F, typename T_y0, typename T_t0, typename T_ts, typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args, require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd, T_abs_tol_bwd> * , require_any_not_st_arithmetic<T_y0, T_t0, T_ts, T_Args...> * >
  auto ode_adjoint_impl(const char * function_name, F && f, const T_y0 & y0, const T_t0 & t0, const std::vector<T_ts> & ts, double relative_tolerance_forward, const T_abs_tol_fwd & absolute_tolerance_forward, double relative_tolerance_backward, const T_abs_tol_bwd & absolute_tolerance_backward, double relative_tolerance_quadrature, double absolute_tolerance_quadrature, long max_num_steps, long num_steps_between_checkpoints, int interpolation_polynomial, int solver_forward, int solver_backward, std::ostream * msgs, const T_Args &... args);

template <typename F, typename T_y0, typename T_t0, typename T_ts, typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args, require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd, T_abs_tol_bwd> * , require_all_st_arithmetic<T_y0, T_t0, T_ts, T_Args...> * >
  std::vector<Eigen::VectorXd> ode_adjoint_impl(const char * function_name, F && f, const T_y0 & y0, const T_t0 & t0, const std::vector<T_ts> & ts, double relative_tolerance_forward, const T_abs_tol_fwd & absolute_tolerance_forward, double relative_tolerance_backward, const T_abs_tol_bwd & absolute_tolerance_backward, double relative_tolerance_quadrature, double absolute_tolerance_quadrature, long max_num_steps, long num_steps_between_checkpoints, int interpolation_polynomial, int solver_forward, int solver_backward, std::ostream * msgs, const T_Args &... args);

template <typename F, typename T_y0, typename T_t0, typename T_ts, typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args, require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd, T_abs_tol_bwd> * >
  auto ode_adjoint_tol_ctl(F && f, const T_y0 & y0, const T_t0 & t0, const std::vector<T_ts> & ts, double relative_tolerance_forward, const T_abs_tol_fwd & absolute_tolerance_forward, double relative_tolerance_backward, const T_abs_tol_bwd & absolute_tolerance_backward, double relative_tolerance_quadrature, double absolute_tolerance_quadrature, long max_num_steps, long num_steps_between_checkpoints, int interpolation_polynomial, int solver_forward, int solver_backward, std::ostream * msgs, const T_Args &... args);

template <typename T, require_stan_scalar_or_eigen_t<T> * >
 inline  auto std_normal_log_qf(const var_value<T> & log_p);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> cholesky_corr_constrain(const T & y, int K);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> cholesky_corr_constrain(const T & y, int K, scalar_type_t<T> & lp);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> cholesky_factor_constrain(const T & x, int M, int N);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> cholesky_factor_constrain(const T & x, int M, int N, scalar_type_t<T> & lp);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> corr_matrix_constrain(const T & x, Eigen::Index k);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> corr_matrix_constrain(const T & x, Eigen::Index k, scalar_type_t<T> & lp);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> cov_matrix_constrain(const T & x, Eigen::Index K);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> cov_matrix_constrain(const T & x, Eigen::Index K, scalar_type_t<T> & lp);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> cov_matrix_constrain_lkj(const T & x, size_t k);

template <typename T, require_var_vector_t<T> * >
  var_value<Eigen::MatrixXd> cov_matrix_constrain_lkj(const T & x, size_t k, scalar_type_t<T> & lp);

template <typename T, typename... Types, require_eigen_vt<std::is_arithmetic, T> * , require_any_var_matrix_t<T, Types...> * >
 inline  auto identity_constrain(T && x, Types &&... );

template <typename T, typename... Types, require_eigen_vt<is_var, T> * , require_any_var_matrix_t<T, Types...> * >
 inline  auto identity_constrain(T && x, Types &&... );

template <typename T, typename... Types, require_var_matrix_t<T> * , require_any_var_matrix_t<T, Types...> * >
 inline  auto identity_constrain(T && x, Types &&... );

template <typename T, typename L, require_all_stan_scalar_t<T, L> * , require_any_var_t<T, L> * >
 inline  auto lb_constrain(const T & x, const L & lb);

template <typename T, typename L, require_all_stan_scalar_t<T, L> * , require_any_var_t<T, L> * >
 inline  auto lb_constrain(const T & x, const L & lb, var & lp);

template <typename T, typename L, require_matrix_t<T> * , require_stan_scalar_t<L> * , require_any_st_var<T, L> * >
 inline  auto lb_constrain(const T & x, const L & lb);

template <typename T, typename L, require_matrix_t<T> * , require_stan_scalar_t<L> * , require_any_st_var<T, L> * >
 inline  auto lb_constrain(const T & x, const L & lb, return_type_t<T, L> & lp);

template <typename T, typename L, require_all_matrix_t<T, L> * , require_any_st_var<T, L> * >
 inline  auto lb_constrain(const T & x, const L & lb);

template <typename T, typename L, require_all_matrix_t<T, L> * , require_any_st_var<T, L> * >
 inline  auto lb_constrain(const T & x, const L & lb, return_type_t<T, L> & lp);

template <typename T, typename U, require_all_stan_scalar_t<T, U> * , require_any_var_t<T, U> * >
 inline  auto ub_constrain(const T & x, const U & ub);

template <typename T, typename U, require_all_stan_scalar_t<T, U> * , require_any_var_t<T, U> * >
 inline  auto ub_constrain(const T & x, const U & ub, return_type_t<T, U> & lp);

template <typename T, typename U, require_matrix_t<T> * , require_stan_scalar_t<U> * , require_any_st_var<T, U> * >
 inline  auto ub_constrain(const T & x, const U & ub);

template <typename T, typename U, require_matrix_t<T> * , require_stan_scalar_t<U> * , require_any_st_var<T, U> * >
 inline  auto ub_constrain(const T & x, const U & ub, return_type_t<T, U> & lp);

template <typename T, typename U, require_all_matrix_t<T, U> * , require_any_st_var<T, U> * >
 inline  auto ub_constrain(const T & x, const U & ub);

template <typename T, typename U, require_all_matrix_t<T, U> * , require_any_st_var<T, U> * >
 inline  auto ub_constrain(const T & x, const U & ub, return_type_t<T, U> & lp);

template <typename T, typename L, typename U, require_all_stan_scalar_t<T, L, U> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub);

template <typename T, typename L, typename U, require_all_stan_scalar_t<T, L, U> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub, return_type_t<T, L, U> & lp);

template <typename T, typename L, typename U, require_matrix_t<T> * , require_all_stan_scalar_t<L, U> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub);

template <typename T, typename L, typename U, require_matrix_t<T> * , require_all_stan_scalar_t<L, U> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub, return_type_t<T, L, U> & lp);

template <typename T, typename L, typename U, require_all_matrix_t<T, L> * , require_stan_scalar_t<U> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub);

template <typename T, typename L, typename U, require_all_matrix_t<T, L> * , require_stan_scalar_t<U> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub, std::decay_t<return_type_t<T, L, U>> & lp);

template <typename T, typename L, typename U, require_all_matrix_t<T, U> * , require_stan_scalar_t<L> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub);

template <typename T, typename L, typename U, require_all_matrix_t<T, U> * , require_stan_scalar_t<L> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub, std::decay_t<return_type_t<T, L, U>> & lp);

template <typename T, typename L, typename U, require_all_matrix_t<T, L, U> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub);

template <typename T, typename L, typename U, require_all_matrix_t<T, L, U> * , require_var_t<return_type_t<T, L, U>> * >
 inline  auto lub_constrain(const T & x, const L & lb, const U & ub, return_type_t<T, L, U> & lp);

template <typename T, require_rev_col_vector_t<T> * >
 inline  auto ordered_constrain(const T & x);

template <typename VarVec, require_var_col_vector_t<VarVec> * >
  auto ordered_constrain(const VarVec & x, scalar_type_t<VarVec> & lp);

template <typename T, require_rev_col_vector_t<T> * >
 inline  auto positive_ordered_constrain(const T & x);

template <typename T, require_rev_col_vector_t<T> * >
 inline  auto simplex_constrain(const T & y);

template <typename T, require_rev_col_vector_t<T> * >
  auto simplex_constrain(const T & y, scalar_type_t<T> & lp);

template <typename T, require_rev_matrix_t<T> * >
 inline  plain_type_t<T> stochastic_column_constrain(const T & y);

template <typename T, require_rev_matrix_t<T> * >
 inline  plain_type_t<T> stochastic_column_constrain(const T & y, scalar_type_t<T> & lp);

template <typename T, require_rev_matrix_t<T> * >
 inline  plain_type_t<T> stochastic_row_constrain(const T & y);

template <typename T, require_rev_matrix_t<T> * >
 inline  plain_type_t<T> stochastic_row_constrain(const T & y, scalar_type_t<T> & lp);

template <typename T, require_rev_col_vector_t<T> * >
 inline  auto sum_to_zero_constrain(T && y);

template <typename T, require_rev_col_vector_t<T> * >
 inline  auto sum_to_zero_constrain(T && y, scalar_type_t<T> & lp);

template <typename T, require_rev_col_vector_t<T> * >
 inline  auto unit_vector_constrain(const T & y);

template <typename T, require_eigen_col_vector_vt<is_var, T> * >
 inline  auto unit_vector_constrain(const T & y, var & lp);

template <typename T, require_var_col_vector_t<T> * >
 inline  auto unit_vector_constrain(const T & y, var & lp);

template <typename T>
 inline  fvar<T> abs(const fvar<T> & x);

template <typename T>
 inline  fvar<T> abs(const std::complex<fvar<T>> & z);

template <typename T, require_fvar_t<T> * >
 inline  auto sum(const std::vector<T> & m);

template <typename T, require_eigen_vt<is_fvar, T> * >
 inline  value_type_t<T> sum(const T & m);

template <typename T>
 inline  fvar<T> acos(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> acos(const std::complex<fvar<T>> & x);

template <typename T>
 inline  fvar<T> acosh(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> acosh(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> asin(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> asin(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> arg(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> asinh(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> asinh(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> atan(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> atan(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> atan2(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> atan2(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> atan2(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> atanh(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> atanh(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> bessel_first_kind(int v, const fvar<T> & z);

template <typename T>
 inline  fvar<T> bessel_second_kind(int v, const fvar<T> & z);

template <typename T>
 inline  fvar<T> beta(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> beta(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> beta(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> binary_log_loss(int y, const fvar<T> & y_hat);

template <typename T>
 inline  fvar<T> cbrt(const fvar<T> & x);

template <typename T>
 inline  fvar<T> ceil(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> conj(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> cos(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> cos(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> exp(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> exp(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> cosh(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> cosh(const std::complex<fvar<T>> & z);

template <typename EigMat, require_eigen_vt<is_fvar, EigMat> * >
 inline  value_type_t<EigMat> determinant(const EigMat & m);

template <typename T>
 inline  fvar<T> digamma(const fvar<T> & x);

template <typename T>
 inline  fvar<T> erf(const fvar<T> & x);

template <typename T>
 inline  fvar<T> erfc(const fvar<T> & x);

template <typename T>
 inline  fvar<T> exp2(const fvar<T> & x);

template <typename T>
 inline  fvar<T> expm1(const fvar<T> & x);

template <typename T>
 inline  fvar<T> fabs(const fvar<T> & x);

template <typename T>
 inline  fvar<T> falling_factorial(const fvar<T> & x, int n);

template <typename T>
 inline  fvar<T> fdim(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  fvar<T> fdim(const fvar<T> & x, double y);

template <typename T>
 inline  fvar<T> fdim(double x, const fvar<T> & y);

template <typename T>
 inline  fvar<T> floor(const fvar<T> & x);

template <typename T1, typename T2, typename T3, require_all_stan_scalar_t<T1, T2, T3> * >
 inline  fvar<return_type_t<T1, T2, T3>> fma(const fvar<T1> & x1, const fvar<T2> & x2, const fvar<T3> & x3);

template <typename T1, typename T2, typename T3, require_all_stan_scalar_t<T1, T2, T3> * >
 inline  fvar<return_type_t<T1, T2, T3>> fma(const T1 & x1, const fvar<T2> & x2, const fvar<T3> & x3);

template <typename T1, typename T2, typename T3, require_all_stan_scalar_t<T1, T2, T3> * >
 inline  fvar<return_type_t<T1, T2, T3>> fma(const fvar<T1> & x1, const T2 & x2, const fvar<T3> & x3);

template <typename T1, typename T2, typename T3, require_all_stan_scalar_t<T1, T2, T3> * >
 inline  fvar<return_type_t<T1, T2, T3>> fma(const fvar<T1> & x1, const fvar<T2> & x2, const T3 & x3);

template <typename T1, typename T2, typename T3, require_all_stan_scalar_t<T1, T2, T3> * >
 inline  fvar<return_type_t<T1, T2, T3>> fma(const T1 & x1, const T2 & x2, const fvar<T3> & x3);

template <typename T1, typename T2, typename T3, require_all_stan_scalar_t<T1, T2, T3> * >
 inline  fvar<return_type_t<T1, T2, T3>> fma(const fvar<T1> & x1, const T2 & x2, const T3 & x3);

template <typename T1, typename T2, typename T3, require_all_stan_scalar_t<T1, T2, T3> * >
 inline  fvar<return_type_t<T1, T2, T3>> fma(const T1 & x1, const fvar<T2> & x2, const T3 & x3);

template <typename T>
 inline  fvar<T> fmax(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> fmax(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> fmax(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> fmin(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> fmin(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> fmin(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> fmod(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> fmod(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> fmod(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> gamma_p(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> gamma_p(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> gamma_p(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> gamma_q(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> gamma_q(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> gamma_q(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> inv(const fvar<T> & x);

template <typename T>
 inline  fvar<T> inv_sqrt(const fvar<T> & x);

template <typename T>
 inline  fvar<T> inv_square(const fvar<T> & x);

template <typename T>
 inline  fvar<T> log(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> log(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> sqrt(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> sqrt(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> square(const fvar<T> & x);

template <typename T1, typename T2, require_any_fvar_t<base_type_t<T1>, base_type_t<T2>> * , require_all_stan_scalar_t<T1, T2> * >
 inline  auto pow(const T1 & x1, const T2 & x2);

template <typename T1, typename T2, require_any_container_t<T1, T2> * , require_all_not_matrix_st<is_var, T1, T2> * , require_any_fvar_t<base_type_t<T1>, base_type_t<T2>> * >
 inline  auto pow(const T1 & a, const T2 & b);

template <typename T>
 inline  fvar<T> inc_beta(const fvar<T> & a, const fvar<T> & b, const fvar<T> & x);

template <typename T>
 inline  fvar<T> inc_beta(double a, const fvar<T> & b, const fvar<T> & x);

template <typename T>
 inline  fvar<T> inc_beta(const fvar<T> & a, double b, const fvar<T> & x);

template <typename T>
 inline  fvar<T> inc_beta(const fvar<T> & a, const fvar<T> & b, double x);

template <typename T>
 inline  fvar<T> inc_beta(double a, double b, const fvar<T> & x);

template <typename T>
 inline  fvar<T> inc_beta(const fvar<T> & a, double b, double x);

template <typename T>
 inline  fvar<T> inc_beta(double a, const fvar<T> & b, double x);

template <typename T>
 inline  fvar<T> log1m(const fvar<T> & x);

template <typename T>
 inline  T value_of(const fvar<T> & v);

template <typename T>
  void grad_inc_beta(fvar<T> & g1, fvar<T> & g2, fvar<T> a, fvar<T> b, fvar<T> z);

template <typename Ta, typename Tz, typename FvarT, require_all_stan_scalar_t<Ta, Tz> * , require_any_fvar_t<Ta, Tz> * >
  FvarT hypergeometric_1f0(const Ta & a, const Tz & z);

template <typename Ta1, typename Ta2, typename Tb, typename Tz, require_all_stan_scalar_t<Ta1, Ta2, Tb, Tz> * , require_any_fvar_t<Ta1, Ta2, Tb, Tz> * >
 inline  return_type_t<Ta1, Ta1, Tb, Tz> hypergeometric_2F1(const Ta1 & a1, const Ta2 & a2, const Tb & b, const Tz & z);

template <typename Ta, typename Tb, typename Tz, typename FvarT, bool grad_a, bool grad_b, bool grad_z, require_all_vector_t<Ta, Tb> * , require_fvar_t<FvarT> * >
 inline  FvarT hypergeometric_pFq(const Ta & a, const Tb & b, const Tz & z);

template <typename T>
 inline  fvar<T> hypot(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> hypot(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> hypot(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> inv_erfc(const fvar<T> & x);

template <typename T>
 inline  fvar<T> inv_Phi(const fvar<T> & p);

template <typename T>
 inline  fvar<T> inv_cloglog(const fvar<T> & x);

template <typename T1, typename T2, typename T3, require_all_stan_scalar_t<T1, T2, T3> * , require_any_fvar_t<T1, T2, T3> * >
 inline  fvar<partials_return_t<T1, T2, T3>> inv_inc_beta(const T1 & a, const T2 & b, const T3 & p);

template <typename T>
 inline  fvar<T> inv_logit(const fvar<T> & x);

template <typename Mat1, typename Mat2, require_all_eigen_vt<is_fvar, Mat1, Mat2> * , require_vt_same<Mat1, Mat2> * , require_not_eigen_row_and_col_t<Mat1, Mat2> * >
 inline  auto multiply(const Mat1 & m1, const Mat2 & m2);

template <typename Mat1, typename Mat2, require_eigen_vt<is_fvar, Mat1> * , require_eigen_vt<std::is_floating_point, Mat2> * , require_not_eigen_row_and_col_t<Mat1, Mat2> * >
 inline  auto multiply(const Mat1 & m1, const Mat2 & m2);

template <typename Mat1, typename Mat2, require_eigen_vt<std::is_floating_point, Mat1> * , require_eigen_vt<is_fvar, Mat2> * , require_not_eigen_row_and_col_t<Mat1, Mat2> * >
 inline  auto multiply(const Mat1 & m1, const Mat2 & m2);

template <typename T, require_stan_scalar_t<T> * , require_not_fvar_t<T> * >
 inline  fvar<T> to_fvar(const T & x);

template <typename T, require_fvar_t<scalar_type_t<T>> * >
 inline  T && to_fvar(T && x);

template <typename T>
 inline  std::vector<fvar<T>> to_fvar(const std::vector<T> & v);

template <typename T>
 inline  std::vector<fvar<T>> to_fvar(const std::vector<T> & v, const std::vector<T> & d);

template <typename T, require_eigen_t<T> * , require_not_eigen_vt<is_fvar, T> * >
 inline  promote_scalar_t<fvar<value_type_t<T>>, T> to_fvar(const T & m);

template <typename T1, typename T2, require_all_eigen_t<T1, T2> * , require_vt_same<T1, T2> * >
 inline  promote_scalar_t<fvar<value_type_t<T1>>, T1> to_fvar(const T1 & val, const T2 & deriv);

template <typename EigMat, require_eigen_vt<is_fvar, EigMat> * >
 inline  Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime, EigMat::ColsAtCompileTime> inverse(const EigMat & m);

template <typename T>
 inline  int is_inf(const fvar<T> & x);

template <typename T, require_fvar_t<T> * >
 inline  bool is_nan(T && x);

template <typename T>
 inline  fvar<T> lambert_w0(const fvar<T> & x);

template <typename T>
 inline  fvar<T> lambert_wm1(const fvar<T> & x);

template <typename T>
 inline  fvar<T> lbeta(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> lbeta(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> lbeta(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> ldexp(const fvar<T> & a, int b);

template <typename T>
 inline  fvar<T> lgamma(const fvar<T> & x);

template <typename T>
 inline  fvar<return_type_t<T, int>> lmgamma(int x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> lmultiply(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> lmultiply(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> lmultiply(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> log10(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> log10(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> log1m_exp(const fvar<T> & x);

template <typename T>
 inline  fvar<T> log1m_inv_logit(const fvar<T> & x);

template <typename T>
 inline  fvar<T> log1p(const fvar<T> & x);

template <typename T>
 inline  fvar<T> log1p_exp(const fvar<T> & x);

template <typename T>
 inline  fvar<T> log2(const fvar<T> & x);

template <typename EigMat, require_eigen_vt<is_fvar, EigMat> * >
 inline  value_type_t<EigMat> log_determinant(const EigMat & m);

template <typename T>
 inline  fvar<T> log_diff_exp(const fvar<T> & x1, const fvar<T> & x2);

template <typename T1, typename T2, require_arithmetic_t<T1> * >
 inline  fvar<T2> log_diff_exp(const T1 & x1, const fvar<T2> & x2);

template <typename T1, typename T2, require_arithmetic_t<T2> * >
 inline  fvar<T1> log_diff_exp(const fvar<T1> & x1, const T2 & x2);

template <typename T>
 inline  fvar<T> log_falling_factorial(const fvar<T> & x, const fvar<T> & n);

template <typename T>
 inline  fvar<T> log_falling_factorial(double x, const fvar<T> & n);

template <typename T>
 inline  fvar<T> log_falling_factorial(const fvar<T> & x, double n);

template <typename T>
 inline  fvar<T> log_inv_logit(const fvar<T> & x);

template <typename T>
 inline  fvar<T> log_inv_logit_diff(const fvar<T> & x, const fvar<T> & y);

template <typename T>
 inline  fvar<T> log_inv_logit_diff(const fvar<T> & x, double y);

template <typename T>
 inline  fvar<T> log_inv_logit_diff(double x, const fvar<T> & y);

template <typename T_theta, typename T_lambda1, typename T_lambda2, int N>
 inline  void log_mix_partial_helper(const T_theta & theta, const T_lambda1 & lambda1, const T_lambda2 & lambda2, std::array<promote_args_t<T_theta, T_lambda1, T_lambda2>, N> & partials_array);

template <typename T>
 inline  fvar<T> log_mix(const fvar<T> & theta, const fvar<T> & lambda1, const fvar<T> & lambda2);

template <typename T, typename P, require_all_arithmetic_t<P> * >
 inline  fvar<T> log_mix(const fvar<T> & theta, const fvar<T> & lambda1, P lambda2);

template <typename T, typename P, require_all_arithmetic_t<P> * >
 inline  fvar<T> log_mix(const fvar<T> & theta, P lambda1, const fvar<T> & lambda2);

template <typename T, typename P, require_all_arithmetic_t<P> * >
 inline  fvar<T> log_mix(P theta, const fvar<T> & lambda1, const fvar<T> & lambda2);

template <typename T, typename P1, typename P2, require_all_arithmetic_t<P1, P2> * >
 inline  fvar<T> log_mix(const fvar<T> & theta, P1 lambda1, P2 lambda2);

template <typename T, typename P1, typename P2, require_all_arithmetic_t<P1, P2> * >
 inline  fvar<T> log_mix(P1 theta, const fvar<T> & lambda1, P2 lambda2);

template <typename T, typename P1, typename P2, require_all_arithmetic_t<P1, P2> * >
 inline  fvar<T> log_mix(P1 theta, P2 lambda1, const fvar<T> & lambda2);

template <typename T>
 inline  fvar<T> log_rising_factorial(const fvar<T> & x, const fvar<T> & n);

template <typename T>
 inline  fvar<T> log_rising_factorial(const fvar<T> & x, double n);

template <typename T>
 inline  fvar<T> log_rising_factorial(double x, const fvar<T> & n);

template <typename ColVec, require_eigen_col_vector_vt<is_fvar, ColVec> * >
 inline  auto softmax(const ColVec & alpha);

template <typename T, require_vector_st<is_fvar, T> * >
 inline  auto log_softmax(const T & x);

template <typename T>
 inline  fvar<T> log_sum_exp(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> log_sum_exp(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> log_sum_exp(const fvar<T> & x1, double x2);

template <typename T, require_container_st<is_fvar, T> * >
 inline  auto log_sum_exp(const T & x);

template <typename T>
 inline  fvar<T> logit(const fvar<T> & x);

template <typename T1, typename T2, require_all_eigen_vt<is_fvar, T1, T2> * , require_vt_same<T1, T2> * >
 inline  Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime, T2::ColsAtCompileTime> mdivide_left(const T1 & A, const T2 & b);

template <typename T1, typename T2, require_eigen_vt<std::is_arithmetic, T1> * , require_eigen_vt<is_fvar, T2> * >
 inline  Eigen::Matrix<value_type_t<T2>, T1::RowsAtCompileTime, T2::ColsAtCompileTime> mdivide_left(const T1 & A, const T2 & b);

template <typename T1, typename T2, require_eigen_vt<is_fvar, T1> * , require_eigen_vt<std::is_arithmetic, T2> * >
 inline  Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime, T2::ColsAtCompileTime> mdivide_left(const T1 & A, const T2 & b);

template <typename T, typename EigMat, require_eigen_vt<std::is_arithmetic, T> * , require_eigen_vt<is_fvar, EigMat> * >
 inline  Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, EigMat::ColsAtCompileTime> mdivide_left_ldlt(LDLT_factor<T> & A, const EigMat & b);

template <typename T1, typename T2, require_all_eigen_vt<is_fvar, T1, T2> * , require_vt_same<T1, T2> * >
 inline  Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime, T2::ColsAtCompileTime> mdivide_left_tri_low(const T1 & A, const T2 & b);

template <typename T1, typename T2, require_eigen_t<T1> * , require_vt_same<double, T1> * , require_eigen_vt<is_fvar, T2> * >
 inline  Eigen::Matrix<value_type_t<T2>, T1::RowsAtCompileTime, T2::ColsAtCompileTime> mdivide_left_tri_low(const T1 & A, const T2 & b);

template <typename T1, typename T2, require_eigen_vt<is_fvar, T1> * , require_eigen_t<T2> * , require_vt_same<double, T2> * >
 inline  Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime, T2::ColsAtCompileTime> mdivide_left_tri_low(const T1 & A, const T2 & b);

template <typename EigMat1, typename EigMat2, require_all_eigen_vt<is_fvar, EigMat1, EigMat2> * , require_vt_same<EigMat1, EigMat2> * >
 inline  Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime> mdivide_right(const EigMat1 & A, const EigMat2 & b);

template <typename EigMat1, typename EigMat2, require_eigen_vt<is_fvar, EigMat1> * , require_eigen_vt<std::is_arithmetic, EigMat2> * >
 inline  Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime> mdivide_right(const EigMat1 & A, const EigMat2 & b);

template <typename EigMat1, typename EigMat2, require_eigen_vt<std::is_arithmetic, EigMat1> * , require_eigen_vt<is_fvar, EigMat2> * >
 inline  Eigen::Matrix<value_type_t<EigMat2>, EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime> mdivide_right(const EigMat1 & A, const EigMat2 & b);

template <typename EigMat1, typename EigMat2, require_all_eigen_vt<is_fvar, EigMat1, EigMat2> * , require_vt_same<EigMat1, EigMat2> * >
 inline  Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime> mdivide_right_tri_low(const EigMat1 & A, const EigMat2 & b);

template <typename EigMat1, typename EigMat2, require_eigen_vt<is_fvar, EigMat1> * , require_eigen_vt<std::is_arithmetic, EigMat2> * >
 inline  Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime> mdivide_right_tri_low(const EigMat1 & A, const EigMat2 & b);

template <typename EigMat1, typename EigMat2, require_eigen_vt<std::is_arithmetic, EigMat1> * , require_eigen_vt<is_fvar, EigMat2> * >
 inline  Eigen::Matrix<value_type_t<EigMat2>, EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime> mdivide_right_tri_low(const EigMat1 & A, const EigMat2 & b);

template <typename T>
 inline  fvar<T> modified_bessel_first_kind(int v, const fvar<T> & z);

template <typename T>
 inline  fvar<T> modified_bessel_second_kind(int v, const fvar<T> & z);

template <typename T>
 inline  fvar<T> multiply_log(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> multiply_log(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> multiply_log(const fvar<T> & x1, double x2);

template <typename EigMat, require_eigen_vt<is_fvar, EigMat> * >
 inline  Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime, EigMat::RowsAtCompileTime> multiply_lower_tri_self_transpose(const EigMat & m);

template <typename T>
 inline  fvar<T> norm(const std::complex<fvar<T>> & z);

template <typename Container, require_eigen_vt<is_fvar, Container> * >
 inline  auto norm1(const Container & x);

template <typename Container, require_eigen_vt<is_fvar, Container> * >
 inline  auto norm2(const Container & x);

template <typename T>
 inline  fvar<T> owens_t(const fvar<T> & x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> owens_t(double x1, const fvar<T> & x2);

template <typename T>
 inline  fvar<T> owens_t(const fvar<T> & x1, double x2);

template <typename T>
 inline  fvar<T> Phi(const fvar<T> & x);

template <typename T>
 inline  fvar<T> Phi_approx(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> polar(const fvar<T> & r, const fvar<T> & theta);

template <typename T, typename U>
 inline  std::complex<fvar<T>> polar(const fvar<T> & r, U theta);

template <typename T, typename U>
 inline  std::complex<fvar<T>> polar(U r, const fvar<T> & theta);

template <typename T>
 inline  double primitive_value(const fvar<T> & v);

template <typename T>
 inline  std::complex<fvar<T>> proj(const std::complex<fvar<T>> & z);

template <typename EigMat1, typename EigMat2, require_all_eigen_t<EigMat1, EigMat2> * , require_not_eigen_col_vector_t<EigMat2> * , require_any_vt_fvar<EigMat1, EigMat2> * >
 inline  promote_scalar_t<return_type_t<EigMat1, EigMat2>, EigMat2> quad_form(const EigMat1 & A, const EigMat2 & B);

template <typename EigMat, typename ColVec, require_eigen_t<EigMat> * , require_eigen_col_vector_t<ColVec> * , require_any_vt_fvar<EigMat, ColVec> * >
 inline  return_type_t<EigMat, ColVec> quad_form(const EigMat & A, const ColVec & B);

template <typename EigMat1, typename EigMat2, require_all_eigen_t<EigMat1, EigMat2> * , require_not_eigen_col_vector_t<EigMat2> * , require_any_vt_fvar<EigMat1, EigMat2> * >
 inline  promote_scalar_t<return_type_t<EigMat1, EigMat2>, EigMat2> quad_form_sym(const EigMat1 & A, const EigMat2 & B);

template <typename EigMat, typename ColVec, require_eigen_t<EigMat> * , require_eigen_col_vector_t<ColVec> * , require_any_vt_fvar<EigMat, ColVec> * >
 inline  return_type_t<EigMat, ColVec> quad_form_sym(const EigMat & A, const ColVec & B);

template <typename T>
 inline  fvar<T> rising_factorial(const fvar<T> & x, int n);

template <typename T>
 inline  fvar<T> round(const fvar<T> & x);

template <typename T>
 inline  fvar<T> sin(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> sin(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> sinh(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> sinh(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> tan(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> tan(const std::complex<fvar<T>> & z);

template <typename T>
 inline  fvar<T> tanh(const fvar<T> & x);

template <typename T>
 inline  std::complex<fvar<T>> tanh(const std::complex<fvar<T>> & z);

template <typename EigMat, require_eigen_vt<is_fvar, EigMat> * >
 inline  Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime, EigMat::RowsAtCompileTime> tcrossprod(const EigMat & m);

template <typename T>
 inline  fvar<T> tgamma(const fvar<T> & x);

template <typename EigMat1, typename EigMat2, require_all_eigen_t<EigMat1, EigMat2> * , require_any_vt_fvar<EigMat1, EigMat2> * >
 inline  return_type_t<EigMat1, EigMat2> trace_quad_form(const EigMat1 & A, const EigMat2 & B);

template <typename T>
 inline  fvar<T> trigamma(const fvar<T> & u);

template <typename T>
 inline  fvar<T> trunc(const fvar<T> & x);

template <typename EigMat, require_eigen_col_vector_vt<is_fvar, EigMat> * >
 inline  auto unit_vector_constrain(const EigMat & y);

template <typename EigMat, typename T, require_eigen_vt<is_fvar, EigMat> * , require_stan_scalar_t<T> * >
 inline  auto unit_vector_constrain(const EigMat & y, T & lp);

template <typename T>
 inline  double value_of_rec(const fvar<T> & v);

template <typename T, typename F>
  void gradient(const F & f, const Eigen::Matrix<T, Eigen::Dynamic, 1> & x, T & fx, Eigen::Matrix<T, Eigen::Dynamic, 1> & grad_fx);

template <typename F, typename... TArgs, require_any_st_fvar<TArgs...> * >
 inline  auto finite_diff(const F & func, const TArgs &... args);

template <typename F, typename... TArgs, require_all_not_st_fvar<TArgs...> * >
 inline  auto finite_diff(const F & func, const TArgs &... args);

template <typename T, typename F>
  void hessian(const F & f, const Eigen::Matrix<T, Eigen::Dynamic, 1> & x, T & fx, Eigen::Matrix<T, Eigen::Dynamic, 1> & grad, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & H);

template <typename F, typename T_a, typename T_b, typename... Args, require_any_st_fvar<T_a, T_b, Args...> * >
 inline  return_type_t<T_a, T_b, Args...> integrate_1d_impl(const F & f, const T_a & a, const T_b & b, double relative_tolerance, std::ostream * msgs, const Args &... args);

template <typename F, typename T_a, typename T_b, typename T_theta, require_any_fvar_t<T_a, T_b, T_theta> * >
 inline  return_type_t<T_a, T_b, T_theta> integrate_1d(const F & f, const T_a & a, const T_b & b, const std::vector<T_theta> & theta, const std::vector<double> & x_r, const std::vector<int> & x_i, std::ostream * msgs, const double relative_tolerance);

template <typename T, typename F>
  void jacobian(const F & f, const Eigen::Matrix<T, Eigen::Dynamic, 1> & x, Eigen::Matrix<T, Eigen::Dynamic, 1> & fx, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & J);

template <typename T>
 inline  fvar<T> std_normal_log_qf(const fvar<T> & p);

// Classes Parsed: 60
// Functions Parsed: 891

} // namespace math
} // namespace stan
#endif // STAN_MATH_FORWARD_DECL_HPP

