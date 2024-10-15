
#ifndef STAN_MATH_FORWARD_DECL_HPP
#define STAN_MATH_FORWARD_DECL_HPP
#include <stan/math/manual_forward_decls.hpp>
namespace stan {
namespace math {

class stack_alloc;

template <typename ChainableT, typename ChainableAllocT>
struct AutodiffStackSingleton;

class chainable_alloc;

class vari_base;

template <typename T>
struct arena_allocator;

class chainable_alloc;

template <typename T>
class chainable_object;

template <typename T>
class unsafe_chainable_object;

template <typename MatrixType, typename>
class arena_matrix;

class vari_base;

template <typename T, typename>
class vari_view;

class gevv_vvv_vari;

template <typename EigRev, typename EigVari, typename EigDbl>
class vi_val_adj_functor;

template <typename EigRev, typename EigDbl>
class val_adj_functor;

template <typename EigVar, typename EigVari>
class vi_val_functor;

template <typename EigVar, typename EigVari>
class vi_adj_functor;

class ad_tape_observer;

class op_ddv_vari;

class op_dv_vari;

class op_dvd_vari;

class op_dvv_vari;

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

template <int Lmm, typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args>
class cvodes_integrator;

template <typename dae_type>
struct idas_service;

class idas_integrator;

template <typename F, typename Tyy, typename Typ, typename... T_par>
class dae_system;

template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args>
class cvodes_integrator_adjoint_vari;

// Classes Parsed: 56
template <typename T>
bool is_aligned(T *ptr, unsigned int bytes_aligned);

template <typename T>
auto make_chainable_ptr(T &&obj);

template <typename T>
auto make_unsafe_chainable_ptr(T &&obj);

static inline bool empty_nested();

static inline size_t nested_size();

static inline void grad();

template <typename Vari>
void grad(Vari *vi);

template <typename F>
inline void reverse_pass_callback(F &&functor);

template <typename Vari>
void grad(Vari *vi);

template <typename EigVar, typename EigVari, typename EigDbl>
inline void read_vi_val_adj(const EigVar &VarMat, EigVari &VariMat,
                            EigDbl &ValMat, EigDbl &AdjMat);

template <typename EigRev, typename EigDbl>
inline void read_val_adj(const EigRev &VarMat, EigDbl &ValMat, EigDbl &AdjMat);

template <typename EigVar, typename EigVari, typename EigDbl>
inline void read_vi_val(const EigVar &VarMat, EigVari &VariMat, EigDbl &ValMat);

template <typename EigVar, typename EigVari, typename EigDbl>
inline void read_vi_adj(const EigVar &VarMat, EigVari &VariMat, EigDbl &AdjMat);

template <typename T>
inline fvar<T> operator+(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> operator+(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> operator+(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> operator/(const fvar<T> &x1, const fvar<T> &x2);

template <typename T, typename U, require_arithmetic_t<U> * = nullptr>
inline fvar<T> operator/(const fvar<T> &x1, U x2);

template <typename T, typename U, require_arithmetic_t<U> * = nullptr>
inline fvar<T> operator/(U x1, const fvar<T> &x2);

template <typename T>
inline std::complex<fvar<T>> operator/(const std::complex<fvar<T>> &x1,
                                       const std::complex<fvar<T>> &x2);

template <typename T, typename U, require_arithmetic_t<U> * = nullptr>
inline std::complex<fvar<T>> operator/(const std::complex<fvar<T>> &x1,
                                       const std::complex<U> &x2);

template <typename T>
inline std::complex<fvar<T>> operator/(const std::complex<fvar<T>> &x1,
                                       const fvar<T> &x2);

template <typename T, typename U, require_arithmetic_t<U> * = nullptr>
inline std::complex<fvar<T>> operator/(const std::complex<fvar<T>> &x1, U x2);

template <typename T, typename U, require_arithmetic_t<U> * = nullptr>
inline std::complex<fvar<T>> operator/(const std::complex<U> &x1,
                                       const std::complex<fvar<T>> &x2);

template <typename T, typename U, require_arithmetic_t<U> * = nullptr>
inline std::complex<fvar<T>> operator/(const std::complex<U> &x1,
                                       const fvar<T> &x2);

template <typename T>
inline std::complex<fvar<T>> operator/(const fvar<T> &x1,
                                       const std::complex<fvar<T>> &x2);

template <typename T, typename U,
          std::enable_if_t<std::is_arithmetic<U>::value> * = nullptr>
inline std::complex<fvar<T>> operator/(const fvar<T> &x1,
                                       const std::complex<U> &x2);

template <typename T, typename U, require_arithmetic_t<U> * = nullptr>
inline std::complex<fvar<T>> operator/(U x1, const std::complex<fvar<T>> &x2);

template <typename T>
inline bool operator==(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline bool operator==(const fvar<T> &x, double y);

template <typename T>
inline bool operator==(double x, const fvar<T> &y);

template <typename T>
inline bool operator>(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline bool operator>(const fvar<T> &x, double y);

template <typename T>
inline bool operator>(double x, const fvar<T> &y);

template <typename T>
inline bool operator>=(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline bool operator>=(const fvar<T> &x, double y);

template <typename T>
inline bool operator>=(double x, const fvar<T> &y);

template <typename T>
inline bool operator<(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline bool operator<(double x, const fvar<T> &y);

template <typename T>
inline bool operator<(const fvar<T> &x, double y);

template <typename T>
inline bool operator<=(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline bool operator<=(const fvar<T> &x, double y);

template <typename T>
inline bool operator<=(double x, const fvar<T> &y);

template <typename T>
inline bool operator&&(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline bool operator&&(const fvar<T> &x, double y);

template <typename T>
inline bool operator&&(double x, const fvar<T> &y);

template <typename T>
inline bool operator||(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline bool operator||(const fvar<T> &x, double y);

template <typename T>
inline bool operator||(double x, const fvar<T> &y);

template <typename T>
inline fvar<T> operator*(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline fvar<T> operator*(double x, const fvar<T> &y);

template <typename T>
inline fvar<T> operator*(const fvar<T> &x, double y);

template <typename T>
inline std::complex<stan::math::fvar<T>> operator*(
    const std::complex<stan::math::fvar<T>> &x,
    const std::complex<stan::math::fvar<T>> &y);

template <typename T>
inline std::complex<stan::math::fvar<T>> operator*(
    const std::complex<double> &x, const std::complex<stan::math::fvar<T>> &y);

template <typename T>
inline std::complex<stan::math::fvar<T>> operator*(
    const std::complex<stan::math::fvar<T>> &x, const std::complex<double> &y);

template <typename T>
inline bool operator!=(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline bool operator!=(const fvar<T> &x, double y);

template <typename T>
inline bool operator!=(double x, const fvar<T> &y);

template <typename T>
inline fvar<T> operator-(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> operator-(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> operator-(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> operator-(const fvar<T> &x);

template <typename T>
inline bool operator!(const fvar<T> &x);

template <typename T>
inline fvar<T> operator+(const fvar<T> &x);

template <typename... Pargs>
inline double *accumulate_adjoints(double *dest, const stan::math::var &x,
                                   Pargs &&...args);

template <typename VarVec,
          require_std_vector_vt<stan::is_var, VarVec> * = nullptr,
          typename... Pargs>
inline double *accumulate_adjoints(double *dest, VarVec &&x, Pargs &&...args);

template <typename VecContainer,
          require_std_vector_st<stan::is_var, VecContainer> * = nullptr,
          require_std_vector_vt<stan::is_container, VecContainer> * = nullptr,
          typename... Pargs>
inline double *accumulate_adjoints(double *dest, VecContainer &&x,
                                   Pargs &&...args);

template <typename EigT, require_eigen_vt<stan::is_var, EigT> * = nullptr,
          typename... Pargs>
inline double *accumulate_adjoints(double *dest, EigT &&x, Pargs &&...args);

template <typename Arith, require_st_arithmetic<Arith> * = nullptr,
          typename... Pargs>
inline double *accumulate_adjoints(double *dest, Arith &&x, Pargs &&...args);

inline double *accumulate_adjoints(double *dest);

template <typename... Pargs>
inline double *accumulate_adjoints(double *dest, const stan::math::var &x,
                                   Pargs &&...args);

template <typename VarVec,
          require_std_vector_vt<stan::is_var, VarVec> * = nullptr,
          typename... Pargs>
inline double *accumulate_adjoints(double *dest, VarVec &&x, Pargs &&...args);

template <typename VecContainer,
          require_std_vector_st<stan::is_var, VecContainer> * = nullptr,
          require_std_vector_vt<stan::is_container, VecContainer> * = nullptr,
          typename... Pargs>
inline double *accumulate_adjoints(double *dest, VecContainer &&x,
                                   Pargs &&...args);

template <typename EigT, require_eigen_vt<stan::is_var, EigT> * = nullptr,
          typename... Pargs>
inline double *accumulate_adjoints(double *dest, EigT &&x, Pargs &&...args);

template <typename Arith, require_st_arithmetic<Arith> * = nullptr,
          typename... Pargs>
inline double *accumulate_adjoints(double *dest, Arith &&x, Pargs &&...args);

inline double *accumulate_adjoints(double *dest);

template <int R, int C>
stan::math::vari **build_vari_array(const Eigen::Matrix<var, R, C> &x);

template <typename... Pargs>
inline size_t count_vars(Pargs &&...args);

template <typename T, typename F>
internal::callback_vari<plain_type_t<T>, F> *make_callback_vari(T &&value,
                                                                F &&functor);

template <typename T, typename F>
var_value<plain_type_t<T>> make_callback_var(T &&value, F &&functor);

template <typename Arith,
          require_arithmetic_t<scalar_type_t<Arith>> * = nullptr>
inline Arith deep_copy_vars(Arith &&arg);

inline stan::math::var deep_copy_vars(const stan::math::var &arg);

template <typename VarVec,
          require_std_vector_vt<stan::is_var, VarVec> * = nullptr>
inline auto deep_copy_vars(VarVec &&arg);

template <typename VecContainer,
          require_std_vector_st<stan::is_var, VecContainer> * = nullptr,
          require_std_vector_vt<stan::is_container, VecContainer> * = nullptr>
inline auto deep_copy_vars(VecContainer &&arg);

template <typename EigT, require_eigen_vt<stan::is_var, EigT> * = nullptr>
inline auto deep_copy_vars(EigT &&arg);

static inline void recover_memory_nested();

static inline void set_zero_all_adjoints_nested();

static inline void start_nested();

inline stan::math::var operator+(const stan::math::var &a,
                                 const stan::math::var &b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline stan::math::var operator+(const stan::math::var &a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline stan::math::var operator+(Arith a, const stan::math::var &b);

template <typename VarMat1, typename VarMat2,
          require_all_rev_matrix_t<VarMat1, VarMat2> * = nullptr>
inline auto add(VarMat1 &&a, VarMat2 &&b);

template <typename Arith, typename VarMat,
          require_st_arithmetic<Arith> * = nullptr,
          require_rev_matrix_t<VarMat> * = nullptr>
inline auto add(VarMat &&a, const Arith &b);

template <typename Arith, typename VarMat,
          require_st_arithmetic<Arith> * = nullptr,
          require_rev_matrix_t<VarMat> * = nullptr>
inline auto add(const Arith &a, VarMat &&b);

template <typename Var, typename EigMat,
          require_var_vt<std::is_arithmetic, Var> * = nullptr,
          require_eigen_vt<std::is_arithmetic, EigMat> * = nullptr>
inline auto add(const Var &a, const EigMat &b);

template <typename EigMat, typename Var,
          require_eigen_vt<std::is_arithmetic, EigMat> * = nullptr,
          require_var_vt<std::is_arithmetic, Var> * = nullptr>
inline auto add(const EigMat &a, const Var &b);

template <typename Var, typename VarMat,
          require_var_vt<std::is_arithmetic, Var> * = nullptr,
          require_rev_matrix_t<VarMat> * = nullptr>
inline auto add(const Var &a, VarMat &&b);

template <typename Var, typename VarMat,
          require_var_vt<std::is_arithmetic, Var> * = nullptr,
          require_rev_matrix_t<VarMat> * = nullptr>
inline auto add(VarMat &&a, const Var &b);

template <typename T1, typename T2,
          require_any_var_vt<std::is_arithmetic, T1, T2> * = nullptr,
          require_any_arithmetic_t<T1, T2> * = nullptr>
inline auto add(const T1 &a, const T2 &b);

template <typename T1, typename T2,
          require_all_var_vt<std::is_arithmetic, T1, T2> * = nullptr>
inline auto add(const T1 &a, const T2 &b);

template <typename VarMat1, typename VarMat2,
          require_any_var_matrix_t<VarMat1, VarMat2> * = nullptr>
inline auto operator+(VarMat1 &&a, VarMat2 &&b);

inline stan::math::var operator-(const stan::math::var &a,
                                 const stan::math::var &b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline stan::math::var operator-(const stan::math::var &a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline stan::math::var operator-(Arith a, const stan::math::var &b);

template <typename VarMat1, typename VarMat2,
          require_all_rev_matrix_t<VarMat1, VarMat2> * = nullptr>
inline auto subtract(const VarMat1 &a, const VarMat2 &b);

template <typename Arith, typename VarMat,
          require_st_arithmetic<Arith> * = nullptr,
          require_rev_matrix_t<VarMat> * = nullptr>
inline auto subtract(const VarMat &a, const Arith &b);

template <typename Arith, typename VarMat,
          require_st_arithmetic<Arith> * = nullptr,
          require_rev_matrix_t<VarMat> * = nullptr>
inline auto subtract(const Arith &a, const VarMat &b);

template <typename Var, typename EigMat,
          require_var_vt<std::is_arithmetic, Var> * = nullptr,
          require_eigen_vt<std::is_arithmetic, EigMat> * = nullptr>
inline auto subtract(const Var &a, const EigMat &b);

template <typename EigMat, typename Var,
          require_eigen_vt<std::is_arithmetic, EigMat> * = nullptr,
          require_var_vt<std::is_arithmetic, Var> * = nullptr>
inline auto subtract(const EigMat &a, const Var &b);

template <typename Var, typename VarMat,
          require_var_vt<std::is_arithmetic, Var> * = nullptr,
          require_rev_matrix_t<VarMat> * = nullptr>
inline auto subtract(const Var &a, const VarMat &b);

template <typename Var, typename VarMat,
          require_rev_matrix_t<VarMat> * = nullptr,
          require_var_vt<std::is_arithmetic, Var> * = nullptr>
inline auto subtract(const VarMat &a, const Var &b);

template <typename T1, typename T2,
          require_any_var_vt<std::is_arithmetic, T1, T2> * = nullptr,
          require_any_arithmetic_t<T1, T2> * = nullptr>
inline auto subtract(const T1 &a, const T2 &b);

template <typename T1, typename T2,
          require_all_var_vt<std::is_arithmetic, T1, T2> * = nullptr>
inline auto subtract(const T1 &a, const T2 &b);

template <typename VarMat1, typename VarMat2,
          require_any_var_matrix_t<VarMat1, VarMat2> * = nullptr>
inline auto operator-(const VarMat1 &a, const VarMat2 &b);

inline stan::math::var operator*(const stan::math::var &a,
                                 const stan::math::var &b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline stan::math::var operator*(const stan::math::var &a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline stan::math::var operator*(Arith a, const stan::math::var &b);

inline std::complex<stan::math::var> operator*(
    const std::complex<stan::math::var> &x,
    const std::complex<stan::math::var> &y);

inline std::complex<stan::math::var> operator*(
    const std::complex<double> &x, const std::complex<stan::math::var> &y);

inline std::complex<stan::math::var> operator*(
    const std::complex<stan::math::var> &x, const std::complex<double> &y);

template <typename T, require_not_arena_t<T> * = nullptr,
          require_not_container_t<T> * = nullptr,
          require_not_matrix_cl_t<T> * = nullptr>
inline arena_t<T> to_arena(T &&a);

template <typename T, require_arena_t<T> * = nullptr,
          require_not_matrix_cl_t<T> * = nullptr,
          require_not_std_vector_t<T> * = nullptr>
inline std::remove_reference_t<T> to_arena(T &&a);

template <typename T, require_eigen_t<T> * = nullptr,
          require_not_arena_t<T> * = nullptr>
inline arena_t<T> to_arena(const T &a);

template <typename T>
inline std::vector<T, arena_allocator<T>> to_arena(
    const std::vector<T, arena_allocator<T>> &a);

template <typename T, require_same_t<T, arena_t<T>> * = nullptr>
inline arena_t<std::vector<T>> to_arena(const std::vector<T> &a);

template <typename T, require_not_same_t<T, arena_t<T>> * = nullptr>
inline arena_t<std::vector<T>> to_arena(const std::vector<T> &a);

template <bool Condition, typename T, std::enable_if_t<!Condition> * = nullptr>
inline T to_arena_if(T &&a);

template <bool Condition, typename T, std::enable_if_t<Condition> * = nullptr>
inline arena_t<T> to_arena_if(const T &a);

template <typename T>
inline auto &value_of(const var_value<T> &v);

inline stan::math::var operator/(const stan::math::var &dividend,
                                 const stan::math::var &divisor);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline stan::math::var operator/(const stan::math::var &dividend,
                                 Arith divisor);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline stan::math::var operator/(Arith dividend,
                                 const stan::math::var &divisor);

template <typename Scalar, typename Mat, require_matrix_t<Mat> * = nullptr,
          require_stan_scalar_t<Scalar> * = nullptr,
          require_all_st_var_or_arithmetic<Scalar, Mat> * = nullptr,
          require_any_st_var<Scalar, Mat> * = nullptr>
inline auto divide(const Mat &m, Scalar c);

template <typename Scalar, typename Mat, require_matrix_t<Mat> * = nullptr,
          require_stan_scalar_t<Scalar> * = nullptr,
          require_all_st_var_or_arithmetic<Scalar, Mat> * = nullptr,
          require_any_st_var<Scalar, Mat> * = nullptr>
inline auto divide(Scalar c, const Mat &m);

template <
    typename Mat1, typename Mat2,
    require_all_matrix_st<stan::is_var_or_arithmetic, Mat1, Mat2> * = nullptr,
    require_any_matrix_st<stan::is_var, Mat1, Mat2> * = nullptr>
inline auto divide(const Mat1 &m1, const Mat2 &m2);

template <typename T1, typename T2,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto operator/(const T1 &dividend, const T2 &divisor);

inline std::complex<var> operator/(const std::complex<var> &x1,
                                   const std::complex<var> &x2);

inline bool operator==(const stan::math::var &a, const stan::math::var &b);

inline bool operator==(int a, const stan::math::var &b);

inline bool operator==(const stan::math::var &a, int b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator==(const stan::math::var &a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator==(Arith a, const stan::math::var &b);

inline bool operator==(const stan::math::var &x, const std::complex<var> &z);

inline bool operator==(const std::complex<var> &z, const stan::math::var &y);

inline bool operator>(const stan::math::var &a, const stan::math::var &b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator>(const stan::math::var &a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator>(Arith a, const stan::math::var &b);

inline bool operator>=(const stan::math::var &a, const stan::math::var &b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator>=(const stan::math::var &a, Arith b);

template <typename Arith, typename Var, require_arithmetic_t<Arith> * = nullptr>
inline bool operator>=(Arith a, const stan::math::var &b);

inline bool operator<(const stan::math::var &a, const stan::math::var &b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator<(const stan::math::var &a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator<(Arith a, const stan::math::var &b);

inline bool operator<=(const stan::math::var &a, const stan::math::var &b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator<=(const stan::math::var &a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator<=(Arith a, const stan::math::var &b);

inline bool operator&&(const stan::math::var &x, const stan::math::var &y);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator&&(const stan::math::var &x, Arith y);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator&&(Arith x, const stan::math::var &y);

inline bool operator||(const stan::math::var &x, const stan::math::var &y);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator||(const stan::math::var &x, Arith y);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator||(Arith x, const stan::math::var &y);

inline bool operator!=(const stan::math::var &a, const stan::math::var &b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator!=(const stan::math::var &a, Arith b);

template <typename Arith, require_arithmetic_t<Arith> * = nullptr>
inline bool operator!=(Arith a, const stan::math::var &b);

inline bool operator!=(const stan::math::var &x, const std::complex<var> &z);

inline bool operator!=(const std::complex<var> &z, const stan::math::var &y);

inline stan::math::var &operator--(stan::math::var &a);

inline stan::math::var operator--(stan::math::var &a, int);

inline stan::math::var &operator++(stan::math::var &a);

inline stan::math::var operator++(stan::math::var &a, int);

inline stan::math::var operator-(const stan::math::var &a);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto operator-(const T &a);

inline bool operator!(const stan::math::var &x);

inline stan::math::var operator+(const stan::math::var &a);

template <typename T>
inline void dims(const var_value<T> &x, std::vector<int> &result);

template <typename T>
inline void dims(const vari_value<T> &x, std::vector<int> &result);

template <typename Arith, typename VecVar, typename VecArith,
          typename... ContainerOperands, typename... ContainerGradients>
inline stan::math::var precomputed_gradients(
    Arith value, const VecVar &operands, const VecArith &gradients,
    const std::tuple<ContainerOperands...> &container_operands,
    const std::tuple<ContainerGradients...> &container_gradients);

template <typename Arith, typename VecVar, typename VecArith>
inline stan::math::var precomputed_gradients(Arith value,
                                             const VecVar &operands,
                                             const VecArith &gradients);

inline void print_stack(std::ostream &o);

static inline void recover_memory();

static inline void set_zero_all_adjoints();

template <typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest,
                                     const stan::math::var &x, Pargs &&...args);

template <typename VarVec,
          require_std_vector_vt<stan::is_var, VarVec> * = nullptr,
          typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest, VarVec &&x,
                                     Pargs &&...args);

template <typename VecContainer,
          require_std_vector_st<stan::is_var, VecContainer> * = nullptr,
          require_std_vector_vt<stan::is_container, VecContainer> * = nullptr,
          typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest, VecContainer &&x,
                                     Pargs &&...args);

template <typename EigT, require_eigen_vt<stan::is_var, EigT> * = nullptr,
          typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest, EigT &&x,
                                     Pargs &&...args);

template <typename Arith, require_st_arithmetic<Arith> * = nullptr,
          typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest, Arith &&x,
                                     Pargs &&...args);

inline stan::math::vari **save_varis(stan::math::vari **dest);

template <typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest,
                                     const stan::math::var &x, Pargs &&...args);

template <typename VarVec,
          require_std_vector_vt<stan::is_var, VarVec> * = nullptr,
          typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest, VarVec &&x,
                                     Pargs &&...args);

template <typename VecContainer,
          require_std_vector_st<stan::is_var, VecContainer> * = nullptr,
          require_std_vector_vt<stan::is_container, VecContainer> * = nullptr,
          typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest, VecContainer &&x,
                                     Pargs &&...args);

template <typename EigT, require_eigen_vt<stan::is_var, EigT> * = nullptr,
          typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest, EigT &&x,
                                     Pargs &&...args);

template <typename Arith, require_st_arithmetic<Arith> * = nullptr,
          typename... Pargs>
inline stan::math::vari **save_varis(stan::math::vari **dest, Arith &&x,
                                     Pargs &&...args);

inline stan::math::vari **save_varis(stan::math::vari **dest);

inline void zero_adjoints() noexcept;

template <typename T, require_st_arithmetic<T> * = nullptr>
inline void zero_adjoints(T &x) noexcept;

inline void zero_adjoints(stan::math::var &x);

template <typename EigMat,
          require_eigen_vt<stan::is_autodiff, EigMat> * = nullptr>
inline void zero_adjoints(EigMat &x);

template <typename StdVec,
          require_std_vector_st<stan::is_autodiff, StdVec> * = nullptr>
inline void zero_adjoints(StdVec &x);

template <typename EigFvar, typename EigOut>
inline void read_fvar(const EigFvar &FvarMat, EigOut &ValMat, EigOut &DMat);

template <typename T, typename F>
void derivative(const F &f, const T &x, T &fx, T &dfx_dx);

template <typename F>
void hessian(const F &f, const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
             double &fx, Eigen::Matrix<double, Eigen::Dynamic, 1> &grad,
             Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &H);

template <typename F>
void finite_diff_grad_hessian_auto(const F &f, const Eigen::VectorXd &x,
                                   double &fx, Eigen::MatrixXd &hess,
                                   std::vector<Eigen::MatrixXd> &grad_hess_fx);

template <typename F>
void grad_hessian(
    const F &f, const Eigen::Matrix<double, Eigen::Dynamic, 1> &x, double &fx,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &H,
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> &grad_H);

template <typename T1, typename T2, typename F>
void gradient_dot_vector(const F &f,
                         const Eigen::Matrix<T1, Eigen::Dynamic, 1> &x,
                         const Eigen::Matrix<T2, Eigen::Dynamic, 1> &v, T1 &fx,
                         T1 &grad_fx_dot_v);

template <typename F>
void grad_tr_mat_times_hessian(
    const F &f, const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &M,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &grad_tr_MH);

template <typename F>
void hessian_times_vector(const F &f,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &v,
                          double &fx,
                          Eigen::Matrix<double, Eigen::Dynamic, 1> &Hv);

template <typename T, typename F>
void hessian_times_vector(const F &f,
                          const Eigen::Matrix<T, Eigen::Dynamic, 1> &x,
                          const Eigen::Matrix<T, Eigen::Dynamic, 1> &v, T &fx,
                          Eigen::Matrix<T, Eigen::Dynamic, 1> &Hv);

template <typename T, typename F>
void partial_derivative(const F &f,
                        const Eigen::Matrix<T, Eigen::Dynamic, 1> &x, int n,
                        T &fx, T &dfx_dxn);

template <typename Mat1, typename Mat2,
          require_all_eigen_vt<stan::is_fvar, Mat1, Mat2> * = nullptr,
          require_vt_same<Mat1, Mat2> * = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2> * = nullptr>
inline auto multiply(const Mat1 &m1, const Mat2 &m2);

template <typename Mat1, typename Mat2,
          require_eigen_vt<stan::is_fvar, Mat1> * = nullptr,
          require_eigen_vt<std::is_floating_point, Mat2> * = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2> * = nullptr>
inline auto multiply(const Mat1 &m1, const Mat2 &m2);

template <typename Mat1, typename Mat2,
          require_eigen_vt<std::is_floating_point, Mat1> * = nullptr,
          require_eigen_vt<stan::is_fvar, Mat2> * = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2> * = nullptr>
inline auto multiply(const Mat1 &m1, const Mat2 &m2);

template <typename EigMat, require_eigen_vt<stan::is_fvar, EigMat> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::RowsAtCompileTime>
tcrossprod(const EigMat &m);

template <typename T>
inline fvar<T> inv_sqrt(const fvar<T> &x);

template <typename T>
inline fvar<T> sqrt(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> sqrt(const std::complex<fvar<T>> &z);

template <typename EigMat,
          require_eigen_col_vector_vt<stan::is_fvar, EigMat> * = nullptr>
inline auto unit_vector_constrain(const EigMat &y);

template <typename EigMat, typename T,
          require_eigen_vt<stan::is_fvar, EigMat> * = nullptr,
          require_stan_scalar_t<T> * = nullptr>
inline auto unit_vector_constrain(const EigMat &y, T &lp);

template <typename T>
inline fvar<T> abs(const fvar<T> &x);

template <typename T>
inline fvar<T> abs(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> sum(const std::vector<fvar<T>> &m);

template <typename T, require_eigen_vt<stan::is_fvar, T> * = nullptr>
inline value_type_t<T> sum(const T &m);

template <typename T>
inline fvar<T> acos(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> acos(const std::complex<fvar<T>> &x);

template <typename T>
inline fvar<T> acosh(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> acosh(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> asin(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> asin(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> arg(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> asinh(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> asinh(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> atan(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> atan(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> atan2(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> atan2(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> atan2(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> atanh(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> atanh(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> bessel_first_kind(int v, const fvar<T> &z);

template <typename T>
inline fvar<T> bessel_second_kind(int v, const fvar<T> &z);

template <typename T>
inline fvar<T> beta(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> beta(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> beta(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> binary_log_loss(int y, const fvar<T> &y_hat);

template <typename T>
inline fvar<T> cbrt(const fvar<T> &x);

template <typename T>
inline fvar<T> ceil(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> conj(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> cos(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> cos(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> exp(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> exp(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> cosh(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> cosh(const std::complex<fvar<T>> &z);

template <typename EigMat, require_eigen_vt<stan::is_fvar, EigMat> * = nullptr>
inline value_type_t<EigMat> determinant(const EigMat &m);

template <typename T>
inline fvar<T> digamma(const fvar<T> &x);

template <typename T>
inline fvar<T> erf(const fvar<T> &x);

template <typename T>
inline fvar<T> erfc(const fvar<T> &x);

template <typename T>
inline fvar<T> exp2(const fvar<T> &x);

template <typename T>
inline fvar<T> expm1(const fvar<T> &x);

template <typename T>
inline fvar<T> fabs(const fvar<T> &x);

template <typename T>
inline fvar<T> falling_factorial(const fvar<T> &x, int n);

template <typename T>
inline fvar<T> fdim(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline fvar<T> fdim(const fvar<T> &x, double y);

template <typename T>
inline fvar<T> fdim(double x, const fvar<T> &y);

template <typename T>
inline fvar<T> floor(const fvar<T> &x);

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3> * = nullptr>
inline fvar<return_type_t<T1, T2, T3>> fma(const fvar<T1> &x1,
                                           const fvar<T2> &x2,
                                           const fvar<T3> &x3);

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3> * = nullptr>
inline fvar<return_type_t<T1, T2, T3>> fma(const T1 &x1, const fvar<T2> &x2,
                                           const fvar<T3> &x3);

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3> * = nullptr>
inline fvar<return_type_t<T1, T2, T3>> fma(const fvar<T1> &x1, const T2 &x2,
                                           const fvar<T3> &x3);

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3> * = nullptr>
inline fvar<return_type_t<T1, T2, T3>> fma(const fvar<T1> &x1,
                                           const fvar<T2> &x2, const T3 &x3);

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3> * = nullptr>
inline fvar<return_type_t<T1, T2, T3>> fma(const T1 &x1, const T2 &x2,
                                           const fvar<T3> &x3);

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3> * = nullptr>
inline fvar<return_type_t<T1, T2, T3>> fma(const fvar<T1> &x1, const T2 &x2,
                                           const T3 &x3);

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3> * = nullptr>
inline fvar<return_type_t<T1, T2, T3>> fma(const T1 &x1, const fvar<T2> &x2,
                                           const T3 &x3);

template <typename T>
inline fvar<T> fmax(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> fmax(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> fmax(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> fmin(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> fmin(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> fmin(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> fmod(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> fmod(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> fmod(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> gamma_p(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> gamma_p(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> gamma_p(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> gamma_q(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> gamma_q(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> gamma_q(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> inv(const fvar<T> &x);

template <typename T>
inline fvar<T> inv_square(const fvar<T> &x);

template <typename T>
inline fvar<T> pow(const fvar<T> &x1, const fvar<T> &x2);

template <typename T, typename U, require_arithmetic_t<U> * = nullptr>
inline fvar<T> pow(U x1, const fvar<T> &x2);

template <typename T, typename U, require_arithmetic_t<U> * = nullptr>
inline fvar<T> pow(const fvar<T> &x1, U x2);

template <typename V>
inline std::complex<fvar<V>> pow(const std::complex<fvar<V>> &x,
                                 const std::complex<fvar<V>> &y);

template <typename V, typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<fvar<V>> pow(const std::complex<fvar<V>> &x,
                                 const std::complex<T> &y);

template <typename V>
inline std::complex<fvar<V>> pow(const std::complex<fvar<V>> &x,
                                 const fvar<V> &y);

template <typename V, typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<fvar<V>> pow(const std::complex<fvar<V>> &x, const T &y);

template <typename V, typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<fvar<V>> pow(const std::complex<T> &x,
                                 const std::complex<fvar<V>> &y);

template <typename V, typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<fvar<V>> pow(const std::complex<T> &x, const fvar<V> &y);

template <typename V>
inline std::complex<fvar<V>> pow(const fvar<V> &x,
                                 const std::complex<fvar<V>> &y);

template <typename V, typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<fvar<V>> pow(const fvar<V> &x, const std::complex<T> &y);

template <typename T, typename V, require_arithmetic_t<T> * = nullptr>
inline std::complex<fvar<V>> pow(T x, const std::complex<fvar<V>> &y);

template <typename T>
inline std::complex<fvar<T>> pow(const std::complex<fvar<T>> &x, int y);

template <typename T>
inline fvar<T> inc_beta(const fvar<T> &a, const fvar<T> &b, const fvar<T> &x);

template <typename T>
inline fvar<T> inc_beta(double a, const fvar<T> &b, const fvar<T> &x);

template <typename T>
inline fvar<T> inc_beta(const fvar<T> &a, double b, const fvar<T> &x);

template <typename T>
inline fvar<T> inc_beta(const fvar<T> &a, const fvar<T> &b, double x);

template <typename T>
inline fvar<T> inc_beta(double a, double b, const fvar<T> &x);

template <typename T>
inline fvar<T> inc_beta(const fvar<T> &a, double b, double x);

template <typename T>
inline fvar<T> inc_beta(double a, const fvar<T> &b, double x);

template <typename T>
inline fvar<T> log(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> log(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> log1m(const fvar<T> &x);

template <typename T>
inline T value_of(const fvar<T> &v);

template <typename T>
void grad_inc_beta(fvar<T> &g1, fvar<T> &g2, fvar<T> a, fvar<T> b, fvar<T> z);

template <typename Ta, typename Tz, typename FvarT,
          require_all_stan_scalar_t<Ta, Tz> * = nullptr,
          require_any_fvar_t<Ta, Tz> * = nullptr>
FvarT hypergeometric_1f0(const Ta &a, const Tz &z);

template <typename Ta1, typename Ta2, typename Tb, typename Tz,
          require_all_stan_scalar_t<Ta1, Ta2, Tb, Tz> * = nullptr,
          require_any_fvar_t<Ta1, Ta2, Tb, Tz> * = nullptr>
inline return_type_t<Ta1, Ta1, Tb, Tz> hypergeometric_2F1(const Ta1 &a1,
                                                          const Ta2 &a2,
                                                          const Tb &b,
                                                          const Tz &z);

template <typename Ta, typename Tb, typename Tz, typename FvarT, bool grad_a,
          bool grad_b, bool grad_z, require_all_vector_t<Ta, Tb> * = nullptr,
          require_fvar_t<FvarT> * = nullptr>
inline FvarT hypergeometric_pFq(const Ta &a, const Tb &b, const Tz &z);

template <typename T>
inline fvar<T> hypot(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> hypot(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> hypot(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> square(const fvar<T> &x);

template <typename T>
inline fvar<T> inv_erfc(const fvar<T> &x);

template <typename T>
inline fvar<T> inv_Phi(const fvar<T> &p);

template <typename T>
inline fvar<T> inv_cloglog(const fvar<T> &x);

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3> * = nullptr,
          require_any_fvar_t<T1, T2, T3> * = nullptr>
inline fvar<partials_return_t<T1, T2, T3>> inv_inc_beta(const T1 &a,
                                                        const T2 &b,
                                                        const T3 &p);

template <typename T>
inline fvar<T> inv_logit(const fvar<T> &x);

template <typename T, require_stan_scalar_t<T> * = nullptr,
          require_not_fvar_t<T> * = nullptr>
inline fvar<T> to_fvar(const T &x);

template <typename T, require_fvar_t<scalar_type_t<T>> * = nullptr>
inline T &&to_fvar(T &&x);

template <typename T>
inline std::vector<fvar<T>> to_fvar(const std::vector<T> &v);

template <typename T>
inline std::vector<fvar<T>> to_fvar(const std::vector<T> &v,
                                    const std::vector<T> &d);

template <typename T, require_eigen_t<T> * = nullptr,
          require_not_eigen_vt<stan::is_fvar, T> * = nullptr>
inline promote_scalar_t<fvar<value_type_t<T>>, T> to_fvar(const T &m);

template <typename T1, typename T2, require_all_eigen_t<T1, T2> * = nullptr,
          require_vt_same<T1, T2> * = nullptr>
inline promote_scalar_t<fvar<value_type_t<T1>>, T1> to_fvar(const T1 &val,
                                                            const T2 &deriv);

template <typename EigMat, require_eigen_vt<stan::is_fvar, EigMat> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
inverse(const EigMat &m);

template <typename T>
inline int is_inf(const fvar<T> &x);

template <typename T, require_fvar_t<T> * = nullptr>
inline bool is_nan(T &&x);

template <typename T>
inline fvar<T> lambert_w0(const fvar<T> &x);

template <typename T>
inline fvar<T> lambert_wm1(const fvar<T> &x);

template <typename T>
inline fvar<T> lbeta(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> lbeta(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> lbeta(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> ldexp(const fvar<T> &a, int b);

template <typename T>
inline fvar<T> lgamma(const fvar<T> &x);

template <typename T>
inline fvar<return_type_t<T, int>> lmgamma(int x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> lmultiply(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> lmultiply(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> lmultiply(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> log10(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> log10(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> log1m_exp(const fvar<T> &x);

template <typename T>
inline fvar<T> log1m_inv_logit(const fvar<T> &x);

template <typename T>
inline fvar<T> log1p(const fvar<T> &x);

template <typename T>
inline fvar<T> log1p_exp(const fvar<T> &x);

template <typename T>
inline fvar<T> log2(const fvar<T> &x);

template <typename EigMat, require_eigen_vt<stan::is_fvar, EigMat> * = nullptr>
inline value_type_t<EigMat> log_determinant(const EigMat &m);

template <typename T>
inline fvar<T> log_diff_exp(const fvar<T> &x1, const fvar<T> &x2);

template <typename T1, typename T2, require_arithmetic_t<T1> * = nullptr>
inline fvar<T2> log_diff_exp(const T1 &x1, const fvar<T2> &x2);

template <typename T1, typename T2, require_arithmetic_t<T2> * = nullptr>
inline fvar<T1> log_diff_exp(const fvar<T1> &x1, const T2 &x2);

template <typename T>
inline fvar<T> log_falling_factorial(const fvar<T> &x, const fvar<T> &n);

template <typename T>
inline fvar<T> log_falling_factorial(double x, const fvar<T> &n);

template <typename T>
inline fvar<T> log_falling_factorial(const fvar<T> &x, double n);

template <typename T>
inline fvar<T> log_inv_logit(const fvar<T> &x);

template <typename T>
inline fvar<T> log_inv_logit_diff(const fvar<T> &x, const fvar<T> &y);

template <typename T>
inline fvar<T> log_inv_logit_diff(const fvar<T> &x, double y);

template <typename T>
inline fvar<T> log_inv_logit_diff(double x, const fvar<T> &y);

template <typename T_theta, typename T_lambda1, typename T_lambda2, int N>
inline void log_mix_partial_helper(
    const T_theta &theta, const T_lambda1 &lambda1, const T_lambda2 &lambda2,
    std::array<promote_args_t<T_theta, T_lambda1, T_lambda2>, N>
        &partials_array);

template <typename T>
inline fvar<T> log_mix(const fvar<T> &theta, const fvar<T> &lambda1,
                       const fvar<T> &lambda2);

template <typename T, typename P, require_all_arithmetic_t<P> * = nullptr>
inline fvar<T> log_mix(const fvar<T> &theta, const fvar<T> &lambda1, P lambda2);

template <typename T, typename P, require_all_arithmetic_t<P> * = nullptr>
inline fvar<T> log_mix(const fvar<T> &theta, P lambda1, const fvar<T> &lambda2);

template <typename T, typename P, require_all_arithmetic_t<P> * = nullptr>
inline fvar<T> log_mix(P theta, const fvar<T> &lambda1, const fvar<T> &lambda2);

template <typename T, typename P1, typename P2,
          require_all_arithmetic_t<P1, P2> * = nullptr>
inline fvar<T> log_mix(const fvar<T> &theta, P1 lambda1, P2 lambda2);

template <typename T, typename P1, typename P2,
          require_all_arithmetic_t<P1, P2> * = nullptr>
inline fvar<T> log_mix(P1 theta, const fvar<T> &lambda1, P2 lambda2);

template <typename T, typename P1, typename P2,
          require_all_arithmetic_t<P1, P2> * = nullptr>
inline fvar<T> log_mix(P1 theta, P2 lambda1, const fvar<T> &lambda2);

template <typename T>
inline fvar<T> log_rising_factorial(const fvar<T> &x, const fvar<T> &n);

template <typename T>
inline fvar<T> log_rising_factorial(const fvar<T> &x, double n);

template <typename T>
inline fvar<T> log_rising_factorial(double x, const fvar<T> &n);

template <typename ColVec,
          require_eigen_col_vector_vt<stan::is_fvar, ColVec> * = nullptr>
inline auto softmax(const ColVec &alpha);

template <typename T, require_vector_st<stan::is_fvar, T> * = nullptr>
inline auto log_softmax(const T &x);

template <typename T>
inline fvar<T> log_sum_exp(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> log_sum_exp(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> log_sum_exp(const fvar<T> &x1, double x2);

template <typename T, require_container_st<stan::is_fvar, T> * = nullptr>
inline auto log_sum_exp(const T &x);

template <typename T>
inline fvar<T> logit(const fvar<T> &x);

template <typename T1, typename T2,
          require_all_eigen_vt<stan::is_fvar, T1, T2> * = nullptr,
          require_vt_same<T1, T2> * = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left(const T1 &A, const T2 &b);

template <typename T1, typename T2,
          require_eigen_vt<std::is_arithmetic, T1> * = nullptr,
          require_eigen_vt<stan::is_fvar, T2> * = nullptr>
inline Eigen::Matrix<value_type_t<T2>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left(const T1 &A, const T2 &b);

template <typename T1, typename T2,
          require_eigen_vt<stan::is_fvar, T1> * = nullptr,
          require_eigen_vt<std::is_arithmetic, T2> * = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left(const T1 &A, const T2 &b);

template <typename T, typename EigMat,
          require_eigen_vt<std::is_arithmetic, T> * = nullptr,
          require_eigen_vt<stan::is_fvar, EigMat> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic,
                     EigMat::ColsAtCompileTime>
mdivide_left_ldlt(LDLT_factor<T> &A, const EigMat &b);

template <typename T1, typename T2,
          require_all_eigen_vt<stan::is_fvar, T1, T2> * = nullptr,
          require_vt_same<T1, T2> * = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1 &A, const T2 &b);

template <typename T1, typename T2, require_eigen_t<T1> * = nullptr,
          require_vt_same<double, T1> * = nullptr,
          require_eigen_vt<stan::is_fvar, T2> * = nullptr>
inline Eigen::Matrix<value_type_t<T2>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1 &A, const T2 &b);

template <
    typename T1, typename T2, require_eigen_vt<stan::is_fvar, T1> * = nullptr,
    require_eigen_t<T2> * = nullptr, require_vt_same<double, T2> * = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1 &A, const T2 &b);

template <typename EigMat1, typename EigMat2,
          require_all_eigen_vt<stan::is_fvar, EigMat1, EigMat2> * = nullptr,
          require_vt_same<EigMat1, EigMat2> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_right(const EigMat1 &A, const EigMat2 &b);

template <typename EigMat1, typename EigMat2,
          require_eigen_vt<stan::is_fvar, EigMat1> * = nullptr,
          require_eigen_vt<std::is_arithmetic, EigMat2> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_right(const EigMat1 &A, const EigMat2 &b);

template <typename EigMat1, typename EigMat2,
          require_eigen_vt<std::is_arithmetic, EigMat1> * = nullptr,
          require_eigen_vt<stan::is_fvar, EigMat2> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat2>, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_right(const EigMat1 &A, const EigMat2 &b);

template <typename EigMat1, typename EigMat2,
          require_all_eigen_vt<stan::is_fvar, EigMat1, EigMat2> * = nullptr,
          require_vt_same<EigMat1, EigMat2> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_right_tri_low(const EigMat1 &A, const EigMat2 &b);

template <typename EigMat1, typename EigMat2,
          require_eigen_vt<stan::is_fvar, EigMat1> * = nullptr,
          require_eigen_vt<std::is_arithmetic, EigMat2> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_right_tri_low(const EigMat1 &A, const EigMat2 &b);

template <typename EigMat1, typename EigMat2,
          require_eigen_vt<std::is_arithmetic, EigMat1> * = nullptr,
          require_eigen_vt<stan::is_fvar, EigMat2> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat2>, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_right_tri_low(const EigMat1 &A, const EigMat2 &b);

template <typename T>
inline fvar<T> modified_bessel_first_kind(int v, const fvar<T> &z);

template <typename T>
inline fvar<T> modified_bessel_second_kind(int v, const fvar<T> &z);

template <typename T>
inline fvar<T> multiply_log(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> multiply_log(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> multiply_log(const fvar<T> &x1, double x2);

template <typename EigMat, require_eigen_vt<stan::is_fvar, EigMat> * = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::RowsAtCompileTime>
multiply_lower_tri_self_transpose(const EigMat &m);

template <typename T>
inline fvar<T> norm(const std::complex<fvar<T>> &z);

template <typename Container,
          require_eigen_vt<stan::is_fvar, Container> * = nullptr>
inline auto norm1(const Container &x);

template <typename Container,
          require_eigen_vt<stan::is_fvar, Container> * = nullptr>
inline auto norm2(const Container &x);

template <typename T>
inline fvar<T> owens_t(const fvar<T> &x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> owens_t(double x1, const fvar<T> &x2);

template <typename T>
inline fvar<T> owens_t(const fvar<T> &x1, double x2);

template <typename T>
inline fvar<T> Phi(const fvar<T> &x);

template <typename T>
inline fvar<T> Phi_approx(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> polar(const fvar<T> &r, const fvar<T> &theta);

template <typename T, typename U>
inline std::complex<fvar<T>> polar(const fvar<T> &r, U theta);

template <typename T, typename U>
inline std::complex<fvar<T>> polar(U r, const fvar<T> &theta);

template <typename T>
inline double primitive_value(const fvar<T> &v);

template <typename T>
inline std::complex<fvar<T>> proj(const std::complex<fvar<T>> &z);

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2> * = nullptr,
          require_not_eigen_col_vector_t<EigMat2> * = nullptr,
          require_any_vt_fvar<EigMat1, EigMat2> * = nullptr>
inline promote_scalar_t<return_type_t<EigMat1, EigMat2>, EigMat2> quad_form(
    const EigMat1 &A, const EigMat2 &B);

template <typename EigMat, typename ColVec, require_eigen_t<EigMat> * = nullptr,
          require_eigen_col_vector_t<ColVec> * = nullptr,
          require_any_vt_fvar<EigMat, ColVec> * = nullptr>
inline return_type_t<EigMat, ColVec> quad_form(const EigMat &A,
                                               const ColVec &B);

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2> * = nullptr,
          require_not_eigen_col_vector_t<EigMat2> * = nullptr,
          require_any_vt_fvar<EigMat1, EigMat2> * = nullptr>
inline promote_scalar_t<return_type_t<EigMat1, EigMat2>, EigMat2> quad_form_sym(
    const EigMat1 &A, const EigMat2 &B);

template <typename EigMat, typename ColVec, require_eigen_t<EigMat> * = nullptr,
          require_eigen_col_vector_t<ColVec> * = nullptr,
          require_any_vt_fvar<EigMat, ColVec> * = nullptr>
inline return_type_t<EigMat, ColVec> quad_form_sym(const EigMat &A,
                                                   const ColVec &B);

template <typename T>
inline fvar<T> rising_factorial(const fvar<T> &x, int n);

template <typename T>
inline fvar<T> round(const fvar<T> &x);

template <typename T>
inline fvar<T> sin(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> sin(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> sinh(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> sinh(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> tan(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> tan(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> tanh(const fvar<T> &x);

template <typename T>
inline std::complex<fvar<T>> tanh(const std::complex<fvar<T>> &z);

template <typename T>
inline fvar<T> tgamma(const fvar<T> &x);

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2> * = nullptr,
          require_any_vt_fvar<EigMat1, EigMat2> * = nullptr>
inline return_type_t<EigMat1, EigMat2> trace_quad_form(const EigMat1 &A,
                                                       const EigMat2 &B);

template <typename T>
inline fvar<T> trigamma(const fvar<T> &u);

template <typename T>
inline fvar<T> trunc(const fvar<T> &x);

template <typename T>
inline double value_of_rec(const fvar<T> &v);

template <typename T, typename F>
void gradient(const F &f, const Eigen::Matrix<T, Eigen::Dynamic, 1> &x, T &fx,
              Eigen::Matrix<T, Eigen::Dynamic, 1> &grad_fx);

template <typename F, typename... TArgs,
          require_any_st_fvar<TArgs...> * = nullptr>
inline auto finite_diff(const F &func, const TArgs &...args);

template <typename F, typename... TArgs,
          require_all_not_st_fvar<TArgs...> * = nullptr>
inline auto finite_diff(const F &func, const TArgs &...args);

template <typename T, typename F>
void hessian(const F &f, const Eigen::Matrix<T, Eigen::Dynamic, 1> &x, T &fx,
             Eigen::Matrix<T, Eigen::Dynamic, 1> &grad,
             Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &H);

template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_fvar<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_impl(
    const F &f, const T_a &a, const T_b &b, double relative_tolerance,
    std::ostream *msgs, const Args &...args);

template <typename F, typename T_a, typename T_b, typename T_theta,
          require_any_fvar_t<T_a, T_b, T_theta> * = nullptr>
inline return_type_t<T_a, T_b, T_theta> integrate_1d(
    const F &f, const T_a &a, const T_b &b, const std::vector<T_theta> &theta,
    const std::vector<double> &x_r, const std::vector<int> &x_i,
    std::ostream *msgs, const double relative_tolerance);

template <typename T, typename F>
void jacobian(const F &f, const Eigen::Matrix<T, Eigen::Dynamic, 1> &x,
              Eigen::Matrix<T, Eigen::Dynamic, 1> &fx,
              Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &J);

template <typename T>
inline fvar<T> std_normal_log_qf(const fvar<T> &p);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> cholesky_corr_constrain(const T &y, int K);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> cholesky_corr_constrain(const T &y, int K,
                                                   scalar_type_t<T> &lp);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> cholesky_factor_constrain(const T &x, int M, int N);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> cholesky_factor_constrain(const T &x, int M, int N,
                                                     scalar_type_t<T> &lp);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> corr_matrix_constrain(const T &x, Eigen::Index k);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> corr_matrix_constrain(const T &x, Eigen::Index k,
                                                 scalar_type_t<T> &lp);

template <typename T1, typename T2, require_all_vector_t<T1, T2> * = nullptr,
          require_not_complex_t<return_type_t<T1, T2>> * = nullptr,
          require_all_not_std_vector_t<T1, T2> * = nullptr,
          require_any_st_var<T1, T2> * = nullptr>
inline stan::math::var dot_product(const T1 &v1, const T2 &v2);

template <typename T, require_eigen_vector_vt<stan::is_var, T> * = nullptr>
inline stan::math::var dot_self(const T &v);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline stan::math::var dot_self(const T &v);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline auto multiply_lower_tri_self_transpose(const T &L);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> cov_matrix_constrain(const T &x, Eigen::Index K);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> cov_matrix_constrain(const T &x, Eigen::Index K,
                                                scalar_type_t<T> &lp);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> cov_matrix_constrain_lkj(const T &x, size_t k);

template <typename T, require_var_vector_t<T> * = nullptr>
var_value<Eigen::MatrixXd> cov_matrix_constrain_lkj(const T &x, size_t k,
                                                    scalar_type_t<T> &lp);

template <typename T, require_eigen_vt<stan::is_var, T> * = nullptr>
inline var_value<
    Eigen::Matrix<double, T::RowsAtCompileTime, T::ColsAtCompileTime>>
to_var_value(const T &a);

template <typename T, require_var_t<T> * = nullptr>
inline T to_var_value(T &&a);

template <typename T>
inline auto to_var_value(const std::vector<T> &a);

template <typename T, typename... Types,
          require_eigen_vt<std::is_arithmetic, T> * = nullptr,
          require_any_var_matrix_t<T, Types...> * = nullptr>
inline auto identity_constrain(T &&x, Types &&...);

template <typename T, typename... Types,
          require_eigen_vt<stan::is_var, T> * = nullptr,
          require_any_var_matrix_t<T, Types...> * = nullptr>
inline auto identity_constrain(T &&x, Types &&...);

template <typename T, typename... Types, require_var_matrix_t<T> * = nullptr,
          require_any_var_matrix_t<T, Types...> * = nullptr>
inline auto identity_constrain(T &&x, Types &&...);

template <typename Alloc>
inline stan::math::var sum(const std::vector<var, Alloc> &m);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline stan::math::var sum(T &&x);

template <typename T, typename L, require_all_stan_scalar_t<T, L> * = nullptr,
          require_any_var_t<T, L> * = nullptr>
inline auto lb_constrain(const T &x, const L &lb);

template <typename T, typename L, require_all_stan_scalar_t<T, L> * = nullptr,
          require_any_var_t<T, L> * = nullptr>
inline auto lb_constrain(const T &x, const L &lb, stan::math::var &lp);

template <typename T, typename L, require_matrix_t<T> * = nullptr,
          require_stan_scalar_t<L> * = nullptr,
          require_any_st_var<T, L> * = nullptr>
inline auto lb_constrain(const T &x, const L &lb);

template <typename T, typename L, require_matrix_t<T> * = nullptr,
          require_stan_scalar_t<L> * = nullptr,
          require_any_st_var<T, L> * = nullptr>
inline auto lb_constrain(const T &x, const L &lb, return_type_t<T, L> &lp);

template <typename T, typename L, require_all_matrix_t<T, L> * = nullptr,
          require_any_st_var<T, L> * = nullptr>
inline auto lb_constrain(const T &x, const L &lb);

template <typename T, typename L, require_all_matrix_t<T, L> * = nullptr,
          require_any_st_var<T, L> * = nullptr>
inline auto lb_constrain(const T &x, const L &lb, return_type_t<T, L> &lp);

template <typename T, typename U, require_all_stan_scalar_t<T, U> * = nullptr,
          require_any_var_t<T, U> * = nullptr>
inline auto ub_constrain(const T &x, const U &ub);

template <typename T, typename U, require_all_stan_scalar_t<T, U> * = nullptr,
          require_any_var_t<T, U> * = nullptr>
inline auto ub_constrain(const T &x, const U &ub, return_type_t<T, U> &lp);

template <typename T, typename U, require_matrix_t<T> * = nullptr,
          require_stan_scalar_t<U> * = nullptr,
          require_any_st_var<T, U> * = nullptr>
inline auto ub_constrain(const T &x, const U &ub);

template <typename T, typename U, require_matrix_t<T> * = nullptr,
          require_stan_scalar_t<U> * = nullptr,
          require_any_st_var<T, U> * = nullptr>
inline auto ub_constrain(const T &x, const U &ub, return_type_t<T, U> &lp);

template <typename T, typename U, require_all_matrix_t<T, U> * = nullptr,
          require_any_st_var<T, U> * = nullptr>
inline auto ub_constrain(const T &x, const U &ub);

template <typename T, typename U, require_all_matrix_t<T, U> * = nullptr,
          require_any_st_var<T, U> * = nullptr>
inline auto ub_constrain(const T &x, const U &ub, return_type_t<T, U> &lp);

template <typename T, typename L, typename U,
          require_all_stan_scalar_t<T, L, U> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub);

template <typename T, typename L, typename U,
          require_all_stan_scalar_t<T, L, U> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub,
                          return_type_t<T, L, U> &lp);

template <typename T, typename L, typename U, require_matrix_t<T> * = nullptr,
          require_all_stan_scalar_t<L, U> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub);

template <typename T, typename L, typename U, require_matrix_t<T> * = nullptr,
          require_all_stan_scalar_t<L, U> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub,
                          return_type_t<T, L, U> &lp);

template <typename T, typename L, typename U,
          require_all_matrix_t<T, L> * = nullptr,
          require_stan_scalar_t<U> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub);

template <typename T, typename L, typename U,
          require_all_matrix_t<T, L> * = nullptr,
          require_stan_scalar_t<U> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub,
                          std::decay_t<return_type_t<T, L, U>> &lp);

template <typename T, typename L, typename U,
          require_all_matrix_t<T, U> * = nullptr,
          require_stan_scalar_t<L> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub);

template <typename T, typename L, typename U,
          require_all_matrix_t<T, U> * = nullptr,
          require_stan_scalar_t<L> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub,
                          std::decay_t<return_type_t<T, L, U>> &lp);

template <typename T, typename L, typename U,
          require_all_matrix_t<T, L, U> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub);

template <typename T, typename L, typename U,
          require_all_matrix_t<T, L, U> * = nullptr,
          require_var_t<return_type_t<T, L, U>> * = nullptr>
inline auto lub_constrain(const T &x, const L &lb, const U &ub,
                          return_type_t<T, L, U> &lp);

template <typename T, require_rev_col_vector_t<T> * = nullptr>
inline auto ordered_constrain(const T &x);

template <typename VarVec, require_var_col_vector_t<VarVec> * = nullptr>
auto ordered_constrain(const VarVec &x, scalar_type_t<VarVec> &lp);

template <typename T, require_rev_col_vector_t<T> * = nullptr>
inline auto positive_ordered_constrain(const T &x);

template <typename T, require_rev_col_vector_t<T> * = nullptr>
inline auto simplex_constrain(const T &y);

template <typename T, require_rev_col_vector_t<T> * = nullptr>
auto simplex_constrain(const T &y, scalar_type_t<T> &lp);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline plain_type_t<T> stochastic_column_constrain(const T &y);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline plain_type_t<T> stochastic_column_constrain(const T &y,
                                                   scalar_type_t<T> &lp);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline plain_type_t<T> stochastic_row_constrain(const T &y);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline plain_type_t<T> stochastic_row_constrain(const T &y,
                                                scalar_type_t<T> &lp);

template <typename T, require_rev_col_vector_t<T> * = nullptr>
inline auto sum_to_zero_constrain(T &&y);

template <typename T, require_rev_col_vector_t<T> * = nullptr>
inline auto sum_to_zero_constrain(T &&y, scalar_type_t<T> &lp);

template <typename T, require_rev_col_vector_t<T> * = nullptr>
inline auto unit_vector_constrain(const T &y);

template <typename T, require_eigen_col_vector_vt<stan::is_var, T> * = nullptr>
inline auto unit_vector_constrain(const T &y, stan::math::var &lp);

template <typename T, require_var_col_vector_t<T> * = nullptr>
inline auto unit_vector_constrain(const T &y, stan::math::var &lp);

inline stan::math::var Phi(const stan::math::var &a);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto Phi(const T &a);

inline stan::math::var Phi_approx(const stan::math::var &a);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto Phi_approx(const T &a);

inline stan::math::var fabs(const stan::math::var &a);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto fabs(const VarMat &a);

inline stan::math::var hypot(const stan::math::var &a,
                             const stan::math::var &b);

inline stan::math::var hypot(const stan::math::var &a, double b);

inline stan::math::var hypot(double a, const stan::math::var &b);

template <typename T>
inline auto abs(const var_value<T> &a);

inline stan::math::var abs(const std::complex<var> &z);

inline stan::math::var atan2(const stan::math::var &a,
                             const stan::math::var &b);

inline stan::math::var atan2(const stan::math::var &a, double b);

inline stan::math::var atan2(double a, const stan::math::var &b);

template <typename Mat1, typename Mat2,
          require_any_var_matrix_t<Mat1, Mat2> * = nullptr,
          require_all_matrix_t<Mat1, Mat2> * = nullptr>
inline auto atan2(const Mat1 &a, const Mat2 &b);

template <typename Scalar, typename VarMat,
          require_var_matrix_t<VarMat> * = nullptr,
          require_stan_scalar_t<Scalar> * = nullptr>
inline auto atan2(const Scalar &a, const VarMat &b);

template <typename VarMat, typename Scalar,
          require_var_matrix_t<VarMat> * = nullptr,
          require_stan_scalar_t<Scalar> * = nullptr>
inline auto atan2(const VarMat &a, const Scalar &b);

inline stan::math::var arg(const std::complex<var> &z);

template <typename T>
inline auto &value_of_rec(const var_value<T> &v);

inline int is_inf(const stan::math::var &v);

inline bool is_nan(const stan::math::var &v);

inline stan::math::var sinh(const stan::math::var &a);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto sinh(const VarMat &a);

inline std::complex<var> sinh(const std::complex<var> &z);

inline stan::math::var cos(stan::math::var a);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto cos(const VarMat &a);

inline std::complex<var> cos(const std::complex<var> &z);

inline stan::math::var sin(const stan::math::var &a);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto sin(const VarMat &a);

inline std::complex<var> sin(const std::complex<var> &z);

inline stan::math::var exp(const stan::math::var &a);

inline std::complex<var> exp(const std::complex<var> &z);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto exp(const T &x);

inline stan::math::var cosh(const stan::math::var &a);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto cosh(const VarMat &a);

inline std::complex<var> cosh(const std::complex<var> &z);

inline stan::math::var square(const stan::math::var &x);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto square(const T &x);

inline stan::math::var norm(const std::complex<var> &z);

inline stan::math::var sqrt(const stan::math::var &a);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto sqrt(const T &a);

inline std::complex<var> sqrt(const std::complex<var> &z);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto log(const var_value<T> &a);

inline std::complex<var> log(const std::complex<var> &z);

inline std::complex<var> polar(const stan::math::var &r,
                               const stan::math::var &theta);

template <typename T>
inline std::complex<var> polar(T r, const stan::math::var &theta);

template <typename T>
inline std::complex<var> polar(const stan::math::var &r, T theta);

inline stan::math::var asinh(const stan::math::var &x);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto asinh(const VarMat &x);

inline std::complex<var> asinh(const std::complex<var> &z);

inline stan::math::var asin(const stan::math::var &x);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto asin(const VarMat &x);

inline std::complex<var> asin(const std::complex<var> &z);

inline stan::math::var acos(const stan::math::var &x);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto acos(const VarMat &x);

inline std::complex<var> acos(const std::complex<var> &x);

inline stan::math::var acosh(const stan::math::var &x);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto acosh(const VarMat &x);

inline std::complex<var> acosh(const std::complex<var> &z);

template <typename T1, typename T2,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto append_col(const T1 &A, const T2 &B);

template <typename Scal, typename RowVec,
          require_stan_scalar_t<Scal> * = nullptr,
          require_t<is_eigen_row_vector<RowVec>> * = nullptr>
inline auto append_col(const Scal &A, const var_value<RowVec> &B);

template <typename RowVec, typename Scal,
          require_t<is_eigen_row_vector<RowVec>> * = nullptr,
          require_stan_scalar_t<Scal> * = nullptr>
inline auto append_col(const var_value<RowVec> &A, const Scal &B);

template <typename T1, typename T2,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto append_row(const T1 &A, const T2 &B);

template <typename Scal, typename ColVec,
          require_stan_scalar_t<Scal> * = nullptr,
          require_t<is_eigen_col_vector<ColVec>> * = nullptr>
inline auto append_row(const Scal &A, const var_value<ColVec> &B);

template <typename ColVec, typename Scal,
          require_t<is_eigen_col_vector<ColVec>> * = nullptr,
          require_stan_scalar_t<Scal> * = nullptr>
inline auto append_row(const var_value<ColVec> &A, const Scal &B);

inline int as_bool(const stan::math::var &v);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto as_array_or_scalar(T &&v);

template <typename T, require_var_col_vector_t<T> * = nullptr>
inline auto as_column_vector_or_scalar(T &&a);

template <typename T, require_var_row_vector_t<T> * = nullptr,
          require_not_var_col_vector_t<T> * = nullptr>
inline auto as_column_vector_or_scalar(T &&a);

inline stan::math::var atanh(const stan::math::var &x);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto atanh(const VarMat &x);

inline std::complex<var> atanh(const std::complex<var> &z);

inline stan::math::var atan(const stan::math::var &x);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto atan(const VarMat &x);

inline std::complex<var> atan(const std::complex<var> &z);

inline stan::math::var bessel_first_kind(int v, const stan::math::var &a);

template <typename T1, typename T2, require_st_integral<T1> * = nullptr,
          require_eigen_t<T2> * = nullptr>
inline auto bessel_first_kind(const T1 &v, const var_value<T2> &a);

inline stan::math::var bessel_second_kind(int v, const stan::math::var &a);

template <typename T1, typename T2, require_st_integral<T1> * = nullptr,
          require_eigen_t<T2> * = nullptr>
inline auto bessel_second_kind(const T1 &v, const var_value<T2> &a);

inline stan::math::var beta(const stan::math::var &a, const stan::math::var &b);

inline stan::math::var beta(const stan::math::var &a, double b);

inline stan::math::var beta(double a, const stan::math::var &b);

template <typename Mat1, typename Mat2,
          require_any_var_matrix_t<Mat1, Mat2> * = nullptr,
          require_all_matrix_t<Mat1, Mat2> * = nullptr>
inline auto beta(const Mat1 &a, const Mat2 &b);

template <typename Scalar, typename VarMat,
          require_var_matrix_t<VarMat> * = nullptr,
          require_stan_scalar_t<Scalar> * = nullptr>
inline auto beta(const Scalar &a, const VarMat &b);

template <typename VarMat, typename Scalar,
          require_var_matrix_t<VarMat> * = nullptr,
          require_stan_scalar_t<Scalar> * = nullptr>
inline auto beta(const VarMat &a, const Scalar &b);

inline stan::math::var binary_log_loss(int y, const stan::math::var &y_hat);

template <typename Mat, require_eigen_t<Mat> * = nullptr>
inline auto binary_log_loss(int y, const var_value<Mat> &y_hat);

template <typename StdVec, typename Mat, require_eigen_t<Mat> * = nullptr,
          require_st_integral<StdVec> * = nullptr>
inline auto binary_log_loss(const StdVec &y, const var_value<Mat> &y_hat);

inline stan::math::var cbrt(const stan::math::var &a);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto cbrt(const VarMat &a);

inline stan::math::var ceil(const stan::math::var &a);

template <typename T, require_matrix_t<T> * = nullptr>
inline auto ceil(const var_value<T> &a);

template <typename EigMat, require_eigen_vt<stan::is_var, EigMat> * = nullptr>
inline auto cholesky_decompose(const EigMat &A);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto cholesky_decompose(const T &A);

template <typename EigVec, require_rev_vector_t<EigVec> * = nullptr>
inline auto cumulative_sum(const EigVec &x);

template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2> * = nullptr,
          require_any_eigen_vt<stan::is_var, Mat1, Mat2> * = nullptr>
inline Eigen::Matrix<return_type_t<Mat1, Mat2>, 1, Mat1::ColsAtCompileTime>
columns_dot_product(const Mat1 &v1, const Mat2 &v2);

template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2> * = nullptr,
          require_any_var_matrix_t<Mat1, Mat2> * = nullptr>
inline auto columns_dot_product(const Mat1 &v1, const Mat2 &v2);

template <typename Mat, require_eigen_vt<stan::is_var, Mat> * = nullptr>
inline Eigen::Matrix<var, 1, Mat::ColsAtCompileTime> columns_dot_self(
    const Mat &x);

template <typename Mat, require_var_matrix_t<Mat> * = nullptr>
inline auto columns_dot_self(const Mat &x);

inline std::complex<var> conj(const std::complex<var> &z);

template <typename T, require_var_t<T> * = nullptr>
auto &adjoint_of(const T &x);

template <typename T, require_not_var_t<T> * = nullptr>
internal::nonexisting_adjoint adjoint_of(const T &x);

template <typename T_x, typename T_sigma,
          require_st_arithmetic<T_x> * = nullptr,
          require_stan_scalar_t<T_sigma> * = nullptr>
inline Eigen::Matrix<var, -1, -1> gp_exp_quad_cov(
    const std::vector<T_x> &x, const T_sigma sigma,
    const stan::math::var length_scale);

template <typename T_x,
          require_arithmetic_t<typename scalar_type<T_x>::type> * = nullptr>
inline Eigen::Matrix<var, -1, -1> cov_exp_quad(const std::vector<T_x> &x,
                                               const stan::math::var &sigma,
                                               const stan::math::var &l);

template <typename T_x,
          require_arithmetic_t<typename scalar_type<T_x>::type> * = nullptr>
inline Eigen::Matrix<var, -1, -1> cov_exp_quad(const std::vector<T_x> &x,
                                               double sigma,
                                               const stan::math::var &l);

template <int Options, typename VarMatrix, typename Vec1, typename Vec2,
          require_var_t<VarMatrix> * = nullptr,
          require_eigen_dense_base_t<value_type_t<VarMatrix>> * = nullptr,
          require_all_std_vector_vt<std::is_integral, Vec1, Vec2> * = nullptr>
inline auto to_soa_sparse_matrix(int m, int n, VarMatrix &&w, Vec1 &&u,
                                 Vec2 &&v);

template <int Options, typename MatrixVar, typename Vec1, typename Vec2,
          require_eigen_dense_base_vt<stan::is_var, MatrixVar> * = nullptr,
          require_all_std_vector_vt<std::is_integral, Vec1, Vec2> * = nullptr>
inline auto to_soa_sparse_matrix(int m, int n, MatrixVar &&w, Vec1 &&u,
                                 Vec2 &&v);

template <int Options, typename Mat, typename Vec1, typename Vec2,
          require_eigen_dense_base_vt<std::is_arithmetic, Mat> * = nullptr,
          require_all_std_vector_vt<std::is_integral, Vec1, Vec2> * = nullptr>
inline auto to_soa_sparse_matrix(int m, int n, Mat &&w, Vec1 &&u, Vec2 &&v);

template <typename T1, typename T2,
          require_any_rev_matrix_t<T1, T2> * = nullptr>
inline auto csr_matrix_times_vector(int m, int n, const T1 &w,
                                    const std::vector<int> &v,
                                    const std::vector<int> &u, const T2 &b);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline stan::math::var determinant(const T &m);

template <typename T1, typename T2, require_vector_t<T1> * = nullptr,
          require_matrix_t<T2> * = nullptr,
          require_any_st_var<T1, T2> * = nullptr>
auto diag_pre_multiply(const T1 &m1, const T2 &m2);

template <typename T1, typename T2, require_matrix_t<T1> * = nullptr,
          require_vector_t<T2> * = nullptr,
          require_any_st_var<T1, T2> * = nullptr>
auto diag_post_multiply(const T1 &m1, const T2 &m2);

inline stan::math::var digamma(const stan::math::var &a);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto digamma(const T &a);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline auto eigendecompose_sym(const T &m);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline auto eigenvalues_sym(const T &m);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline auto eigenvectors_sym(const T &m);

template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2> * = nullptr,
          require_any_rev_matrix_t<Mat1, Mat2> * = nullptr>
auto elt_divide(const Mat1 &m1, const Mat2 &m2);

template <typename Scal, typename Mat, require_stan_scalar_t<Scal> * = nullptr,
          require_var_matrix_t<Mat> * = nullptr>
auto elt_divide(Scal s, const Mat &m);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
          require_return_type_t<stan::is_var, T1, T2> * = nullptr,
          require_not_row_and_col_vector_t<T1, T2> * = nullptr>
inline auto multiply(T1 &&A, T2 &&B);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
          require_return_type_t<stan::is_var, T1, T2> * = nullptr,
          require_row_and_col_vector_t<T1, T2> * = nullptr>
inline stan::math::var multiply(const T1 &A, const T2 &B);

template <typename T1, typename T2, require_not_matrix_t<T1> * = nullptr,
          require_matrix_t<T2> * = nullptr,
          require_return_type_t<stan::is_var, T1, T2> * = nullptr,
          require_not_row_and_col_vector_t<T1, T2> * = nullptr>
inline auto multiply(const T1 &a, T2 &&B);

template <typename T1, typename T2, require_matrix_t<T1> * = nullptr,
          require_not_matrix_t<T2> * = nullptr,
          require_any_st_var<T1, T2> * = nullptr,
          require_not_complex_t<value_type_t<T1>> * = nullptr,
          require_not_complex_t<value_type_t<T2>> * = nullptr,
          require_not_row_and_col_vector_t<T1, T2> * = nullptr>
inline auto multiply(T1 &&A, T2 &&B);

template <typename T1, typename T2,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto operator*(T1 &&a, T2 &&b);

template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2> * = nullptr,
          require_any_rev_matrix_t<Mat1, Mat2> * = nullptr>
auto elt_multiply(const Mat1 &m1, const Mat2 &m2);

inline stan::math::var erf(const stan::math::var &a);

template <typename T, require_matrix_t<T> * = nullptr>
inline auto erf(const var_value<T> &a);

inline stan::math::var erfc(const stan::math::var &a);

template <typename T, require_matrix_t<T> * = nullptr>
inline auto erfc(const var_value<T> &a);

inline stan::math::var exp2(const stan::math::var &a);

template <typename T, require_eigen_t<T> * = nullptr>
inline auto exp2(const var_value<T> &a);

inline stan::math::var expm1(const stan::math::var &a);

template <typename T, require_eigen_t<T> * = nullptr>
inline auto expm1(const var_value<T> &a);

inline stan::math::var falling_factorial(const stan::math::var &a, int b);

template <typename T1, typename T2, require_eigen_t<T1> * = nullptr,
          require_st_integral<T2> * = nullptr>
inline auto falling_factorial(const var_value<T1> &a, const T2 &b);

inline stan::math::var fdim(const stan::math::var &a, const stan::math::var &b);

inline stan::math::var fdim(double a, const stan::math::var &b);

inline stan::math::var fdim(const stan::math::var &a, double b);

template <typename V, require_eigen_vector_vt<stan::is_complex, V> * = nullptr,
          require_var_t<base_type_t<value_type_t<V>>> * = nullptr>
inline plain_type_t<V> fft(const V &x);

template <typename V, require_eigen_vector_vt<stan::is_complex, V> * = nullptr,
          require_var_t<base_type_t<value_type_t<V>>> * = nullptr>
inline plain_type_t<V> inv_fft(const V &y);

template <typename M,
          require_eigen_dense_dynamic_vt<stan::is_complex, M> * = nullptr,
          require_var_t<base_type_t<value_type_t<M>>> * = nullptr>
inline plain_type_t<M> fft2(const M &x);

template <typename M,
          require_eigen_dense_dynamic_vt<stan::is_complex, M> * = nullptr,
          require_var_t<base_type_t<value_type_t<M>>> * = nullptr>
inline plain_type_t<M> inv_fft2(const M &y);

template <typename VarMat, typename S, require_var_matrix_t<VarMat> * = nullptr,
          require_var_t<S> * = nullptr>
inline void fill(VarMat &x, const S &y);

template <typename VarMat, typename S, require_var_matrix_t<VarMat> * = nullptr,
          require_arithmetic_t<S> * = nullptr>
inline void fill(VarMat &x, const S &y);

inline stan::math::var floor(const stan::math::var &a);

template <typename T, require_eigen_t<T> * = nullptr>
inline auto floor(const var_value<T> &a);

inline stan::math::var fma(const stan::math::var &x, const stan::math::var &y,
                           const stan::math::var &z);

template <typename Tc, require_arithmetic_t<Tc> * = nullptr>
inline stan::math::var fma(const stan::math::var &x, const stan::math::var &y,
                           Tc &&z);

template <typename Tb, require_arithmetic_t<Tb> * = nullptr>
inline stan::math::var fma(const stan::math::var &x, Tb &&y,
                           const stan::math::var &z);

template <typename Tb, typename Tc,
          require_all_arithmetic_t<Tb, Tc> * = nullptr>
inline stan::math::var fma(const stan::math::var &x, Tb &&y, Tc &&z);

template <typename Ta, typename Tc,
          require_all_arithmetic_t<Ta, Tc> * = nullptr>
inline stan::math::var fma(Ta &&x, const stan::math::var &y, Tc &&z);

template <typename Ta, typename Tb,
          require_all_arithmetic_t<Ta, Tb> * = nullptr>
inline stan::math::var fma(Ta &&x, Tb &&y, const stan::math::var &z);

template <typename Ta, require_arithmetic_t<Ta> * = nullptr>
inline stan::math::var fma(Ta &&x, const stan::math::var &y,
                           const stan::math::var &z);

template <typename T1, typename T2, typename T3,
          require_any_matrix_t<T1, T2, T3> * = nullptr,
          require_var_t<return_type_t<T1, T2, T3>> * = nullptr>
inline auto fma(const T1 &x, const T2 &y, const T3 &z);

inline stan::math::var fmax(const stan::math::var &a, const stan::math::var &b);

inline stan::math::var fmax(const stan::math::var &a, double b);

inline stan::math::var fmax(double a, const stan::math::var &b);

inline stan::math::var fmin(const stan::math::var &a, const stan::math::var &b);

inline stan::math::var fmin(const stan::math::var &a, double b);

inline stan::math::var fmin(double a, const stan::math::var &b);

inline stan::math::var fmod(const stan::math::var &a, const stan::math::var &b);

inline stan::math::var fmod(const stan::math::var &a, double b);

inline stan::math::var fmod(double a, const stan::math::var &b);

template <typename T, require_var_matrix_t<T> * = nullptr>
Eigen::Matrix<var, T::RowsAtCompileTime, T::ColsAtCompileTime> from_var_value(
    const T &a);

template <
    typename T,
    require_any_t<conjunction<is_eigen<T>, is_var<scalar_type_t<T>>>,
                  std::is_same<std::decay_t<T>, var>,
                  bool_constant<!std::is_same<scalar_type_t<T>, var>::value>>
        * = nullptr>
T from_var_value(T &&a);

template <typename T>
auto from_var_value(const std::vector<T> &a);

inline stan::math::var gamma_p(const stan::math::var &a,
                               const stan::math::var &b);

inline stan::math::var gamma_p(const stan::math::var &a, double b);

inline stan::math::var gamma_p(double a, const stan::math::var &b);

inline stan::math::var gamma_q(const stan::math::var &a,
                               const stan::math::var &b);

inline stan::math::var gamma_q(const stan::math::var &a, double b);

inline stan::math::var gamma_q(double a, const stan::math::var &b);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline auto inverse(const T &m);

template <typename VarMat, require_rev_matrix_t<VarMat> * = nullptr>
inline auto generalized_inverse(const VarMat &G);

template <typename T_x, typename T_sigma,
          require_st_arithmetic<T_x> * = nullptr,
          require_stan_scalar_t<T_sigma> * = nullptr>
inline Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> gp_periodic_cov(
    const std::vector<T_x> &x, const T_sigma sigma, const stan::math::var l,
    const stan::math::var p);

inline void grad(stan::math::var &v, Eigen::Matrix<var, Eigen::Dynamic, 1> &x,
                 Eigen::VectorXd &g);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto inv(const var_value<T> &a);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto inv_sqrt(const var_value<T> &a);

inline stan::math::var inv_square(const stan::math::var &a);

template <typename Scal1, typename Scal2,
          require_any_st_var<Scal1, Scal2> * = nullptr,
          require_all_stan_scalar_t<Scal1, Scal2> * = nullptr>
inline stan::math::var pow(const Scal1 &base, const Scal2 &exponent);

template <typename Mat1, typename Mat2,
          require_all_st_var_or_arithmetic<Mat1, Mat2> * = nullptr,
          require_any_matrix_st<stan::is_var, Mat1, Mat2> * = nullptr,
          require_all_not_stan_scalar_t<Mat1, Mat2> * = nullptr>
inline auto pow(const Mat1 &base, const Mat2 &exponent);

template <typename Mat1, typename Scal1,
          require_all_st_var_or_arithmetic<Mat1, Scal1> * = nullptr,
          require_all_matrix_st<stan::is_var, Mat1> * = nullptr,
          require_stan_scalar_t<Scal1> * = nullptr>
inline auto pow(const Mat1 &base, const Scal1 &exponent);

template <typename Scal1, typename Mat1,
          require_all_st_var_or_arithmetic<Scal1, Mat1> * = nullptr,
          require_stan_scalar_t<Scal1> * = nullptr,
          require_all_matrix_st<stan::is_var, Mat1> * = nullptr>
inline auto pow(Scal1 base, const Mat1 &exponent);

inline std::complex<var> pow(const std::complex<var> &x,
                             const std::complex<var> &y);

template <typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<var> pow(const std::complex<var> &x,
                             const std::complex<T> y);

inline std::complex<var> pow(const std::complex<var> &x,
                             const stan::math::var &y);

template <typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<var> pow(const std::complex<var> &x, T y);

template <typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<var> pow(std::complex<T> x, const std::complex<var> &y);

template <typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<var> pow(std::complex<T> x, const stan::math::var &y);

inline std::complex<var> pow(const stan::math::var &x,
                             const std::complex<var> &y);

template <typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<var> pow(const stan::math::var &x, std::complex<T> y);

template <typename T, require_arithmetic_t<T> * = nullptr>
inline std::complex<var> pow(T x, const std::complex<var> &y);

inline std::complex<var> pow(const std::complex<var> &x, int y);

inline stan::math::var inc_beta(const stan::math::var &a,
                                const stan::math::var &b,
                                const stan::math::var &c);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto lgamma(const var_value<T> &a);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto log1m(const var_value<T> &a);

template <typename Ta1, typename Ta2, typename Tb, typename Tz,
          require_all_stan_scalar_t<Ta1, Ta2, Tb, Tz> * = nullptr,
          require_any_var_t<Ta1, Ta2, Tb, Tz> * = nullptr>
inline return_type_t<Ta1, Ta1, Tb, Tz> hypergeometric_2F1(const Ta1 &a1,
                                                          const Ta2 &a2,
                                                          const Tb &b,
                                                          const Tz &z);

inline void grad_inc_beta(stan::math::var &g1, stan::math::var &g2,
                          const stan::math::var &a, const stan::math::var &b,
                          const stan::math::var &z);

template <typename Ta, typename Tz,
          require_all_stan_scalar_t<Ta, Tz> * = nullptr,
          require_any_var_t<Ta, Tz> * = nullptr>
stan::math::var hypergeometric_1f0(const Ta &a, const Tz &z);

template <typename Ta, typename Tb, typename Tz, bool grad_a, bool grad_b,
          bool grad_z, require_all_matrix_t<Ta, Tb> * = nullptr,
          require_return_type_t<stan::is_var, Ta, Tb, Tz> * = nullptr>
inline stan::math::var hypergeometric_pFq(const Ta &a, const Tb &b,
                                          const Tz &z);

inline stan::math::var if_else(bool c, const stan::math::var &y_true,
                               const stan::math::var &y_false);

inline stan::math::var if_else(bool c, double y_true,
                               const stan::math::var &y_false);

inline stan::math::var if_else(bool c, const stan::math::var &y_true,
                               double y_false);

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3> * = nullptr,
          require_any_var_t<T1, T2, T3> * = nullptr>
inline stan::math::var inv_inc_beta(const T1 &a, const T2 &b, const T3 &p);

template <typename VarMat, typename S, require_var_matrix_t<VarMat> * = nullptr,
          require_stan_scalar_t<S> * = nullptr>
inline void initialize_fill(VarMat &x, const S &y);

inline void initialize_variable(stan::math::var &variable,
                                const stan::math::var &value);

template <int R, int C>
inline void initialize_variable(Eigen::Matrix<var, R, C> &matrix,
                                const stan::math::var &value);

template <typename T>
inline void initialize_variable(std::vector<T> &variables,
                                const stan::math::var &value);

inline stan::math::var inv_Phi(const stan::math::var &p);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto inv_Phi(const T &p);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto inv_cloglog(const var_value<T> &a);

inline stan::math::var inv_erfc(const stan::math::var &a);

template <typename T, require_matrix_t<T> * = nullptr>
inline auto inv_erfc(const var_value<T> &a);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto inv_logit(const var_value<T> &a);

inline bool is_uninitialized(stan::math::var x);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto lambert_w0(const var_value<T> &a);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto lambert_wm1(const var_value<T> &a);

inline stan::math::var lbeta(const stan::math::var &a,
                             const stan::math::var &b);

inline stan::math::var lbeta(const stan::math::var &a, double b);

inline stan::math::var lbeta(double a, const stan::math::var &b);

inline stan::math::var ldexp(const stan::math::var &a, int b);

inline stan::math::var lmgamma(int a, const stan::math::var &b);

inline stan::math::var lmultiply(const stan::math::var &a,
                                 const stan::math::var &b);

inline stan::math::var lmultiply(const stan::math::var &a, double b);

inline stan::math::var lmultiply(double a, const stan::math::var &b);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto lmultiply(const T1 &a, const T2 &b);

template <typename T1, typename T2, require_var_matrix_t<T1> * = nullptr,
          require_stan_scalar_t<T2> * = nullptr>
inline auto lmultiply(const T1 &a, const T2 &b);

template <typename T1, typename T2, require_stan_scalar_t<T1> * = nullptr,
          require_var_matrix_t<T2> * = nullptr>
inline auto lmultiply(const T1 &a, const T2 &b);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto log10(const var_value<T> &a);

inline std::complex<var> log10(const std::complex<var> &z);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto log1m_exp(const var_value<T> &x);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto log1m_inv_logit(const var_value<T> &u);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto log1p(const var_value<T> &a);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto log1p_exp(const var_value<T> &a);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto log2(const var_value<T> &a);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline stan::math::var log_determinant(const T &m);

template <typename T, require_rev_matrix_t<T> * = nullptr>
stan::math::var log_determinant_ldlt(LDLT_factor<T> &A);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline stan::math::var log_determinant_spd(const T &M);

inline stan::math::var log_diff_exp(const stan::math::var &a,
                                    const stan::math::var &b);

inline stan::math::var log_diff_exp(const stan::math::var &a, double b);

inline stan::math::var log_diff_exp(double a, const stan::math::var &b);

inline stan::math::var log_falling_factorial(const stan::math::var &a,
                                             double b);

inline stan::math::var log_falling_factorial(const stan::math::var &a,
                                             const stan::math::var &b);

inline stan::math::var log_falling_factorial(double a,
                                             const stan::math::var &b);

template <typename T, require_arithmetic_t<T> * = nullptr>
inline auto log_inv_logit(const var_value<T> &u);

template <typename T, require_eigen_t<T> * = nullptr>
inline auto log_inv_logit(const var_value<T> &u);

inline stan::math::var log_inv_logit_diff(const stan::math::var &a, double b);

inline stan::math::var log_inv_logit_diff(const stan::math::var &a,
                                          const stan::math::var &b);

inline stan::math::var log_inv_logit_diff(double a, const stan::math::var &b);

inline void log_mix_partial_helper(
    double theta_val, double lambda1_val, double lambda2_val,
    double &one_m_exp_lam2_m_lam1, double &one_m_t_prod_exp_lam2_m_lam1,
    double &one_d_t_plus_one_m_t_prod_exp_lam2_m_lam1);

template <typename T_theta, typename T_lambda1, typename T_lambda2,
          require_any_var_t<T_theta, T_lambda1, T_lambda2> * = nullptr>
inline return_type_t<T_theta, T_lambda1, T_lambda2> log_mix(
    const T_theta &theta, const T_lambda1 &lambda1, const T_lambda2 &lambda2);

inline stan::math::var log_rising_factorial(const stan::math::var &a, double b);

inline stan::math::var log_rising_factorial(const stan::math::var &a,
                                            const stan::math::var &b);

inline stan::math::var log_rising_factorial(double a, const stan::math::var &b);

template <typename T, require_eigen_st<stan::is_var, T> * = nullptr>
auto log_softmax(const T &x);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto log_softmax(const T &x);

template <typename T, require_std_vector_st<stan::is_var, T> * = nullptr>
inline auto log_softmax(const T &x);

inline stan::math::var log_sum_exp(const stan::math::var &a,
                                   const stan::math::var &b);

inline stan::math::var log_sum_exp(const stan::math::var &a, double b);

inline stan::math::var log_sum_exp(double a, const stan::math::var &b);

template <typename T, require_eigen_st<stan::is_var, T> * = nullptr,
          require_not_var_matrix_t<T> * = nullptr>
inline stan::math::var log_sum_exp(const T &v);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline stan::math::var log_sum_exp(const T &x);

template <typename T, require_std_vector_st<stan::is_var, T> * = nullptr>
inline auto log_sum_exp(const T &x);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto logit(const var_value<T> &u);

template <typename Ta, typename Tb, require_all_eigen_t<Ta, Tb> * = nullptr,
          require_any_st_autodiff<Ta, Tb> * = nullptr>
inline Eigen::Matrix<return_type_t<Ta, Tb>, -1, Tb::ColsAtCompileTime>
matrix_exp_multiply(const Ta &A, const Tb &B);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline plain_type_t<T> matrix_power(const T &M, const int n);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
          require_any_st_var<T1, T2> * = nullptr>
inline auto mdivide_left(const T1 &A, const T2 &B);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
          require_any_st_var<T1, T2> * = nullptr>
inline auto mdivide_left_ldlt(LDLT_factor<T1> &A, const T2 &B);

template <typename EigMat1, typename EigMat2,
          require_all_eigen_matrix_base_vt<stan::is_var, EigMat1, EigMat2>
              * = nullptr>
inline Eigen::Matrix<var, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_left_spd(const EigMat1 &A, const EigMat2 &b);

template <typename EigMat1, typename EigMat2,
          require_eigen_matrix_base_vt<stan::is_var, EigMat1> * = nullptr,
          require_eigen_matrix_base_vt<std::is_arithmetic, EigMat2> * = nullptr>
inline Eigen::Matrix<var, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_left_spd(const EigMat1 &A, const EigMat2 &b);

template <typename EigMat1, typename EigMat2,
          require_eigen_matrix_base_vt<std::is_arithmetic, EigMat1> * = nullptr,
          require_eigen_matrix_base_vt<stan::is_var, EigMat2> * = nullptr>
inline Eigen::Matrix<var, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_left_spd(const EigMat1 &A, const EigMat2 &b);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto mdivide_left_spd(const T1 &A, const T2 &B);

template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_all_eigen_vt<stan::is_var, T1, T2> * = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left_tri(const T1 &A, const T2 &b);

template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_eigen_vt<std::is_arithmetic, T1> * = nullptr,
          require_eigen_vt<stan::is_var, T2> * = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left_tri(const T1 &A, const T2 &b);

template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_eigen_vt<stan::is_var, T1> * = nullptr,
          require_eigen_vt<std::is_arithmetic, T2> * = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left_tri(const T1 &A, const T2 &b);

template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_all_matrix_t<T1, T2> * = nullptr,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto mdivide_left_tri(const T1 &A, const T2 &B);

inline stan::math::var modified_bessel_first_kind(int v,
                                                  const stan::math::var &a);

inline stan::math::var modified_bessel_second_kind(int v,
                                                   const stan::math::var &a);

inline stan::math::var multiply_log(const stan::math::var &a,
                                    const stan::math::var &b);

inline stan::math::var multiply_log(const stan::math::var &a, double b);

inline stan::math::var multiply_log(double a, const stan::math::var &b);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto multiply_log(const T1 &a, const T2 &b);

template <typename T1, typename T2, require_var_matrix_t<T1> * = nullptr,
          require_stan_scalar_t<T2> * = nullptr>
inline auto multiply_log(const T1 &a, const T2 &b);

template <typename T1, typename T2, require_stan_scalar_t<T1> * = nullptr,
          require_var_matrix_t<T2> * = nullptr>
inline auto multiply_log(const T1 &a, const T2 &b);

template <typename T, require_eigen_vector_vt<stan::is_var, T> * = nullptr>
inline stan::math::var norm1(const T &v);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline stan::math::var norm1(const T &v);

template <typename T, require_eigen_vector_vt<stan::is_var, T> * = nullptr>
inline stan::math::var norm2(const T &v);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline stan::math::var norm2(const T &v);

template <typename Var1, typename Var2,
          require_all_st_var<Var1, Var2> * = nullptr,
          require_all_not_std_vector_t<Var1, Var2> * = nullptr>
inline auto owens_t(const Var1 &h, const Var2 &a);

template <typename Var, typename Arith,
          require_st_arithmetic<Arith> * = nullptr,
          require_all_not_std_vector_t<Var, Arith> * = nullptr,
          require_st_var<Var> * = nullptr>
inline auto owens_t(const Var &h, const Arith &a);

template <typename Arith, typename Var,
          require_st_arithmetic<Arith> * = nullptr,
          require_all_not_std_vector_t<Var, Arith> * = nullptr,
          require_st_var<Var> * = nullptr>
inline auto owens_t(const Arith &h, const Var &a);

inline double primitive_value(const stan::math::var &v);

inline std::complex<var> proj(const std::complex<var> &z);

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2> * = nullptr,
          require_not_eigen_col_vector_t<EigMat2> * = nullptr,
          require_any_vt_var<EigMat1, EigMat2> * = nullptr>
inline promote_scalar_t<stan::math::var, EigMat2> quad_form(const EigMat1 &A,
                                                            const EigMat2 &B,
                                                            bool symmetric);

template <typename EigMat, typename ColVec, require_eigen_t<EigMat> * = nullptr,
          require_eigen_col_vector_t<ColVec> * = nullptr,
          require_any_vt_var<EigMat, ColVec> * = nullptr>
inline stan::math::var quad_form(const EigMat &A, const ColVec &B,
                                 bool symmetric);

template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2> * = nullptr,
          require_not_col_vector_t<Mat2> * = nullptr,
          require_any_var_matrix_t<Mat1, Mat2> * = nullptr>
inline auto quad_form(const Mat1 &A, const Mat2 &B, bool symmetric);

template <typename Mat, typename Vec, require_matrix_t<Mat> * = nullptr,
          require_col_vector_t<Vec> * = nullptr,
          require_any_var_matrix_t<Mat, Vec> * = nullptr>
inline stan::math::var quad_form(const Mat &A, const Vec &B, bool symmetric);

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2> * = nullptr,
          require_any_vt_var<EigMat1, EigMat2> * = nullptr>
inline auto quad_form_sym(const EigMat1 &A, const EigMat2 &B);

template <typename T, require_var_vector_t<T> * = nullptr>
auto read_corr_L(const T &CPCs, size_t K);

template <typename T1, typename T2, require_var_vector_t<T1> * = nullptr,
          require_stan_scalar_t<T2> * = nullptr>
auto read_corr_L(const T1 &CPCs, size_t K, T2 &log_prob);

template <typename T_CPCs, require_var_vector_t<T_CPCs> * = nullptr>
inline var_value<Eigen::MatrixXd> read_corr_matrix(const T_CPCs &CPCs,
                                                   size_t K);

template <typename T_CPCs, require_var_vector_t<T_CPCs> * = nullptr>
inline var_value<Eigen::MatrixXd> read_corr_matrix(
    const T_CPCs &CPCs, size_t K, scalar_type_t<T_CPCs> &log_prob);

template <typename T_CPCs, typename T_sds,
          require_any_var_vector_t<T_CPCs, T_sds> * = nullptr,
          require_vt_same<T_CPCs, T_sds> * = nullptr>
inline auto read_cov_L(const T_CPCs &CPCs, const T_sds &sds,
                       scalar_type_t<T_CPCs> &log_prob);

template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2> * = nullptr,
          require_any_eigen_vt<stan::is_var, Mat1, Mat2> * = nullptr>
inline Eigen::Matrix<var, Mat1::RowsAtCompileTime, 1> rows_dot_product(
    const Mat1 &v1, const Mat2 &v2);

template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2> * = nullptr,
          require_any_var_matrix_t<Mat1, Mat2> * = nullptr>
inline auto rows_dot_product(const Mat1 &v1, const Mat2 &v2);

template <typename T_CPCs, typename T_sds,
          require_all_var_vector_t<T_CPCs, T_sds> * = nullptr>
var_value<Eigen::MatrixXd> read_cov_matrix(const T_CPCs &CPCs, const T_sds &sds,
                                           scalar_type_t<T_CPCs> &log_prob);

template <typename T_CPCs, typename T_sds,
          require_all_var_vector_t<T_CPCs, T_sds> * = nullptr>
inline var_value<Eigen::MatrixXd> read_cov_matrix(const T_CPCs &CPCs,
                                                  const T_sds &sds);

template <typename Ret, typename T, require_var_matrix_t<Ret> * = nullptr,
          require_var_t<T> * = nullptr>
inline auto rep_matrix(const T &x, int m, int n);

template <typename Ret, typename Vec, require_var_matrix_t<Ret> * = nullptr,
          require_var_matrix_t<Vec> * = nullptr>
inline auto rep_matrix(const Vec &x, int n);

template <typename T_ret, require_var_matrix_t<T_ret> * = nullptr,
          require_eigen_row_vector_t<value_type_t<T_ret>> * = nullptr>
inline auto rep_row_vector(stan::math::var x, int n);

template <typename T_ret, require_var_matrix_t<T_ret> * = nullptr,
          require_eigen_col_vector_t<value_type_t<T_ret>> * = nullptr>
inline auto rep_vector(stan::math::var x, int n);

inline stan::math::var rising_factorial(const stan::math::var &a, int b);

inline stan::math::var round(const stan::math::var &a);

template <typename Mat, require_eigen_vt<stan::is_var, Mat> * = nullptr>
inline Eigen::Matrix<var, Mat::RowsAtCompileTime, 1> rows_dot_self(
    const Mat &x);

template <typename Mat, require_var_matrix_t<Mat> * = nullptr>
inline auto rows_dot_self(const Mat &x);

template <typename T, require_eigen_st<stan::is_var, T> * = nullptr>
stan::math::var sd(const T &x);

template <typename T, require_var_matrix_t<T> * = nullptr>
stan::math::var sd(const T &x);

template <typename T, require_std_vector_st<stan::is_var, T> * = nullptr>
auto sd(const T &m);

template <typename EigMat, require_rev_matrix_t<EigMat> * = nullptr>
inline auto singular_values(const EigMat &m);

template <typename EigMat, require_rev_matrix_t<EigMat> * = nullptr>
inline auto svd(const EigMat &m);

template <typename EigMat, require_rev_matrix_t<EigMat> * = nullptr>
inline auto svd_U(const EigMat &m);

template <typename EigMat, require_rev_matrix_t<EigMat> * = nullptr>
inline auto svd_V(const EigMat &m);

template <typename Mat, require_rev_matrix_t<Mat> * = nullptr>
inline auto softmax(const Mat &alpha);

inline stan::math::var squared_distance(const stan::math::var &a,
                                        const stan::math::var &b);

inline stan::math::var squared_distance(const stan::math::var &a, double b);

inline stan::math::var squared_distance(double a, const stan::math::var &b);

template <typename EigVecVar1, typename EigVecVar2,
          require_all_eigen_vector_vt<stan::is_var, EigVecVar1, EigVecVar2>
              * = nullptr>
inline stan::math::var squared_distance(const EigVecVar1 &v1,
                                        const EigVecVar2 &v2);

template <typename EigVecVar, typename EigVecArith,
          require_eigen_vector_vt<stan::is_var, EigVecVar> * = nullptr,
          require_eigen_vector_vt<std::is_arithmetic, EigVecArith> * = nullptr>
inline stan::math::var squared_distance(const EigVecVar &v1,
                                        const EigVecArith &v2);

template <typename EigVecArith, typename EigVecVar,
          require_eigen_vector_vt<std::is_arithmetic, EigVecArith> * = nullptr,
          require_eigen_vector_vt<stan::is_var, EigVecVar> * = nullptr>
inline stan::math::var squared_distance(const EigVecArith &v1,
                                        const EigVecVar &v2);

template <typename T1, typename T2, require_all_vector_t<T1, T2> * = nullptr,
          require_any_var_vector_t<T1, T2> * = nullptr>
inline stan::math::var squared_distance(const T1 &A, const T2 &B);

inline void stan_print(std::ostream *o, const stan::math::var &x);

inline stan::math::var step(const stan::math::var &a);

inline stan::math::var tanh(const stan::math::var &a);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto tanh(const VarMat &a);

inline std::complex<var> tanh(const std::complex<var> &z);

inline stan::math::var tan(const stan::math::var &a);

template <typename VarMat, require_var_matrix_t<VarMat> * = nullptr>
inline auto tan(const VarMat &a);

inline std::complex<var> tan(const std::complex<var> &z);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline auto tcrossprod(const T &M);

inline stan::math::var tgamma(const stan::math::var &a);

template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto tgamma(const T &a);

inline stan::math::var to_var(double x);

inline stan::math::var &to_var(stan::math::var &x);

inline const stan::math::var &to_var(const stan::math::var &x);

inline std::vector<var> to_var(const std::vector<double> &v);

inline const std::vector<var> &to_var(const std::vector<var> &v);

inline std::vector<var> &to_var(std::vector<var> &v);

inline stan::math::matrix_v to_var(const stan::math::matrix_d &m);

inline stan::math::matrix_v &to_var(stan::math::matrix_v &m);

inline const stan::math::matrix_v &to_var(const stan::math::matrix_v &m);

inline stan::math::vector_v to_var(const stan::math::vector_d &v);

inline const stan::math::vector_v &to_var(const stan::math::vector_v &v);

inline stan::math::vector_v &to_var(stan::math::vector_v &v);

inline stan::math::row_vector_v to_var(const stan::math::row_vector_d &rv);

inline const stan::math::row_vector_v &to_var(
    const stan::math::row_vector_v &rv);

inline stan::math::row_vector_v &to_var(stan::math::row_vector_v &rv);

template <typename EigMat, require_eigen_t<EigMat> * = nullptr>
inline auto to_vector(const var_value<EigMat> &x);

template <typename T, require_rev_matrix_t<T> * = nullptr>
inline auto trace(const T &m);

template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
          require_any_st_var<T1, T2> * = nullptr>
inline stan::math::var trace_inv_quad_form_ldlt(LDLT_factor<T1> &A,
                                                const T2 &B);

template <typename Td, typename Ta, typename Tb,
          require_not_col_vector_t<Td> * = nullptr,
          require_all_matrix_t<Td, Ta, Tb> * = nullptr,
          require_any_st_var<Td, Ta, Tb> * = nullptr>
inline stan::math::var trace_gen_inv_quad_form_ldlt(const Td &D,
                                                    LDLT_factor<Ta> &A,
                                                    const Tb &B);

template <typename Td, typename Ta, typename Tb,
          require_col_vector_t<Td> * = nullptr,
          require_all_matrix_t<Ta, Tb> * = nullptr,
          require_any_st_var<Td, Ta, Tb> * = nullptr>
inline stan::math::var trace_gen_inv_quad_form_ldlt(const Td &D,
                                                    const LDLT_factor<Ta> &A,
                                                    const Tb &B);

template <typename Td, typename Ta, typename Tb,
          require_any_var_t<value_type_t<Td>, value_type_t<Ta>,
                            value_type_t<Tb>> * = nullptr,
          require_all_eigen_t<Td, Ta, Tb> * = nullptr>
inline stan::math::var trace_gen_quad_form(const Td &D, const Ta &A,
                                           const Tb &B);

template <typename Td, typename Ta, typename Tb,
          require_any_var_matrix_t<Td, Ta, Tb> * = nullptr,
          require_all_matrix_t<Td, Ta, Tb> * = nullptr>
inline stan::math::var trace_gen_quad_form(const Td &D, const Ta &A,
                                           const Tb &B);

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2> * = nullptr,
          require_any_st_var<EigMat1, EigMat2> * = nullptr>
inline return_type_t<EigMat1, EigMat2> trace_quad_form(const EigMat1 &A,
                                                       const EigMat2 &B);

template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2> * = nullptr,
          require_any_var_matrix_t<Mat1, Mat2> * = nullptr>
inline stan::math::var trace_quad_form(const Mat1 &A, const Mat2 &B);

inline stan::math::var trigamma(const stan::math::var &u);

inline stan::math::var trunc(const stan::math::var &a);

inline stan::math::var variance(const std::vector<var> &v);

template <typename EigMat, require_eigen_vt<stan::is_var, EigMat> * = nullptr>
stan::math::var variance(const EigMat &m);

template <typename Mat, require_var_matrix_t<Mat> * = nullptr>
inline stan::math::var variance(const Mat &x);

template <typename F>
void jacobian(const F &f, const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
              Eigen::Matrix<double, Eigen::Dynamic, 1> &fx,
              Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &J);

template <typename T1, typename T2>
void algebra_solver_check(const Eigen::Matrix<T1, Eigen::Dynamic, 1> &x,
                          const Eigen::Matrix<T2, Eigen::Dynamic, 1> y,
                          const std::vector<double> &dat,
                          const std::vector<int> &dat_int,
                          double function_tolerance, long max_num_steps);

template <typename F, typename T1, typename T2, typename T_u, typename T_f>
Eigen::Matrix<T2, -1, 1> algebra_solver_fp(
    const F &f, const Eigen::Matrix<T1, -1, 1> &x,
    const Eigen::Matrix<T2, -1, 1> &y, const std::vector<double> &x_r,
    const std::vector<int> &x_i, const std::vector<T_u> &u_scale,
    const std::vector<T_f> &f_scale, std::ostream *msgs, double f_tol,
    int max_num_steps);

template <typename F, typename T, typename... Args,
          require_eigen_vector_t<T> * = nullptr>
T &solve_powell_call_solver(const F &f, T &x, std::ostream *const msgs,
                            const double relative_tolerance,
                            const double function_tolerance,
                            const int64_t max_num_steps, const Args &...args);

template <typename F, typename T, typename... Args,
          require_eigen_vector_t<T> * = nullptr,
          require_all_st_arithmetic<Args...> * = nullptr>
Eigen::VectorXd solve_powell_tol(const F &f, const T &x,
                                 const double relative_tolerance,
                                 const double function_tolerance,
                                 const int64_t max_num_steps,
                                 std::ostream *const msgs, const Args &...args);

template <typename F, typename T, typename... T_Args,
          require_eigen_vector_t<T> * = nullptr>
Eigen::Matrix<stan::return_type_t<T_Args...>, Eigen::Dynamic, 1> solve_powell(
    const F &f, const T &x, std::ostream *const msgs, const T_Args &...args);

template <typename F, typename T1, typename T2,
          require_all_eigen_vector_t<T1, T2> * = nullptr>
Eigen::Matrix<value_type_t<T2>, Eigen::Dynamic, 1> algebra_solver(
    const F &f, const T1 &x, const T2 &y, const std::vector<double> &dat,
    const std::vector<int> &dat_int, std::ostream *msgs,
    const double relative_tolerance, const double function_tolerance,
    const int64_t max_num_steps);

template <typename F, typename T, typename... T_Args,
          require_eigen_vector_t<T> * = nullptr,
          require_any_st_var<T_Args...> * = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, 1> solve_powell_tol(
    const F &f, const T &x, const double relative_tolerance,
    const double function_tolerance, const int64_t max_num_steps,
    std::ostream *const msgs, const T_Args &...args);

template <typename F1, typename... Args>
Eigen::VectorXd kinsol_solve(const F1 &f, const Eigen::VectorXd &x,
                             const double scaling_step_tol,
                             const double function_tolerance,
                             const int64_t max_num_steps,
                             const bool custom_jacobian,
                             const int steps_eval_jacobian,
                             const int global_line_search,
                             std::ostream *const msgs, const Args &...args);

template <typename F, typename T, typename... Args,
          require_eigen_vector_t<T> * = nullptr,
          require_all_st_arithmetic<Args...> * = nullptr>
Eigen::VectorXd solve_newton_tol(const F &f, const T &x,
                                 const double scaling_step_size,
                                 const double function_tolerance,
                                 const int64_t max_num_steps,
                                 std::ostream *const msgs, const Args &...args);

template <typename F, typename T, typename... T_Args,
          require_eigen_vector_t<T> * = nullptr,
          require_any_st_var<T_Args...> * = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, 1> solve_newton_tol(
    const F &f, const T &x, const double scaling_step_size,
    const double function_tolerance, const int64_t max_num_steps,
    std::ostream *const msgs, const T_Args &...args);

template <typename F, typename T, typename... T_Args,
          require_eigen_vector_t<T> * = nullptr>
Eigen::Matrix<stan::return_type_t<T_Args...>, Eigen::Dynamic, 1> solve_newton(
    const F &f, const T &x, std::ostream *const msgs, const T_Args &...args);

template <typename F, typename T1, typename T2,
          require_all_eigen_vector_t<T1, T2> * = nullptr>
Eigen::Matrix<scalar_type_t<T2>, Eigen::Dynamic, 1> algebra_solver_newton(
    const F &f, const T1 &x, const T2 &y, const std::vector<double> &dat,
    const std::vector<int> &dat_int, std::ostream *const msgs,
    const double scaling_step_size, const double function_tolerance,
    const long max_num_steps);

template <typename T1, typename T2, typename F,
          require_any_var_matrix_t<T1, T2> * = nullptr,
          require_all_matrix_t<T1, T2> * = nullptr>
inline auto apply_scalar_binary(const T1 &x, const T2 &y, const F &f);

template <typename T1, typename T2, typename F,
          require_any_var_matrix_t<T1, T2> * = nullptr,
          require_any_std_vector_vt<std::is_integral, T1, T2> * = nullptr>
inline auto apply_scalar_binary(const T1 &x, const T2 &y, const F &f);

template <typename T1, typename T2, typename F,
          require_any_std_vector_vt<stan::is_std_vector, T1, T2> * = nullptr,
          require_any_std_vector_st<std::is_integral, T1, T2> * = nullptr,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto apply_scalar_binary(const T1 &x, const T2 &y, const F &f);

template <typename T1, typename T2, typename F,
          require_any_stan_scalar_t<T1, T2> * = nullptr,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto apply_scalar_binary(const T1 &x, const T2 &y, const F &f);

inline void cvodes_set_options(void *cvodes_mem, long max_num_steps);

template <typename F, typename T_y0_t0, typename T_t0, typename T_t,
          typename... Args,
          require_any_autodiff_t<T_y0_t0, T_t0, T_t, scalar_type_t<Args>...>
              * = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, 1> ode_store_sensitivities(
    const F &f, const std::vector<double> &coupled_state,
    const Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1> &y0, const T_t0 &t0,
    const T_t &t, std::ostream *msgs, const Args &...args);

template <typename F>
void gradient(const F &f, const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
              double &fx, Eigen::Matrix<double, Eigen::Dynamic, 1> &grad_fx);

template <typename F, typename EigVec, typename InputIt,
          require_eigen_vector_vt<std::is_arithmetic, EigVec> * = nullptr>
void gradient(const F &f, const EigVec &x, double &fx, InputIt first_grad_fx,
              InputIt last_grad_fx);

template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_var<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_impl(
    const F &f, const T_a &a, const T_b &b, double relative_tolerance,
    std::ostream *msgs, const Args &...args);

template <typename F, typename T_a, typename T_b, typename T_theta,
          require_any_var_t<T_a, T_b, T_theta> * = nullptr>
inline return_type_t<T_a, T_b, T_theta> integrate_1d(
    const F &f, const T_a &a, const T_b &b, const std::vector<T_theta> &theta,
    const std::vector<double> &x_r, const std::vector<int> &x_i,
    std::ostream *msgs, const double relative_tolerance);

template <typename F, typename T_yy, typename T_yp, typename... T_Args,
          require_all_eigen_col_vector_t<T_yy, T_yp> * = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
dae_tol_impl(const char *func, const F &f, const T_yy &yy0, const T_yp &yp0,
             double t0, const std::vector<double> &ts, double rtol, double atol,
             int64_t max_num_steps, std::ostream *msgs, const T_Args &...args);

template <typename F, typename T_yy, typename T_yp, typename... T_Args,
          require_all_eigen_col_vector_t<T_yy, T_yp> * = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
dae_tol(const F &f, const T_yy &yy0, const T_yp &yp0, double t0,
        const std::vector<double> &ts, double rtol, double atol,
        int64_t max_num_steps, std::ostream *msgs, const T_Args &...args);

template <typename F, typename T_yy, typename T_yp, typename... T_Args,
          require_all_eigen_col_vector_t<T_yy, T_yp> * = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
dae(const F &f, const T_yy &yy0, const T_yp &yp0, double t0,
    const std::vector<double> &ts, std::ostream *msgs, const T_Args &...args);

template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args, require_eigen_col_vector_t<T_y0> * = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_adams_tol_impl(const char *function_name, const F &f, const T_y0 &y0,
                   const T_t0 &t0, const std::vector<T_ts> &ts,
                   double relative_tolerance, double absolute_tolerance,
                   long max_num_steps, std::ostream *msgs,
                   const T_Args &...args);

template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args, require_eigen_col_vector_t<T_y0> * = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_adams_tol(const F &f, const T_y0 &y0, const T_t0 &t0,
              const std::vector<T_ts> &ts, double relative_tolerance,
              double absolute_tolerance, long max_num_steps, std::ostream *msgs,
              const T_Args &...args);

template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args, require_eigen_col_vector_t<T_y0> * = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_adams(const F &f, const T_y0 &y0, const T_t0 &t0,
          const std::vector<T_ts> &ts, std::ostream *msgs,
          const T_Args &...args);

template <typename F, typename T_y0, typename T_param, typename T_t0,
          typename T_ts>
std::vector<std::vector<return_type_t<T_y0, T_param, T_t0, T_ts>>>
integrate_ode_adams(const F &f, const std::vector<T_y0> &y0, const T_t0 &t0,
                    const std::vector<T_ts> &ts,
                    const std::vector<T_param> &theta,
                    const std::vector<double> &x, const std::vector<int> &x_int,
                    std::ostream *msgs, double relative_tolerance,
                    double absolute_tolerance, long max_num_steps);

template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args, require_eigen_col_vector_t<T_y0> * = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_bdf_tol_impl(const char *function_name, const F &f, const T_y0 &y0,
                 const T_t0 &t0, const std::vector<T_ts> &ts,
                 double relative_tolerance, double absolute_tolerance,
                 long max_num_steps, std::ostream *msgs, const T_Args &...args);

template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args, require_eigen_col_vector_t<T_y0> * = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_bdf_tol(const F &f, const T_y0 &y0, const T_t0 &t0,
            const std::vector<T_ts> &ts, double relative_tolerance,
            double absolute_tolerance, long max_num_steps, std::ostream *msgs,
            const T_Args &...args);

template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args, require_eigen_col_vector_t<T_y0> * = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_bdf(const F &f, const T_y0 &y0, const T_t0 &t0, const std::vector<T_ts> &ts,
        std::ostream *msgs, const T_Args &...args);

template <typename F, typename T_y0, typename T_param, typename T_t0,
          typename T_ts>
std::vector<std::vector<return_type_t<T_y0, T_param, T_t0, T_ts>>>
integrate_ode_bdf(const F &f, const std::vector<T_y0> &y0, const T_t0 &t0,
                  const std::vector<T_ts> &ts,
                  const std::vector<T_param> &theta,
                  const std::vector<double> &x, const std::vector<int> &x_int,
                  std::ostream *msgs, double relative_tolerance,
                  double absolute_tolerance, long max_num_steps);

template <
    typename F, typename T_y0, typename T_t0, typename T_ts,
    typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args,
    require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd, T_abs_tol_bwd>
        * = nullptr,
    require_any_not_st_arithmetic<T_y0, T_t0, T_ts, T_Args...> * = nullptr>
auto ode_adjoint_impl(const char *function_name, F &&f, const T_y0 &y0,
                      const T_t0 &t0, const std::vector<T_ts> &ts,
                      double relative_tolerance_forward,
                      const T_abs_tol_fwd &absolute_tolerance_forward,
                      double relative_tolerance_backward,
                      const T_abs_tol_bwd &absolute_tolerance_backward,
                      double relative_tolerance_quadrature,
                      double absolute_tolerance_quadrature, long max_num_steps,
                      long num_steps_between_checkpoints,
                      int interpolation_polynomial, int solver_forward,
                      int solver_backward, std::ostream *msgs,
                      const T_Args &...args);

template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args,
          require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd, T_abs_tol_bwd>
              * = nullptr,
          require_all_st_arithmetic<T_y0, T_t0, T_ts, T_Args...> * = nullptr>
std::vector<Eigen::VectorXd> ode_adjoint_impl(
    const char *function_name, F &&f, const T_y0 &y0, const T_t0 &t0,
    const std::vector<T_ts> &ts, double relative_tolerance_forward,
    const T_abs_tol_fwd &absolute_tolerance_forward,
    double relative_tolerance_backward,
    const T_abs_tol_bwd &absolute_tolerance_backward,
    double relative_tolerance_quadrature, double absolute_tolerance_quadrature,
    long max_num_steps, long num_steps_between_checkpoints,
    int interpolation_polynomial, int solver_forward, int solver_backward,
    std::ostream *msgs, const T_Args &...args);

template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args,
          require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd, T_abs_tol_bwd>
              * = nullptr>
auto ode_adjoint_tol_ctl(F &&f, const T_y0 &y0, const T_t0 &t0,
                         const std::vector<T_ts> &ts,
                         double relative_tolerance_forward,
                         const T_abs_tol_fwd &absolute_tolerance_forward,
                         double relative_tolerance_backward,
                         const T_abs_tol_bwd &absolute_tolerance_backward,
                         double relative_tolerance_quadrature,
                         double absolute_tolerance_quadrature,
                         long max_num_steps, long num_steps_between_checkpoints,
                         int interpolation_polynomial, int solver_forward,
                         int solver_backward, std::ostream *msgs,
                         const T_Args &...args);

template <typename T, require_stan_scalar_or_eigen_t<T> * = nullptr>
inline auto std_normal_log_qf(const var_value<T> &log_p);

// Functions Parsed: 913

}  // namespace math
}  // namespace stan
#endif  // STAN_MATH_FORWARD_DECL_HPP
