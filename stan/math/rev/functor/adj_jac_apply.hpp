#ifndef STAN_MATH_REV_FUNCTOR_ADJ_JAC_APPLY_HPP
#define STAN_MATH_REV_FUNCTOR_ADJ_JAC_APPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/meta/conditional_sequence.hpp>
#include <stan/math/rev/meta/var_tuple_filter.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/count_vars.hpp>
#include <stan/math/rev/core/save_varis.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

class AdjJacOp {
public:
  template <typename Derived>
  inline double* allocate_and_save(const Eigen::EigenBase<Derived>& x) const {
    double* x_mem_ = allocate(x);

    Eigen::Map<Eigen::MatrixXd>(x_mem_, x.rows(), x.cols()) = x;

    return x_mem_;
  }

  template <typename Derived>
  inline double* allocate(const Eigen::EigenBase<Derived>& x) const {
    double* x_mem_
      = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(x.size());

    return x_mem_;
  }

  inline Eigen::Map<Eigen::MatrixXd> map_matrix(double* mem, size_t rows, size_t cols) {
    return Eigen::Map<Eigen::MatrixXd>(mem, rows, cols);
  }

  inline Eigen::Map<Eigen::VectorXd> map_vector(double* mem, size_t rows) {
    return Eigen::Map<Eigen::VectorXd>(mem, rows);
  }

  inline Eigen::Map<Eigen::RowVectorXd> map_row_vector(double* mem, size_t cols) {
    return Eigen::Map<Eigen::RowVectorXd>(mem, cols);
  }
};

template <typename T, typename = void>
class adj_op;

template <typename T>
struct adj_op<T, require_var_vt<std::is_floating_point, T>> {
  using value_type = std::decay_t<plain_type_t<T>>;
  using Scalar = value_type_t<T>;
  using eigen_map = Eigen::Map<value_type>;
  static constexpr Eigen::Index RowsAtCompileTime{0};
  static constexpr Eigen::Index ColsAtCompileTime{0};
  std::decay_t<T>::value_type map_;
  template <typename S, require_floating_point_t<S>* = nullptr>
  explicit adj_op(S&& x) : map_(x) {}
  inline auto rows() const {
    return map_.rows();
  }

  inline auto cols() const {
    return map_.cols();
  }

  inline const auto& map() const {
    return map_;
  }

  inline auto& map() {
    return map_;
  }
};
template <typename T>
struct adj_op<T, require_eigen_vt<is_var, T>> {
  using value_type = std::decay_t<plain_type_t<T>>;
  using Scalar = typename value_type_t<T>::value_type;
  using eigen_map = Eigen::Map<value_type>;
  static constexpr Eigen::Index RowsAtCompileTime{value_type::RowsAtCompileTime};
  static constexpr Eigen::Index ColsAtCompileTime{value_type::ColsAtCompileTime};
  value_type* data_;
  eigen_map map_;

  template <typename S, require_eigen_t<S>* = nullptr>
  explicit adj_op(S&& x) : rows_(x.rows()), cols_(x.cols()),
    data_(stan::math::ChainableStack::instance_->memalloc_.alloc_array<Scalar>(rows_ * cols_)),
    map_(data_, rows_, cols_) {}

  inline auto rows() const {
    return map_.rows();
  }

  inline auto cols() const {
    return map_.cols();
  }

  inline const auto& map() const {
    return map_;
  }

  inline auto& map() {
    return map_;
  }

};

template <typename T>
struct adj_op<T, require_eigen_vt<is_arithmetic, T>> {
  using value_type = std::decay_t<plain_type_t<T>>;
  using Scalar = value_type_t<T>;
  using eigen_map = Eigen::Map<value_type>;
  static constexpr Eigen::Index RowsAtCompileTime{0};
  static constexpr Eigen::Index ColsAtCompileTime{0};

  template <typename S, require_eigen_t<S>* = nullptr>
  explicit adj_op(S&& x) {};
  inline Eigen::Index rows() const {
    return 0;
  }

  inline Eigen::Index cols() const {
    return 0;
  }

  inline auto map() const {
    return Eigen::Map<std::decay_t<plain_type_t<T>>>(nullptr, 0, 0);
  }
}
namespace internal {

/**
 * Compute the dimensionality of the given template argument. The
 * definition of dimensionality is deferred to specializations. By
 * default don't have a value (fail to compile)
 */
template <typename T>
struct compute_dims {};

/**
 * Compute the dimensionality of the given template argument. Double
 * types have dimensionality zero.
 */
template <>
struct compute_dims<double> {
  static constexpr size_t value = 0;
};

/**
 * Compute the dimensionality of the given template argument.
 * std::vector has dimension 1
 */
template <typename T>
struct compute_dims<std::vector<T>> {
  static constexpr size_t value = 1;
};

/**
 * Compute the dimensionality of the given template argument.
 * Eigen::Matrix types all have dimension two.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 */
template <typename T, int R, int C>
struct compute_dims<Eigen::Matrix<T, R, C>> {
  static constexpr size_t value = 2;
};
}  // namespace internal

/**
 * Compile time accumulator for std::arrays
 * @tparam T type held in the std::array.
 * @tparam N size in the std::array.
 * @param A std::array to accumulate over.
 * @param i Starting position to begin accumulation for `A`.
 */
template <typename T, size_t N>
constexpr size_t compile_time_accumulator(const std::array<T, N>& A,
                                          const int i = 0) {
  return (i < N) ? A[i] + compile_time_accumulator(A, i + 1) : size_t(0);
}

template <typename... Targs>
struct x_vis_alloc : vari {
};

/**
 * adj_jac_vari interfaces a user supplied functor  with the reverse mode
 * autodiff. It allows someone to implement functions with custom reverse mode
 * autodiff without having to deal with autodiff types.
 *
 * The requirements on the functor F are described in the documentation for
 * adj_jac_apply
 *
 * @tparam F class of functor
 * @tparam Targs types of arguments: can be any mix of double, var, or
 * Eigen::Matrices with double or var scalar components
 */
template <typename F, typename... Targs>
struct adj_jac_vari : public vari {
  // holds whether each input type held a vari_value
  static constexpr std::array<bool, sizeof...(Targs)> is_var_{
    {is_var_value<scalar_type_t<Targs>>::value...}};
  using FReturnType =
    std::decay_t<plain_type_t<std::result_of_t<F(decltype(is_var_), decltype(value_of(std::declval<Targs>()))...)>>>;

  F f_;  // Function to be invoked

  using x_vis_tuple_ = var_to_vari_filter_t<std::decay_t<Targs>...>;
  x_vis_tuple_ x_vis_{};  // tuple holding pointers to input are var mem.

  int M_;

  using y_vari_type = std::conditional_t<is_std_vector<FReturnType>::value,
					 vari**,
					 vari_value<FReturnType>*>;
  y_vari_type y_vi_;  // vari pointer for output.

inline void fill_adj_jac() {}

  template <typename Arith, typename... Pargs,
            require_st_arithmetic<Arith>* = nullptr>
  inline void fill_adj_jac(const Arith& x, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    std::get<t>(x_vis_) = nullptr;
    fill_adj_jac(args...);
  }

  template <typename T, typename... Pargs>
  inline void fill_adj_jac(const var_value<T>& x, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    std::get<t>(x_vis_) = x.vi_;
    fill_adj_jac(args...);
  }

  // std::vector<var_value>
  template <typename... Pargs>
  inline void fill_adj_jac(const std::vector<var>& x, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;

    std::get<t>(x_vis_)
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(x.size());

    M_ = x.size();

    save_varis(std::get<t>(x_vis_), x);

    fill_adj_jac(args...);
  }

  template <int R, int C, typename... Pargs>
  inline void fill_adj_jac(const Eigen::Matrix<var, R, C>& x,
			   const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;

    std::get<t>(x_vis_)
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(x.size());

    save_varis(std::get<t>(x_vis_), x);

    fill_adj_jac(args...);
  }

  explicit adj_jac_vari(const Targs&... args)
    : vari(NOT_A_NUMBER), y_vi_(nullptr), M_(0) {
    fill_adj_jac(args...);
  }

  inline var build_return_varis_and_vars(double val_y) {
    const_cast<double&>(val_) = val_y;
    y_vi_ = static_cast<vari*>(this);//new vari(val_y, false);
    return this;
  }

  template <typename Derived,
	    typename T = Eigen::Matrix<double,
				       Eigen::MatrixBase<Derived>::RowsAtCompileTime,
				       Eigen::MatrixBase<Derived>::ColsAtCompileTime>>
  inline var_value<T>
  build_return_varis_and_vars(const Eigen::MatrixBase<Derived>& val_y) {
    y_vi_ = new vari_value<T>(val_y, false);
    return y_vi_;
  }

  /**
   * Return a std::vector of vars created from newly allocated varis initialized
   * with the values of val_y
   *
   * @param val_y output of F::operator()
   * @return std::vector of vars
   */
  template <typename T>
  inline std::vector<var_value<T>>
  build_return_varis_and_vars(const std::vector<T>& val_y) {
    std::vector<var_value<T>> var_y;
    var_y.reserve(val_y.size());
    y_vi_ = ChainableStack::instance_->memalloc_.alloc_array<vari_base*>
      (val_y.size());
    for (size_t m = 0; m < val_y.size(); ++m) {
      y_vi_[m] = new vari_value<T>(val_y[m], false);
      var_y.emplace_back(static_cast<vari_value<T>*>(y_vi_[m]));
    }
    return var_y;
  }

  /**
   * The adj_jac_vari functor
   *  1. Initializes an instance of the user defined functor F
   *  2. Calls operator() on the F instance with the double values from the
   * input args
   *  3. Saves copies of the varis pointed to by the input vars for subsequent
   * calls to chain
   *  4. Calls build_return_varis_and_vars to construct the appropriate output
   * data structure of vars
   *
   * Each of the arguments can be an Eigen::Matrix with var or double scalar
   * types, a std::vector with var, double, or int scalar types, or a var, a
   * double, or an int.
   *
   * @param args Input arguments
   * @return Output of f_ as vars
   */
  inline auto operator()(const Targs&... args) {
    return build_return_varis_and_vars(f_(is_var_, value_of(args)...));
  }

  template <typename T>
  inline void accumulate_adjoints_in_varis_impl(std::nullptr_t, T&&) {}

  template <int R1, int C1, typename Derived>
  inline void accumulate_adjoints_in_varis_impl
  (vari_value<Eigen::Matrix<double, R1, C1>> *vi, const Eigen::MatrixBase<Derived>& y_adj_jac) {
    vi->adj_ += y_adj_jac;
  }

  template <typename Derived>
  inline void accumulate_adjoints_in_varis_impl
  (vari **vis, const Eigen::MatrixBase<Derived>& y_adj_jac) {
    Eigen::Map<
      Eigen::Matrix<vari*,
		    Eigen::MatrixBase<Derived>::RowsAtCompileTime,
		    Eigen::MatrixBase<Derived>::ColsAtCompileTime>>
      (vis, y_adj_jac.rows(), y_adj_jac.cols()).adj() += y_adj_jac;
  }

  inline void accumulate_adjoints_in_varis_impl(vari **vis,
						const std::vector<double>& y_adj_jac) {
    for (int n = 0; n < y_adj_jac.size(); ++n) {
      vis[n]->adj_ += y_adj_jac[n];
    }
  }

  inline void accumulate_adjoints_in_varis_impl(vari *vi, double y_adj_jac) {
    vi->adj_ += y_adj_jac;
  }

  inline void accumulate_adjoints_in_varis() {}

  template <typename T,
	    typename... Pargs>
  inline void accumulate_adjoints_in_varis(const T& y_adj_jac,
					   const Pargs&... args) {
    static constexpr size_t t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    if(is_var_[t])
      accumulate_adjoints_in_varis_impl(std::get<t>(x_vis_), y_adj_jac);
    accumulate_adjoints_in_varis(args...);
  }

  /**
   * Store the adjoint in y_vi[0] in y_adj
   *
   * @tparam size dimensionality of M
   * @param[in] y_vi pointer to pointer to vari
   * @param[in] M shape of y_adj
   * @param[out] y_adj reference to variable where adjoint is to be stored
   */
  template <typename T>
  inline const auto& build_y_adj(vari_value<T>* y_vi, int M) {
    return y_vi->adj_;
  }

  /**
   * Store the adjoints from y_vi in y_adj
   *
   * @tparam size dimensionality of M
   * @param[in] y_vi pointer to pointers to varis
   * @param[in] M shape of y_adj
   * @param[out] y_adj reference to std::vector where adjoints are to be stored
   */
  inline auto build_y_adj(vari** y_vi, int M) {
    std::vector<double> y_adj;
    y_adj.reserve(M);
    for (size_t m = 0; m < M; ++m) {
      y_adj.emplace_back(static_cast<vari*>(y_vi[m])->adj_);
    }
    return y_adj;
  }

  /**
   * Propagate the adjoints at the output varis (y_vi_) back to the input
   * varis (x_vis_) by:
   * 1. packing the adjoints in an appropriate container using build_y_adj
   * 2. using the multiply_adjoint_jacobian function of the user defined
   * functor to compute what the adjoints on x_vis_ should be
   * 3. accumulating the adjoints into the varis pointed to by elements of
   * x_vis_ using accumulate_adjoints_in_varis
   *
   * This operation may be called multiple times during the life of the vari
   */
  inline void chain() {
    apply(
        [&, this](auto&&... args) {
          this->accumulate_adjoints_in_varis(args...);
        },
        f_.multiply_adjoint_jacobian(is_var_, build_y_adj(y_vi_, M_)));
  }
};

template <typename F, typename... Targs>
constexpr std::array<bool, sizeof...(Targs)> adj_jac_vari<F, Targs...>::is_var_;

/**
 * Return the result of applying the function defined by a nullary
 * construction of F to the specified input argument
 *
 * adj_jac_apply makes it possible to write efficient reverse-mode
 * autodiff code without ever touching Stan's autodiff internals
 *
 * Mathematically, to use a function in reverse mode autodiff, you need to be
 * able to evaluate the function (y = f(x)) and multiply the Jacobian of that
 * function (df(x)/dx) by a vector.
 *
 * As an example, pretend there exists some large, complicated function, L(x1,
 * x2), which contains our smaller function f(x1, x2). The goal of autodiff is
 * to compute the partials dL/dx1 and dL/dx2. If we break the large function
 * into pieces:
 *
 * y = f(x1, x2)
 * L = g(y)
 *
 * If we were given dL/dy we could compute dL/dx1 by the product dL/dy *
 * dy/dx1 or dL/dx2 by the product dL/dy * dy/dx2
 *
 * Because y = f(x1, x2), dy/dx1 is just df(x1, x2)/dx1, the Jacobian of the
 * function we're trying to define with x2 held fixed. A similar thing happens
 * for dy/dx2. In vector form,
 *
 * dL/dx1 = (dL/dy)' * df(x1, x2)/dx1 and
 * dL/dx2 = (dL/dy)' * df(x1, x2)/dx2
 *
 * So implementing f(x1, x2) and the products above are all that is required
 * mathematically to implement reverse-mode autodiff for a function.
 *
 * adj_jac_apply takes as a template argument a functor F that supplies the
 * non-static member functions (leaving exact template arguments off):
 *
 * (required) Eigen::VectorXd operator()(const std::array<bool, size>&
 * needs_adj, const T1& x1..., const T2& x2, ...)
 *
 * where there can be any number of input arguments x1, x2, ... and T1, T2,
 * ... can be either doubles or any Eigen::Matrix type with double scalar
 * values. needs_adj is an array of size equal to the number of input
 * arguments indicating whether or not the adjoint^T Jacobian product must be
 * computed for each input argument. This argument is passed to operator() so
 * that any unnecessary preparatory calculations for multiply_adjoint_jacobian
 * can be avoided if possible.
 *
 * (required) std::tuple<T1, T2, ...> multiply_adjoint_jacobian(const
 * std::array<bool, size>& needs_adj, const Eigen::VectorXd& adj)
 *
 * where T1, T2, etc. are the same types as in operator(), needs_adj is the
 * same as in operator(), and adj is the vector dL/dy.
 *
 * operator() is responsible for computing f(x) and multiply_adjoint_jacobian
 * is responsible for computing the necessary adjoint transpose Jacobian
 * products (which frequently does not require the calculation of the full
 * Jacobian).
 *
 * operator() will be called before multiply_adjoint_jacobian is called, and
 * is only called once in the lifetime of the functor
 * multiply_adjoint_jacobian is called after operator() and may be called
 * multiple times for any single functor
 *
 * The functor supplied to adj_jac_apply must be careful to allocate any
 * variables it defines in the autodiff arena because its destructor will
 * never be called and memory will leak if allocated anywhere else.
 *
 * Targs (the input argument types) can be any mix of doubles, vars, ints,
 * std::vectors with double, var, or int scalar components, or
 * Eigen::Matrix s of any shape with var or double scalar components
 *
 * @tparam F functor to be connected to the autodiff stack
 * @tparam Targs types of arguments to pass to functor
 * @param args input to the functor
 * @return the result of the specified operation wrapped up in vars
 */
template <typename T,
	  require_not_eigen_vt<is_var, T>* = nullptr>
const T& convert_to_whole_matrix(const T& arg) {
  return arg;
}

template <typename T,
	  require_eigen_vt<is_var, T>* = nullptr>
var_value<Eigen::Matrix<double,
			T::RowsAtCompileTime,
			T::ColsAtCompileTime>>
convert_to_whole_matrix(const T& arg) {
  return arg;
}

template <typename F, typename... Targs>
inline auto adj_jac_apply_impl(const Targs&... args) {
  auto vi = new adj_jac_vari<F, Targs...>(args...);

  return (*vi)(args...);
}

template <typename F, typename... Targs>
inline auto adj_jac_apply(const Targs&... args) {
  return adj_jac_apply_impl<F>(args...);
}

}  // namespace math
}  // namespace stan
#endif
