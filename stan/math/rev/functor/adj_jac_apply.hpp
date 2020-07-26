#ifndef STAN_MATH_REV_FUNCTOR_ADJ_JAC_APPLY_HPP
#define STAN_MATH_REV_FUNCTOR_ADJ_JAC_APPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/adj_op.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {
/**
 * Based on the container type return either a vari* or vari**
 */
template <typename T, typename = void>
struct var_return_type;
template <typename T>
struct var_return_type<T, require_arithmetic_t<T>> {
  using type = vari*;
};
template <typename T>
struct var_return_type<T, require_std_vector_t<T>> {
  using type = vari**;
};
template <typename T>
struct var_return_type<T, require_eigen_t<T>> {
  using type = vari**;
};

}  // namespace internal

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
  static constexpr std::array<bool, sizeof...(Targs)> is_var_{
      is_var<scalar_type_t<Targs>>::value...};
  using FReturnType = std::result_of_t<F(plain_type_t<decltype(value_of(plain_type_t<Targs>()))>...)>;
  using x_vis_tuple_ = var_to_vari_filter_t<std::decay_t<Targs>...>;
  F f_; // Struct that with methods for computing forward and reverse pass.
  x_vis_tuple_ x_vis_; // tuple holding pointers to vari
  /**
   * Will be either vari* or vari**
   */
  typename internal::var_return_type<FReturnType>::type y_vi_;
  /**
   * Stores the dimensions of return value of forward pass to reconstruct
   * container for reverse pass.
   */
  std::array<size_t, 2> dims_;
  /**
   * Return the adjoint
   * @tparam RetType Not called by user. Default of the return type of
   * operator() of the functor F allows for SFINAE to choose the correct
   * specialization given the required output type.
   */
  template <typename RetType = FReturnType,
            require_arithmetic_t<RetType>* = nullptr>
  inline auto y_adj() {
    return y_vi_->adj_;
  }

  /**
   * Return a standard vector of the adjoints
   * @tparam RetType Not called by user. Default of the return type of
   * operator() of the functor F allows for SFINAE to choose the correct
   * specialization given the required output type.
   */
  template <typename RetType = FReturnType,
            require_std_vector_t<RetType>* = nullptr>
  inline auto y_adj() {
    std::vector<double> var_y(dims_[0]);
    for (size_t i = 0; i < dims_[0]; ++i) {
      var_y[i] = y_vi_[i]->adj_;
    }
    return var_y;
  }

  /**
   * Return an eigen matrix of the adjoints
   * @tparam RetType Not called by user. Default of the return type of
   * operator() of the functor F allows for SFINAE to choose the correct
   * specialization given the required output type.
   */
  template <typename RetType = FReturnType, require_eigen_t<RetType>* = nullptr>
  inline auto y_adj() {
    Eigen::Matrix<double, FReturnType::RowsAtCompileTime,
                  FReturnType::ColsAtCompileTime>
        var_y(dims_[0], dims_[1]);
    const size_t iter_size{dims_[0] * dims_[1]};
    for (size_t i = 0; i < iter_size; ++i) {
      var_y(i) = y_vi_[i]->adj_;
    }
    return var_y;
  }

  /**
   * prepare_x_vis populates x_vis_ with the varis from each of its
   * input arguments. The vari pointers for argument n are copied into x_vis_ at
   * the index starting at offsets_[n]. For Eigen::Matrix types, this copying is
   * done in with column major ordering.
   *
   * Each of the arguments can be an Eigen::Matrix with var or double scalar
   * types, a std::vector with var, double, or int scalar types, or a var, a
   * double, or an int.
   *
   * @tparam R number of rows, can be Eigen::Dynamic
   * @tparam C number of columns, can be Eigen::Dynamic
   * @tparam Pargs Types of the rest of the arguments to be processed
   *
   * @param x next argument to have its vari pointers copied if necessary
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
  inline auto prepare_x_vis_impl(const EigMat& x) {
    vari** y
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(x.size());
    for (int i = 0; i < x.size(); ++i) {
      y[i] = x(i).vi_;
    }
    return y;
  }

  inline auto prepare_x_vis_impl(const std::vector<var>& x) {
    vari** y
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      y[i] = x[i].vi_;
    }
    return y;
  }

  inline auto prepare_x_vis_impl(const var& x) {
    vari* y = ChainableStack::instance_->memalloc_.alloc_array<vari>(1);
    y = x.vi_;
    return y;
  }

  template <typename Arith, require_st_arithmetic<Arith>* = nullptr>
  inline auto prepare_x_vis_impl(Arith&& x) {
    return nullptr;
  }

  auto prepare_x_vis() { return x_vis_tuple_{}; }

  template <typename... Args>
  inline x_vis_tuple_ prepare_x_vis(Args&&... args) {
    return x_vis_tuple_{prepare_x_vis_impl(args)...};
  }
  /**
   * Specialization for var, stores the values of the forward pass in
   * `adj_jac_vari` `y_vi_`.
   * @param val_y Return value of `F()` stored as the value in `y_vi_`.
   */
  inline auto& collect_forward_pass(const double& val_y) {
    y_vi_ = new vari(val_y, false);
    return y_vi_;
  }
  /**
   * Specialization for std::vector<var>, stores the values of the forward pass in
   * `adj_jac_vari` `y_vi_`.
   * @param val_y Return value of `F()` stored as the value in `y_vi_`.
   */
  template <typename Vec, require_std_vector_vt<std::is_arithmetic, Vec>* = nullptr>
  inline auto& collect_forward_pass(Vec&& val_y) {
    y_vi_
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(val_y.size());
    this->dims_[0] = val_y.size();
    for (size_t i = 0; i < val_y.size(); ++i) {
      y_vi_[i] = new vari(val_y[i], false);
    }
    return y_vi_;
  }

  /**
   * Specialization for Eigen matrices of vars, stores the values of the forward pass in
   * `adj_jac_vari` `y_vi_`.
   * @param val_y Return value of `F()` stored as the value in `y_vi_`.
   */
  template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
  inline auto& collect_forward_pass(EigMat&& val_y) {
    using eig_mat = std::decay_t<EigMat>;
    using ret_mat = Eigen::Matrix<var, eig_mat::RowsAtCompileTime,
                                  eig_mat::ColsAtCompileTime>;
    y_vi_
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(val_y.size());
    this->dims_[0] = val_y.rows();
    this->dims_[1] = val_y.cols();
    for (size_t i = 0; i < val_y.size(); ++i) {
      y_vi_[i] = new vari(val_y(i), false);
    }
    return y_vi_;
  }
  /**
   * Given a parameter pack used in a recursive function, deduce the position
   *  of the corresponding value in the `is_var_` array.
   */
  template <typename... Pargs>
  struct var_idx_ {
    static constexpr size_t value{sizeof...(Targs) - sizeof...(Pargs) - 1};
  };

  /**
   * Given a parameter pack used in a recursive function, deduce the value
   *  for the corresponding type in the `is_var_` array.
   */
  template <typename... Pargs>
  struct is_var_check_ {
    /**
     * If true, the corresponding type is a var or container of vars.
     */
    static constexpr bool value{is_var_[var_idx_<Pargs...>::value]};
  };

  /**
   * Implimentation of for_each for applying adjoints.
   * @note The static cast to void is used in boost::hana's for_each impl
   *  and is used to suppress unused value warnings from the compiler.
   */
  template <typename FF, typename T1, typename T2, size_t... Is>
  constexpr inline auto for_each_adj_impl(FF&& f, T1&& x, T2&& t, std::index_sequence<Is...>) {
    using Swallow = int[];
    static_cast<void>(Swallow{(static_cast<void>(f(std::get<Is>(std::forward<T1>(x)), std::get<Is>(std::forward<T2>(t)), bool_constant<is_var_[Is]>())),0)...});
  }
  /**
   * Apply a function and object to each element of a tuple
   * @tparam F type with a valid `operator()`
   * @tparam T1 The first type to be used in the argument list of `F`
   * @tparam T2 tuple
   * @param f A functor to apply over each element of the tuple.
   * @param x Any object
   * @param t A tuple.
   */
  template <typename FF, typename T1, typename T2>
  constexpr inline auto for_each_adj(FF&& f, T1&& x, T2&& t) {
    return for_each_adj_impl(
        std::forward<FF>(f), std::forward<T1>(x), std::forward<T2>(t),
        std::make_index_sequence<std::tuple_size<std::decay_t<T2>>::value>());
  }

  /**
   * Accumulate, if necessary, the values of y_adj_jac into the
   * adjoints of the varis pointed to by the appropriate elements
   * of x_vis_. Recursively calls accumulate_adjoints on the rest of the
   * arguments.
   *
   * @tparam R number of rows, can be Eigen::Dynamic
   * @tparam C number of columns, can be Eigen::Dynamic
   * @tparam Pargs Types of the rest of adjoints to accumulate
   *
   * @param y_adj_jac set of values to be accumulated in adjoints
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  struct accumulate_adjoint {
    template <typename Toss, typename T, bool needs_adj, std::enable_if_t<!needs_adj>* = nullptr>
    inline void operator()(Toss&&, T&& x_vis, bool_constant<needs_adj>) {
    }

    template <typename EigMat, typename T, bool needs_adj, std::enable_if_t<needs_adj>* = nullptr,
              require_eigen_t<EigMat>* = nullptr>
    inline void operator()(EigMat&& y_adj_jac, T&& x_vis, bool_constant<needs_adj>) {
      for (int n = 0; n < y_adj_jac.size(); ++n) {
        x_vis[n]->adj_ += y_adj_jac(n);
      }
    }

    /**
     * Accumulate, if necessary, the values of y_adj_jac into the
     * adjoints of the varis pointed to by the appropriate elements
     * of x_vis_. Recursively calls accumulate_adjoints on the rest of the
     * arguments.
     *
     * @tparam Pargs Types of the rest of adjoints to accumulate
     * @param y_adj_jac set of values to be accumulated in adjoints
     * @param args the rest of the arguments (that will be iterated through
     * recursively)
     */
    template <typename Vec, typename T, bool needs_adj, std::enable_if_t<needs_adj>* = nullptr,
              require_std_vector_t<Vec>* = nullptr>
    inline void operator()(Vec&& y_adj_jac, T&& x_vis, bool_constant<needs_adj>) {
      for (int n = 0; n < y_adj_jac.size(); ++n) {
        x_vis[n]->adj_ += y_adj_jac[n];
      }
    }

    /**
     * Accumulate, if necessary, the value of y_adj_jac into the
     * adjoint of the vari pointed to by the appropriate element
     * of x_vis_. Recursively calls accumulate_adjoints on the rest of the
     * arguments.
     *
     * @tparam Pargs Types of the rest of adjoints to accumulate
     * @param y_adj_jac next set of adjoints to be accumulated
     * @param args the rest of the arguments (that will be iterated through
     * recursively)
     */
    template <typename T, bool needs_adj, std::enable_if_t<needs_adj>* = nullptr>
    inline void operator()(const double& y_adj_jac, T&& x_vis, bool_constant<needs_adj>) {
      x_vis->adj_ += y_adj_jac;
    }

    inline void operator()() {}
  } accum_jac;
  template <typename RetType = FReturnType,
            require_arithmetic_t<RetType>* = nullptr>
  inline auto y_var() {
    return var(y_vi_);
  }

  template <typename RetType = FReturnType,
            require_std_vector_t<RetType>* = nullptr>
  inline auto y_var() {
    std::vector<var> var_y(dims_[0]);
    for (size_t i = 0; i < dims_[0]; ++i) {
      var_y[i] = y_vi_[i];
    }
    return var_y;
  }

  template <typename RetType = FReturnType, require_eigen_t<RetType>* = nullptr>
  inline auto y_var() {
    Eigen::Matrix<var, FReturnType::RowsAtCompileTime,
                  FReturnType::ColsAtCompileTime>
        var_y(dims_[0], dims_[1]);
    const size_t iter_size{dims_[0] * dims_[1]};
    for (size_t i = 0; i < iter_size; ++i) {
      var_y(i) = y_vi_[i];
    }
    return var_y;
  }

  /**
   * The adj_jac_vari constructor
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
  template <typename... Args>
  adj_jac_vari(Args&&... args)
      : vari(NOT_A_NUMBER),
        f_(args...),
        x_vis_(prepare_x_vis(args...)),
        y_vi_(collect_forward_pass(f_(value_of(args)...))) {}

  /**
   * Propagate the adjoints at the output varis (y_vi_) back to the input
   * varis (x_vis_) by:
   * 1. packing the adjoints in an appropriate container using build_y_adj
   * 2. using the multiply_adjoint_jacobian function of the user defined
   * functor to compute what the adjoints on x_vis_ should be
   * 3. accumulating the adjoints into the varis pointed to by elements of
   * x_vis_ using accumulate_adjoints
   *
   * This operation may be called multiple times during the life of the vari
   */
  inline void chain() {
    for_each_adj(accum_jac, f_.multiply_adjoint_jacobian(this->y_adj()), x_vis_);
/*    apply(
        [this](auto&&... args) {
          this->accumulate_adjoints(std::forward<decltype(args)>(args)...);
        },
      );*/
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
template <typename F, typename... Targs>
inline auto adj_jac_apply(const Targs&... args) {
  auto* vi = new adj_jac_vari<F, Targs...>(args...);

  return vi->y_var();
}

}  // namespace math
}  // namespace stan
#endif
