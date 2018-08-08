#ifndef STAN_MATH_REV_MAT_FUN_ADJ_JAC_APPLY_HPP
#define STAN_MATH_REV_MAT_FUN_ADJ_JAC_APPLY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat.hpp>
#include <stan/math/rev/core.hpp>
#include <functional>
#include <limits>
#include <vector>

namespace stan {
namespace math {

namespace {
/*
 * Placeholder apply function (until c++17 is available for Stan)
 */
template <class F, class Tuple, std::size_t... I>
constexpr auto apply_impl(const F& f, const Tuple& t,
                          std::index_sequence<I...>) {
  return f(std::get<I>(t)...);
}

template <class F, class Tuple>
constexpr auto apply(const F& f, const Tuple& t) {
  return apply_impl(f, t, std::make_index_sequence<std::tuple_size<Tuple>{}>{});
}
}  // namespace

/*
 * adj_jac_vari interfaces a user supplied functor  with the reverse mode
 * autodiff. It allows someone to implement functions with custom reverse mode
 * autodiff without having to deal with autodiff types.
 *
 * The requirements on the functor F are described in the documentation for
 * adj_jac_apply
 *
 * Targs (the input argument types) can be any mix of double, var, or
 * Eigen::Matrices with double or var scalar components
 *
 * @tparam F class of functor
 * @tparam Targs Types of arguments
 */
template <typename F, typename... Targs>
struct adj_jac_vari : public vari {
  F f_;
  std::array<bool, sizeof...(Targs)> is_var_;
  std::array<int, sizeof...(Targs)> offsets_;
  vari** x_vis_;
  int M_;
  vari** y_vi_;

  /*
   * count_memory is a recursive templated function that accumulates in its
   * first argument (count) the number of varis used in the rest of its
   * arguments
   *
   * This is used to figure out how much space to allocate in x_vis_.
   *
   * The array offsets_ is populated with values to indicate where in x_vis_ the
   * vari pointers for each argument will be stored.
   *
   * As an example, take the call:
   *    size_t s = count_memory(0, a1, a2, a3)
   *
   * If a1 is an Eigen::Matrix<var, 5, 5>, a2 is a double, and a3 is a var, the
   * return value of count_memory is 26 and offsets_ is set to {0, 25, 25}
   *
   * Each of the arguments can be an Eigen::Matrix<var, R, C>, an
   * Eigen::Matrix<double, R, C>, a var, or a double.
   *
   * @param count rolling count of number of varis that must be allocated
   * @param x next argument to have its varis counted
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <int R, int C, typename... Pargs>
  size_t count_memory(size_t count, const Eigen::Matrix<var, R, C>& x,
                      const Pargs&... args) {
    constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    count += x.size();
    return count_memory(count, args...);
  }

  template <int R, int C, typename... Pargs>
  size_t count_memory(size_t count, const Eigen::Matrix<double, R, C>& x,
                      const Pargs&... args) {
    constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    return count_memory(count, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const var& x, const Pargs&... args) {
    constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    count += 1;
    return count_memory(count, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const double& x, const Pargs&... args) {
    constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    return count_memory(count, args...);
  }

  size_t count_memory(size_t count) { return count; }

  /*
   * prepare_x_vis recursively populates x_vis_ with the varis from each of its
   * input arguments. The vari pointers for argument n are copied into x_vis_ at
   * the index starting at offsets_[n].
   *
   * Each of the arguments can be an Eigen::Matrix<var, R, C>, an
   * Eigen::Matrix<double, R, C>, a var, or a double.
   *
   * As an example, take the call:
   *    prepare_x_vis(a1, a2, a3)
   *
   * If a1 is an Eigen::Matrix<var, 5, 5>, a2 is a double, and a3 is a var,
   * x_vis_ to x_vis_ + 24 is filled with the vari pointers of a1, and
   * x_vis_[25] is set equal to the vari pointer of a3. is_var_ is set to {true,
   * false, true}
   *
   * @param x next argument to have its vari pointers copied if necessary
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <int R, int C, typename... Pargs>
  void prepare_x_vis(const Eigen::Matrix<var, R, C>& x, const Pargs&... args) {
    constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    for (int i = 0; i < x.size(); ++i) {
      x_vis_[offsets_[t] + i] = x(i).vi_;
    }
    prepare_x_vis(args...);
  }

  template <int R, int C, typename... Pargs>
  void prepare_x_vis(const Eigen::Matrix<double, R, C>& x,
                     const Pargs&... args) {
    prepare_x_vis(args...);
  }

  template <typename... Pargs>
  void prepare_x_vis(const var& x, const Pargs&... args) {
    constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    x_vis_[offsets_[t]] = x.vi_;
    prepare_x_vis(args...);
  }

  template <typename... Pargs>
  void prepare_x_vis(const double& x, const Pargs&... args) {
    prepare_x_vis(args...);
  }

  void prepare_x_vis() {}

  /**
   * The adj_jac_vari constructor
   *  1. Initializes an instance of the user defined functor F
   *  2. Calls operator() on the F instance with the double values from the
   * input args
   *  3. Saves copies of the varis pointed to by the input vars for subsequent
   * calls to chain
   *  4. Allocates varis for the output of the functor F
   *  5. Populates the array is_var_ with whether or not the scalar type in each
   * argument is a var
   *
   * The input argument types can be any mix of double, var, or Eigen::Matrices
   * with double or var scalar components
   *
   * As an example of how is_var_ is initialized, take the call:
   *    prepare_x_vis(a1, a2, a3)
   *
   * If a1 is an Eigen::Matrix<var, 5, 5>, a2 is a double, and a3 is a var,
   * is_var_ is set to {true, false, true}
   *
   * @param args Input arguments
   */
  explicit adj_jac_vari(const Targs&... args)
      : vari(std::numeric_limits<double>::quiet_NaN()),  // The val_ in this
                                                         // vari is unused
        is_var_({is_var<typename scalar_type<Targs>::type>::value...}),
        x_vis_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            count_memory(0, args...))) {
    prepare_x_vis(args...);

    Eigen::Matrix<double, Eigen::Dynamic, 1> val_y
        = f_(is_var_, value_of(args)...);

    M_ = val_y.size();
    y_vi_ = ChainableStack::instance().memalloc_.alloc_array<vari*>(M_);
    for (int m = 0; m < M_; ++m) {
      y_vi_[m] = new vari(val_y(m), false);
    }
  }

  /*
   * accumulate_adjoints recurses through the adjoints returned from the
   * user-supplied adjoint_jacobian_multiply and accumulates the appropriate
   * vari adjoints pointed to by values of x_vis_.
   *
   * The adjoints of the nth argument are accumulated into the adjoints of the
   * varis pointed to in the range x_vis_ + offsets_[n] to x_vis_ + offsets[n] +
   * size of nth argument.
   *
   * If is_var_[n] is false, then the nth argument is ignored, whatever its
   * value is.
   *
   * Each of the arguments can be an Eigen::Matrix<double, R, C> or a double.
   *
   * As an example, take the call:
   *    accumulate_adjoints(a1, a2, a3)
   *
   * If a1 is an Eigen::Matrix<double, 5, 5>, a2 is a double, a3 is a double,
   * and is_var_ is set to {true, false, true}, the 25 values in a1 are
   * accumulated in the adjs of the varis pointed to by the values of x_vis_ to
   * x_vis_ + 24, the value of a2 is ignored, and the adjoint of the vari
   * pointed to by x_vis_[25] is incremented by a3.
   *
   * @param y_adj_jac next set of adjoints to be accumulated
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <int R, int C, typename... Pargs>
  void accumulate_adjoints(const Eigen::Matrix<double, R, C>& y_adj_jac,
                           const Pargs&... args) {
    const int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    if (is_var_[t]) {
      for (int n = 0; n < y_adj_jac.size(); ++n) {
        x_vis_[offsets_[t] + n]->adj_ += y_adj_jac(n);
      }
    }

    accumulate_adjoints(args...);
  }

  template <typename... Pargs>
  void accumulate_adjoints(const double& y_adj_jac, const Pargs&... args) {
    const int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    if (is_var_[t]) {
      x_vis_[offsets_[t]]->adj_ += y_adj_jac;
    }

    accumulate_adjoints(args...);
  }

  void accumulate_adjoints() {}

  /**
   * chain propagates the adjoints at the output varis (y_vi_) back to the input
   * varis (x_vis_) by using the multiply_adjoint_jacobian function of the user
   * defined functor
   *
   * Unlike the constructor, this operation may be called multiple times during
   * the life of the vari
   */
  void chain() {
    Eigen::Matrix<double, Eigen::Dynamic, 1> y_adj(M_);
    for (int m = 0; m < M_; ++m)
      y_adj(m) = y_vi_[m]->adj_;
    auto y_adj_jacs = f_.multiply_adjoint_jacobian(is_var_, y_adj);

    /**
     * NOTE TO REVIEWER: If the cost of this lambda function redirection seems
     * to high, we could make accumulate_adjoints a non member function and just
     * pass it all the info it needs.
     */
    apply([this](auto&&... args) { this->accumulate_adjoints(args...); },
          y_adj_jacs);
  }
};

/*
 * Return the result of applying the function defined by a nullary construction
 * of F to the specified input argument
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
 * If we were given dL/dy we could compute dL/dx1 by the product dL/dy * dy/dx1
 * or dL/dx2 by the product dL/dy * dy/dx2
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
 * where there can be any number of input arguments x1, x2, ... and T1, T2, ...
 * can be either doubles or any Eigen::Matrix type with double scalar values.
 * needs_adj is an array of size equal to the number of input arguments
 * indicating whether or not the adjoint^T Jacobian product must be computed for
 * each input argument. This argument is passed to operator() so that any
 * unnecessary preparatory calculations for multiply_adjoint_jacobian can be
 * avoided if possible.
 *
 * (required) std::tuple<T1, T2, ...> multiply_adjoint_jacobian(const
 * std::array<bool, size>& needs_adj, const Eigen::VectorXd& adj)
 *
 * where T1, T2, etc. are the same types as in operator(), needs_adj is the same
 * as in operator(), and adj is the vector dL/dy.
 *
 * operator() is responsible for computing f(x) and multiply_adjoint_jacobian is
 * responsible for computing the necessary adjoint transpose Jacobian products
 * (which frequently does not require the calculation of the full Jacobian).
 *
 * operator() will be called before multiply_adjoint_jacobian is called, and is
 * only called once in the lifetime of the functor multiply_adjoint_jacobian is
 * called after operator() and may be called multiple times for any single
 * functor
 *
 * The functor supplied to adj_jac_apply must be careful to allocate any
 * variables it defines in the autodiff arena because its destructor will
 * never be called and memory will leak if allocated anywhere else.
 *
 * Targs (the input argument types) can be any mix of double, var, or
 * Eigen::Matrices with double or var scalar components
 *
 * @tparam F functor to be connected to the autodiff stack
 * @tparam Targs types of arguments to pass to functor
 * @param args input to the functor
 * @return the result of the specified operation wrapped up in vars
 */
template <typename F, typename... Targs>
Eigen::Matrix<var, Eigen::Dynamic, 1> adj_jac_apply(const Targs&... args) {
  auto vi = new adj_jac_vari<F, Targs...>(args...);
  Eigen::Matrix<var, Eigen::Dynamic, 1> y(vi->M_);

  for (int m = 0; m < y.size(); ++m)
    y(m) = (vi->y_vi_)[m];
  return y;
}

}  // namespace math
}  // namespace stan
#endif
