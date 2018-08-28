#ifndef STAN_MATH_REV_MAT_FUN_ADJ_JAC_APPLY_HPP
#define STAN_MATH_REV_MAT_FUN_ADJ_JAC_APPLY_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <limits>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace {
/*
 * Invoke the functor f with arguments given in t and indexed in the index
 * sequence I
 *
 * @tparam F Type of functor
 * @tparam Tuple Type of tuple containing arguments
 * @tparam I Index sequence going from 0 to std::tuple_size<T>::value - 1
 * inclusive
 * @param f functor callable
 * @param t tuple of arguments
 * @param i placeholder variable for index sequence
 */
template <class F, class Tuple, std::size_t... I>
constexpr auto apply_impl(const F& f, const Tuple& t,
                          std::index_sequence<I...> i) {
  return f(std::get<I>(t)...);
}

/*
 * Call the functor f with the tuple of arguments t, like:
 *
 * f(std::get<0>(t), std::get<1>(t), ...)
 *
 * TODO: replace this with implementation in C++ std when C++17 is available
 *
 * @tparam F Type of functor
 * @tparam Tuple Type of tuple containing arguments
 * @param f functor callable
 * @param t tuple of arguments
 */
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
   * count_memory returns count (the first argument) + the number of varis used
   * in the second argument + the number of arguments used to encode the
   * variadic tail args.
   *
   * The adj_jac_vari constructor uses this to figure out how much space to
   * allocate in x_vis_.
   *
   * The array offsets_ is populated with values to indicate where in x_vis_ the
   * vari pointers for each argument will be stored.
   *
   * Each of the arguments can be an Eigen::Matrix with var or double scalar
   * types, a std::vector with var, double, or int scalar types, or a var, a
   * double, or an int.
   *
   * @tparam R Eigen Matrix row type
   * @tparam C Eigen Matrix column type
   * @tparam Pargs Types of rest of arguments
   * @param count rolling count of number of varis that must be allocated
   * @param x next argument to have its varis counted
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <int R, int C, typename... Pargs>
  size_t count_memory(size_t count, const Eigen::Matrix<var, R, C>& x,
                      const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    count += x.size();
    return count_memory(count, args...);
  }

  template <int R, int C, typename... Pargs>
  size_t count_memory(size_t count, const Eigen::Matrix<double, R, C>& x,
                      const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    return count_memory(count, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const std::vector<var>& x,
                      const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    count += x.size();
    return count_memory(count, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const std::vector<double>& x,
                      const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    return count_memory(count, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const std::vector<int>& x,
                      const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    return count_memory(count, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const var& x, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    count += 1;
    return count_memory(count, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const double& x, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    return count_memory(count, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const int& x, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    offsets_[t] = count;
    return count_memory(count, args...);
  }

  size_t count_memory(size_t count) { return count; }

  /*
   * prepare_x_vis populates x_vis_ with the varis from each of its
   * input arguments. The vari pointers for argument n are copied into x_vis_ at
   * the index starting at offsets_[n]. For Eigen::Matrix types, this copying is
   * done in with column major ordering.
   *
   * Each of the arguments can be an Eigen::Matrix with var or double scalar
   * types, a std::vector with var, double, or int scalar types, or a var, a
   * double, or an int.
   *
   * @tparam R Eigen Matrix row type
   * @tparam C Eigen Matrix column type
   * @tparam Pargs Types of the rest of the arguments to be processed
   * @param x next argument to have its vari pointers copied if necessary
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <int R, int C, typename... Pargs>
  void prepare_x_vis(const Eigen::Matrix<var, R, C>& x, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
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
  void prepare_x_vis(const std::vector<var>& x, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    for (size_t i = 0; i < x.size(); ++i) {
      x_vis_[offsets_[t] + i] = x[i].vi_;
    }
    prepare_x_vis(args...);
  }

  template <typename... Pargs>
  void prepare_x_vis(const std::vector<double>& x, const Pargs&... args) {
    prepare_x_vis(args...);
  }

  template <typename... Pargs>
  void prepare_x_vis(const std::vector<int>& x, const Pargs&... args) {
    prepare_x_vis(args...);
  }

  template <typename... Pargs>
  void prepare_x_vis(const var& x, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    x_vis_[offsets_[t]] = x.vi_;
    prepare_x_vis(args...);
  }

  template <typename... Pargs>
  void prepare_x_vis(const double& x, const Pargs&... args) {
    prepare_x_vis(args...);
  }

  template <typename... Pargs>
  void prepare_x_vis(const int& x, const Pargs&... args) {
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
   * Each of the arguments can be an Eigen::Matrix with var or double scalar
   * types, a std::vector with var, double, or int scalar types, or a var, a
   * double, or an int.
   *
   * As an example of how is_var_ is initialized, take the adj_jac_vari
   * constructed as: adj_jac_vari(a1, a2, a3, a4, a5)
   *
   * If a1 is an Eigen::Matrix<var, 5, 5>, a2 is a double, a3 is a var, a4 is
   * a length 5 vector of vars, and a5 is a length 3 vector of ints,
   * is_var_ is set to {true, false, true, true, false}
   *
   * @param args Input arguments
   */
  explicit adj_jac_vari(const Targs&... args)
      : vari(std::numeric_limits<double>::quiet_NaN()),  // The val_ in this
                                                         // vari is unused
        is_var_({{is_var<typename scalar_type<Targs>::type>::value...}}),
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
   * accumulate_adjoints accumulates, if necessary, the Eigen Matrix of values
   * in its first argument into the adjoints of the varis pointed to by the
   * appropriate elements of x_vis_ and then recursively calls
   * accumulate_adjoints on the rest of the arguments.
   *
   * @tparam R Eigen Matrix row type
   * @tparam C Eigen Matrix column type
   * @tparam Pargs Types of the rest of adjoints to accumulate
   * @param y_adj_jac set of values to be accumulated in adjoints
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <int R, int C, typename... Pargs>
  void accumulate_adjoints(const Eigen::Matrix<double, R, C>& y_adj_jac,
                           const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    if (is_var_[t]) {
      for (int n = 0; n < y_adj_jac.size(); ++n) {
        x_vis_[offsets_[t] + n]->adj_ += y_adj_jac(n);
      }
    }

    accumulate_adjoints(args...);
  }

  /*
   * accumulate_adjoints accumulates, if necessary, the std::vector of values in
   * its first argument into the adjoints of the varis pointed to by the
   * appropriate elements of x_vis_ and then recursively calls
   * accumulate_adjoints on the rest of the arguments.
   *
   * @tparam Pargs Types of the rest of adjoints to accumulate
   * @param y_adj_jac set of values to be accumulated in adjoints
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <typename... Pargs>
  void accumulate_adjoints(const std::vector<double>& y_adj_jac,
                           const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    if (is_var_[t]) {
      for (size_t n = 0; n < y_adj_jac.size(); ++n)
        x_vis_[offsets_[t] + n]->adj_ += y_adj_jac[n];
    }

    accumulate_adjoints(args...);
  }

  /*
   * There are no adjoints to accumulate for std::vector<int> arguments, so
   * accumulate_adjoints simply recursively calls itself on the rest of the
   * arguments.
   *
   * @tparam Pargs Types of the rest of adjoints to accumulate
   * @param y_adj_jac ignored
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <typename... Pargs>
  void accumulate_adjoints(const std::vector<int>& y_adj_jac,
                           const Pargs&... args) {
    accumulate_adjoints(args...);
  }

  /*
   * accumulate_adjoints accumulates, if necessary, the y_adj_jac into the
   * adjoint of vari pointed to by the appropriate element of x_vis_ and then
   * recursively calls accumulate_adjoints on the rest of the arguments.
   *
   * @tparam Pargs Types of the rest of adjoints to accumulate
   * @param y_adj_jac next set of adjoints to be accumulated
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <typename... Pargs>
  void accumulate_adjoints(const double& y_adj_jac, const Pargs&... args) {
    static constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    if (is_var_[t]) {
      x_vis_[offsets_[t]]->adj_ += y_adj_jac;
    }

    accumulate_adjoints(args...);
  }

  /*
   * There are no adjoints to accumulate for an int argument, so
   * accumulate_adjoints simply recursively calls itself on the rest of the
   * arguments.
   *
   * @tparam Pargs Types of the rest of adjoints to accumulate
   * @param y_adj_jac ignored
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <typename... Pargs>
  void accumulate_adjoints(const int& y_adj_jac, const Pargs&... args) {
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
