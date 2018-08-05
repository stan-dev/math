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
 * Taken from
 * http://aherrmann.github.io/programming/2016/02/28/unpacking-tuples-in-cpp14/
 */
template <class F, class Tuple, std::size_t... Is>
constexpr auto apply_impl(const F& f, const Tuple& t,
                          std::index_sequence<Is...>) {
  return f(std::get<Is>(t)...);
}

template <class F, class Tuple>
constexpr auto apply(const F& f, const Tuple& t) {
  return apply_impl(f, t, std::make_index_sequence<std::tuple_size<Tuple>{}>{});
}

/**
 * build_y_adj takes the adjoint from the vari pointed to
 * by y_vi_[0] and stores it in y_adj
 *
 * @tparam size dimensionality of M_
 * @param y_vi_ pointer to pointer to vari
 * @param M_ ignored in this specialization
 * @param[out] y_adj reference to variable where adjoint is to be stored
 */
template <size_t size>
void build_y_adj(vari** y_vi_, const std::array<int, size>& M_, double& y_adj) {
  y_adj = y_vi_[0]->adj_;
}

/**
 * build_y_adj takes the adjoints from the varis pointed to
 * by y_vi_ and stores them in y_adj
 *
 * @tparam size dimensionality of M_
 * @param y_vi_ pointer to pointers to varis
 * @param M_ shape of y_adj
 * @param[out] y_adj reference to Eigen::Matrix where adjoints are to be stored
 */
template <size_t size, int R, int C>
void build_y_adj(vari** y_vi_, const std::array<int, size>& M_,
                 Eigen::Matrix<double, R, C>& y_adj) {
  y_adj.resize(M_[0], M_[1]);
  for (int m = 0; m < M_[0] * M_[1]; ++m)
    y_adj(m) = y_vi_[m]->adj_;
}

/**
 * compute the dimensionality of the given template argument. By
 * default assume zero.
 */
template <typename T>
struct compute_dims {
  static constexpr size_t value = 0;
};

/**
 * compute the dimensionality of the given template argument.
 * Eigen::VectorXd is treated like a Matrix and given dimension two
 */
template <>
struct compute_dims<Eigen::VectorXd> {
  static constexpr size_t value = 2;
};

/**
 * compute the dimensionality of the given template argument.
 * Eigen::RowVectorXd is treated like a Matrix and given dimension two
 */
template <>
struct compute_dims<Eigen::RowVectorXd> {
  static constexpr size_t value = 2;
};

/**
 * compute the dimensionality of the given template argument.
 * Eigen::MatrixXd has dimension two
 */
template <>
struct compute_dims<Eigen::MatrixXd> {
  static constexpr size_t value = 2;
};
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
  std::array<bool, sizeof...(Targs)> is_var_;
  using FReturnType
      = std::result_of_t<F(decltype(is_var_), decltype(value_of(Targs()))...)>;

  F f_;
  std::array<int, sizeof...(Targs)> offsets_;
  vari** x_vis_;
  std::array<int, compute_dims<FReturnType>::value> M_;
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

  template <int R, int C>
  size_t count_memory(size_t count, const Eigen::Matrix<var, R, C>& x) {
    constexpr int t = sizeof...(Targs) - 1;
    offsets_[t] = count;
    return count + x.size();
  }

  template <int R, int C>
  size_t count_memory(size_t count, const Eigen::Matrix<double, R, C>& x) {
    constexpr int t = sizeof...(Targs) - 1;
    offsets_[t] = count;
    return count;
  }

  size_t count_memory(size_t count, const var& x) {
    constexpr int t = sizeof...(Targs) - 1;
    offsets_[t] = count;
    return count + 1;
  }

  size_t count_memory(size_t count, const double& x) {
    constexpr int t = sizeof...(Targs) - 1;
    offsets_[t] = count;
    return count;
  }

  /*
   * prepare_x_vis recursively populates x_vis_ with the varis from each of its
   * input arguments. The vari pointers for argument n are copied into x_vis_ at
   * the index starting at offsets_[n]. The array is_var_ is populated with
   * whether or not the scalar type in each argument includes a var
   *
   * Each of the arguments can be an Eigen::Matrix<var, R, C>, an
   * Eigen::Matrix<double, R, C>, a var, or a double.
   *
   * As an example, take the call:
   *    prepare_x_vis(a1, a2, a3)
   *
   * If a1 is an Eigen::Matrix<var, 5, 5>, a2 is a double, and a3 is a var, the
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
    is_var_[t] = true;
    for (int i = 0; i < x.size(); ++i) {
      x_vis_[offsets_[t] + i] = x(i).vi_;
    }
    prepare_x_vis(args...);
  }

  template <int R, int C, typename... Pargs>
  void prepare_x_vis(const Eigen::Matrix<double, R, C>& x,
                     const Pargs&... args) {
    constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    is_var_[t] = false;
    prepare_x_vis(args...);
  }

  template <typename... Pargs>
  void prepare_x_vis(const var& x, const Pargs&... args) {
    constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    is_var_[t] = true;
    x_vis_[offsets_[t]] = x.vi_;
    prepare_x_vis(args...);
  }

  template <typename... Pargs>
  void prepare_x_vis(const double& x, const Pargs&... args) {
    constexpr int t = sizeof...(Targs) - sizeof...(Pargs) - 1;
    is_var_[t] = false;
    prepare_x_vis(args...);
  }

  template <int R, int C>
  void prepare_x_vis(const Eigen::Matrix<var, R, C>& x) {
    constexpr int t = sizeof...(Targs) - 1;
    is_var_[t] = true;
    for (int i = 0; i < x.size(); ++i) {
      x_vis_[offsets_[t] + i] = x(i).vi_;
    }
  }

  template <int R, int C>
  void prepare_x_vis(const Eigen::Matrix<double, R, C>& x) {
    constexpr int t = sizeof...(Targs) - 1;
    is_var_[t] = false;
  }

  void prepare_x_vis(const var& x) {
    constexpr int t = sizeof...(Targs) - 1;
    is_var_[t] = true;
    x_vis_[offsets_[t]] = x.vi_;
  }

  void prepare_x_vis(const double& x) {
    constexpr int t = sizeof...(Targs) - 1;
    is_var_[t] = false;
  }
  explicit adj_jac_vari()
      : vari(std::numeric_limits<double>::quiet_NaN()),  // The val_ in this
                                                         // vari is unused
        x_vis_(NULL),
        y_vi_(NULL) {}

  /**
   * build_return_varis_and_vars takes the double output of the
   * F::operator(), allocates a vari with the
   * value of the output, and then builds and returns a var pointing
   * at that vari
   *
   * @param val_y output of F::operator()
   * @return Eigen::Matrix of vars
   */
  var build_return_varis_and_vars(const double& val_y) {
    y_vi_ = ChainableStack::instance().memalloc_.alloc_array<vari*>(1);
    y_vi_[0] = new vari(val_y, false);

    return y_vi_[0];
  }

  /**
   * build_return_varis_and_vars takes the Eigen::Matrix output of the
   * F::operator(), allocates and populates an array of varis with the
   * values of the output, and then builds an Eigen::Matrix of vars
   * of the same shape as the original output and assigns the vari pointers
   * in those vars to the allocated varis
   *
   * @param val_y output of F::operator()
   * @return Eigen::Matrix of vars
   */
  Eigen::Matrix<var, Eigen::Dynamic, 1> build_return_varis_and_vars(
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& val_y) {
    M_[0] = val_y.rows();
    M_[1] = val_y.cols();
    Eigen::Matrix<var, Eigen::Dynamic, 1> var_y(M_[0] * M_[1]);

    y_vi_
        = ChainableStack::instance().memalloc_.alloc_array<vari*>(var_y.size());
    for (int m = 0; m < var_y.size(); ++m) {
      y_vi_[m] = new vari(val_y(m), false);
      var_y(m) = y_vi_[m];
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
   *  4. Calls the appropriate build_return_varis_and_vars to construct the
   * appropriate output data structure of vars
   *
   * The input argument types can be any mix of double, var, or Eigen::Matrices
   * with double or var scalar components
   *
   * @param args Input arguments
   * @return Output of f_ as vars
   */
  auto operator()(const Targs&... args) {
    x_vis_ = ChainableStack::instance().memalloc_.alloc_array<vari*>(
        count_memory(0, args...));
    prepare_x_vis(args...);

    return build_return_varis_and_vars(f_(is_var_, value_of(args)...));
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

  template <int R, int C>
  void accumulate_adjoints(const Eigen::Matrix<double, R, C>& y_adj_jac) {
    const int t = sizeof...(Targs) - 1;
    if (is_var_[t]) {
      for (int n = 0; n < y_adj_jac.size(); ++n) {
        x_vis_[offsets_[t] + n]->adj_ += y_adj_jac(n);
      }
    }
  }

  void accumulate_adjoints(const double& y_adj_jac) {
    const int t = sizeof...(Targs) - 1;
    if (is_var_[t]) {
      x_vis_[offsets_[t]]->adj_ += y_adj_jac;
    }
  }

  /**
   * chain propagates the adjoints at the output varis (y_vi_) back to the input
   * varis (x_vis_) by packing the adjoings in an appropriate container
   * using build_y_adj and then using the multiply_adjoint_jacobian function of
   * the user defined functor to compute what the adjoints on x_vis_ should be
   *
   * Unlike the constructor, this operation may be called multiple times during
   * the life of the vari
   */
  void chain() {
    FReturnType y_adj;
    build_y_adj(y_vi_, M_, y_adj);
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
 * (optional) a nullary constructor
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
auto adj_jac_apply(const Targs&... args) {
  auto vi = new adj_jac_vari<F, Targs...>();

  return (*vi)(args...);
}

}  // namespace math
}  // namespace stan
#endif
