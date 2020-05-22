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
namespace internal {

/**
 * Place the data for the y_adj_ into the appropriate container.
 */
template <typename T, typename = void>
struct build_y_adj;

/**
 * Specialization returns the real adjoint value for arithmetic and eigen
 * types. Since `vari_value<Eigen>`` holds the dimensions we don't need M.
 */
template <typename T>
struct build_y_adj<T,
                   require_t<disjunction<std::is_arithmetic<T>, is_eigen<T>>>> {
  /**
   * Retreive the adjoints for y_adj_
   *
   * @tparam size dimensionality of M
   * @param[in] y_vi pointer to pointer to vari
   * @param[in] M shape of y_adj (not used in this specialization)
   */
  template <typename VariType, size_t size>
  static inline auto& apply(VariType&& y_vi, const std::array<int, size>& M) {
    return y_vi->adj_;
  }
};

/**
 * Specialization returns the real adjoint value for std vector return types.
 */
template <typename T>
struct build_y_adj<T, require_std_vector_t<T>> {
  /**
   * Store the adjoints from y_vi in y_adj
   *
   * @tparam size dimensionality of M
   * @param[in] y_vi pointer to pointers to varis
   * @param[in] M shape of y_adj
   */
  template <typename VariType, size_t size>
  static inline T apply(VariType& y_vi, const std::array<int, size>& M) {
    T y_adj;
    y_adj.reserve(M[0]);
    for (size_t m = 0; m < M[0]; ++m) {
      y_adj.emplace_back(y_vi[m]->adj_);
    }
    return y_adj;
  }
};

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
  using x_vis_tuple_ = var_to_vari_filter_t<std::decay_t<Targs>...>;
  using x_vis_size_ = std::tuple_size<x_vis_tuple_>;
  x_vis_tuple_ x_vis_{};  // tuple holding pointers to input are var mem.

  // Given a parameter pack, count the number that have a var type.
  template <typename... Types>
  using remaining_vars_ = std::tuple_size<
      stan::math::var_to_vari_filter_t<std::decay_t<Types>...>>;

  /**
   * Given the stored vari types, calculate the position of the element in
   * the tuple holding pointers to vari that is associated with in input var.
   */
  template <typename... Types>
  using var_position_ = std::integral_constant<
      size_t,
      x_vis_size_::value - remaining_vars_<std::decay_t<Types>...>::value>;

  x_vis_alloc(const Targs&... args) : vari(NOT_A_NUMBER) {
    fill_adj_jac(x_vis_, args...);
  }

  /**
   * Fill out the x_vis_ argument for each input that is a var.
   */
  template <typename Mem>
  void fill_adj_jac(Mem&) {}
  /**
   * Specialization that is a no-op for arithmetic class argument types.
   * @tparam Mem tuple of vari pointers associated with the adjoint if it's
   * was originally a vari type.
   * @tparam Arith An arithmetic type.
   * @tparam Pargs parameter pack passed along recursivly.
   * @param mem tuple of vari pointers.
   * @param x Arithmetic.
   * @param args forwarded to recursive calls.
   */
  template <typename Mem, typename Arith, typename... Pargs,
            require_st_arithmetic<Arith>* = nullptr>
  void fill_adj_jac(Mem& mem, Arith&& x, Pargs&&... args) {
    fill_adj_jac(mem, args...);
  }

  /**
   * Specialization for `var_value` types.
   * @tparam Mem tuple of vari pointers associated with the adjoint
   * @tparam VarValue Type A `var_value<T>`.
   * @tparam Pargs parameter pack passed along recursivly.
   * @param mem tuple of vari pointers.
   * @param x var_value<T> type.
   * @param args forwarded to recursive calls.
   */
  template <typename Mem, typename VarValue, typename... Pargs,
            require_var_value_t<VarValue>* = nullptr>
  void fill_adj_jac(Mem& mem, VarValue&& x, Pargs&&... args) {
    using ::stan::internal::get_var_vari_value_t;
    using vari_type = get_var_vari_value_t<VarValue>;
    static constexpr size_t t = var_position_<VarValue, Pargs...>::value;
    std::get<t>(mem)
        = ChainableStack::instance_->memalloc_.alloc_array<vari_type>(1);
    std::get<t>(mem) = x.vi_;
    fill_adj_jac(mem, args...);
  }
  // std::vector<container> assumes not ragged
  template <typename Mem, typename Vec, typename... Pargs,
            require_std_vector_vt<is_container, Vec>* = nullptr>
  void fill_adj_jac(Mem& mem, Vec&& x, Pargs&&... args) {
    using ::stan::internal::get_var_vari_value_t;
    using vari_type = get_var_vari_value_t<scalar_type_t<Vec>>;
    static constexpr size_t t = var_position_<Vec, Pargs...>::value;
    std::get<t>(mem)
        = ChainableStack::instance_->memalloc_.alloc_array<vari_type>(
            x.size() * x[0].size());
    save_varis(std::get<t>(mem), x);
    fill_adj_jac(mem, args...);
  }

  // std::vector<var_value>
  template <typename Mem, typename Vec, typename... Pargs,
            require_std_vector_vt<is_var_value, Vec>* = nullptr>
  void fill_adj_jac(Mem& mem, Vec&& x, Pargs&&... args) {
    using ::stan::internal::get_var_vari_value_t;
    using vari_type = get_var_vari_value_t<value_type_t<Vec>>;
    static constexpr size_t t = var_position_<Vec, Pargs...>::value;
    std::get<t>(mem)
        = ChainableStack::instance_->memalloc_.alloc_array<vari_type*>(
            x.size());
    auto&& inner_mem = std::get<t>(mem);
    for (int i = 0; i < x.size(); ++i) {
      inner_mem[i] = x[i].vi_;
    }
    fill_adj_jac(mem, args...);
  }
  // Eigen
  template <typename Mem, typename EigMat, typename... Pargs,
            require_eigen_vt<is_var_value, EigMat>* = nullptr>
  void fill_adj_jac(Mem& mem, EigMat&& x, Pargs&&... args) {
    using ::stan::internal::get_var_vari_value_t;
    using vari_type = get_var_vari_value_t<EigMat>;
    static constexpr size_t t = var_position_<EigMat, Pargs...>::value;
    auto&& local_mem = std::get<t>(mem);
    local_mem = new vari_type(x.val().eval());
    local_mem->adj_ = x.adj().eval();
    using plain_obj = typename vari_type::PlainObject;
    ChainableStack::instance_->var_stack_.push_back(
        new static_to_dynamic_vari<plain_obj>(x.data()[0].vi_, std::get<t>(mem),
                                              x.size()));
    fill_adj_jac(mem, args...);
  }
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
  using FReturnType
      = std::result_of_t<F(decltype(is_var_), decltype(value_of(Targs()))...)>;
  // Given a parameter pack, count the number that have a var type.
  template <typename... Types>
  using remaining_vars_ = std::tuple_size<
      stan::math::var_to_vari_filter_t<std::decay_t<Types>...>>;

  using x_vis_tuple_ = var_to_vari_filter_t<std::decay_t<Targs>...>;
  using x_vis_size_ = std::tuple_size<x_vis_tuple_>;

  // enable functions if their primary argument is a var.
  template <size_t TargSize, size_t PargSize>
  using require_arg_var_t = std::enable_if_t<is_var_[TargSize - PargSize - 1]>;
  template <size_t TargSize, size_t PargSize>
  // enable functions if their primary argument is not a var.
  using require_not_arg_var_t
      = std::enable_if_t<!is_var_[TargSize - PargSize - 1]>;

  /**
   * Given the stored vari types, calculate the position of the element in
   * the tuple holding pointers to vari that is associated with in input var.
   */
  template <typename... Types>
  using var_position_ = std::integral_constant<
      size_t,
      x_vis_size_::value - remaining_vars_<std::decay_t<Types>...>::value>;

  /**
   * y_vi will be a pointer for eigen and arithmetic types but a pointer to
   * pointers for std::vectors.
   */
  template <typename T>
  using y_vi_type_t
      = std::conditional_t<is_std_vector<std::decay_t<T>>::value,
                           vari_value<value_type_t<T>>**, vari_value<T>*>;

  F f_;  // Function to be invoked
  x_vis_alloc<Targs...>* x_vis_alloc_;

  x_vis_tuple_& x_vis_{x_vis_alloc_->x_vis_};
  // dimensions of output matrix.
  std::array<int, internal::compute_dims<FReturnType>::value> M_;

  y_vi_type_t<FReturnType> y_vi_;  // vari pointer for output.

  /**
   * Initializes is_var_ with true if the scalar type in each argument
   * is a var (and false if not)
   */
  adj_jac_vari()
      : vari(NOT_A_NUMBER),  // The val_ in this vari is unused
        y_vi_(nullptr) {}

  adj_jac_vari(x_vis_alloc<Targs...>* x)
      : vari(NOT_A_NUMBER), x_vis_alloc_(x), y_vi_(nullptr){};

  /**
   * Return a var with a new vari holding the given value
   *
   * @param val_y output of F::operator()
   * @return var
   */
  template <typename T, require_arithmetic_t<T>* = nullptr>
  inline auto build_return_varis_and_vars(const T& val_y) {
    using vari_type = vari_value<T>;
    using var_type = var_value<T>;
    y_vi_ = new vari_type(val_y);
    return var_type(y_vi_);
  }

  /**
   * Return a std::vector of vars created from newly allocated varis initialized
   * with the values of val_y
   *
   * @param val_y output of F::operator()
   * @return std::vector of vars
   */
  template <typename T>
  inline std::vector<var_value<T>> build_return_varis_and_vars(
      const std::vector<T>& val_y) {
    using vari_type = vari_value<T>;
    using var_type = var_value<T>;
    M_[0] = val_y.size();
    std::vector<var_type> var_y;
    var_y.reserve(M_[0]);
    y_vi_ = ChainableStack::instance_->memalloc_.alloc_array<vari_type*>(M_[0]);
    for (size_t m = 0; m < M_[0]; ++m) {
      y_vi_[m] = new vari_type(val_y[m], false);
      var_y.emplace_back(y_vi_[m]);
    }
    return var_y;
  }

  /**
   * Return an Eigen::Matrix of vars created from newly allocated varis
   * initialized with the values of val_y. The shape of the new matrix comes
   * from M_
   *
   * @tparam R number of rows, can be Eigen::Dynamic
   * @tparam C number of columns, can be Eigen::Dynamic
   * @param val_y output of F::operator()
   * @return Eigen::Matrix of vars
   */
  template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
  inline auto build_return_varis_and_vars(EigMat&& val_y) {
    using eig_type = typename std::decay_t<EigMat>::PlainObject;
    using var_type = var_value<eig_type>;
    using vari_type = vari_value<eig_type>;
    M_[0] = val_y.rows();
    M_[1] = val_y.cols();
    y_vi_ = new vari_type(val_y);
    return var_type(y_vi_);
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
    auto func_ret = f_(is_var_, value_of(args)...);
    auto return_obj = build_return_varis_and_vars(func_ret);
    return return_obj;
  }

  template <
      typename Mem, typename T, typename... Pargs,
      require_not_arg_var_t<sizeof...(Targs), sizeof...(Pargs)>* = nullptr>
  inline void accumulate_adjoints_in_varis(Mem& varis, T&& x, Pargs&&... args) {
    accumulate_adjoints_in_varis(varis, args...);
  }
  /**
   * Accumulate, if necessary, the values of y_adj_jac into the
   * adjoints of the varis pointed to by the appropriate elements
   * of x_vis_. Recursively calls accumulate_adjoints_in_varis on the rest of
   * the arguments.
   *
   * @tparam R number of rows, can be Eigen::Dynamic
   * @tparam C number of columns, can be Eigen::Dynamic
   * @tparam Pargs Types of the rest of adjoints to accumulate
   *
   * @param y_adj_jac set of values to be accumulated in adjoints
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <typename Mem, typename EigMat, typename... Pargs,
            require_eigen_t<EigMat>* = nullptr,
            require_arg_var_t<sizeof...(Targs), sizeof...(Pargs)>* = nullptr>
  inline void accumulate_adjoints_in_varis(Mem& varis, const EigMat& y_adj_jac,
                                           const Pargs&... args) {
    static constexpr size_t position = sizeof...(Targs) - sizeof...(Pargs) - 1;
    static constexpr size_t ind = compile_time_accumulator(is_var_, position);
    static constexpr size_t t = x_vis_size_::value - ind;
    auto&& local_mem = std::get<t>(varis);
    local_mem->adj_ += y_adj_jac;
    accumulate_adjoints_in_varis(varis, args...);
  }

  /**
   * Accumulate, if necessary, the values of y_adj_jac into the
   * adjoints of the varis pointed to by the appropriate elements
   * of x_vis_. Recursively calls accumulate_adjoints_in_varis on the rest of
   * the arguments.
   *
   * @tparam Pargs Types of the rest of adjoints to accumulate
   * @param y_adj_jac set of values to be accumulated in adjoints
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <typename Mem, typename... Pargs,
            require_arg_var_t<sizeof...(Targs), sizeof...(Pargs)>* = nullptr>
  inline void accumulate_adjoints_in_varis(Mem& varis,
                                           const std::vector<double>& y_adj_jac,
                                           const Pargs&... args) {
    static constexpr size_t position = sizeof...(Targs) - sizeof...(Pargs) - 1;
    static constexpr size_t ind = compile_time_accumulator(is_var_, position);
    static constexpr size_t t = x_vis_size_::value - ind;
    for (int n = 0; n < y_adj_jac.size(); ++n) {
      std::get<t>(varis)[n]->adj_ += y_adj_jac[n];
    }
    accumulate_adjoints_in_varis(varis, args...);
  }

  /**
   * Accumulate, if necessary, the value of y_adj_jac into the
   * adjoint of the vari pointed to by the appropriate element
   * of x_vis_. Recursively calls accumulate_adjoints_in_varis on the rest of
   * the arguments.
   *
   * @tparam Pargs Types of the rest of adjoints to accumulate
   * @param y_adj_jac next set of adjoints to be accumulated
   * @param args the rest of the arguments (that will be iterated through
   * recursively)
   */
  template <typename Mem, typename... Pargs,
            require_arg_var_t<sizeof...(Targs), sizeof...(Pargs)>* = nullptr>
  inline void accumulate_adjoints_in_varis(Mem& varis, const double& y_adj_jac,
                                           const Pargs&... args) {
    static constexpr size_t position = sizeof...(Targs) - sizeof...(Pargs) - 1;
    static constexpr size_t ind = compile_time_accumulator(is_var_, position);
    static constexpr size_t t = x_vis_size_::value - ind;
    std::get<t>(varis)->adj_ += y_adj_jac;
    accumulate_adjoints_in_varis(varis, args...);
  }

  template <typename Mem>
  inline void accumulate_adjoints_in_varis(Mem&) {}

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
    FReturnType blahh = internal::build_y_adj<FReturnType>::apply(y_vi_, M_);
    auto y_adj_jacs = f_.multiply_adjoint_jacobian(is_var_, blahh);

    apply(
        [&, this](auto&&... args) {
          this->accumulate_adjoints_in_varis(x_vis_, args...);
        },
        y_adj_jacs);
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
  auto* x_alloc = new x_vis_alloc<Targs...>(args...);
  auto vi = new adj_jac_vari<F, Targs...>(x_alloc);

  return (*vi)(args...);
}

}  // namespace math
}  // namespace stan
#endif
