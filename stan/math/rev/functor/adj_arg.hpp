#ifndef STAN_MATH_REV_FUNCTOR_ADJ_ARG_HPP
#define STAN_MATH_REV_FUNCTOR_ADJ_ARG_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T, typename = void>
struct adj_argr;

template <typename T>
struct adj_argr<T, require_container_t<T>> {
  using FunctorReturnType = plain_type_t<decltype(value_of(std::declval<T>()))>;
  // If the type is an `std::vector` this stores an `Eigen::VectorXd`
  using RetType = std::conditional_t<is_std_vector<FunctorReturnType>::value,
                                     Eigen::Matrix<double, -1, 1>, FunctorReturnType>;
  // The Eigen map type used for storage
  using type = Eigen::Map<RetType>;
};

template <typename T>
struct adj_argr<T, require_stan_scalar_t<T>> {
  using type = double;
};

template <typename T>
using adj_arg_t = typename adj_argr<std::decay_t<T>>::type;


template <typename T, bool SaveValue = true>
auto setup_adj_arg(size_t n) {
  double* mem = ChainableStack::instance_->memalloc_.alloc_array<double>(n);
  return adj_arg_t<T>(mem, n);
}

template <typename T, bool SaveValue, require_t<bool_constant<!SaveValue>>* = nullptr>
auto setup_adj_arg(size_t n) {
  return adj_arg_t<T>(nullptr, 0);
}

template <typename T, bool SaveValue = true>
auto setup_adj_arg(size_t n, size_t m) {
  double* mem = ChainableStack::instance_->memalloc_.alloc_array<double>(n * m);
  return adj_arg_t<T>(mem, n, m);
}

template <typename T, bool SaveValue, require_t<bool_constant<!SaveValue>>* = nullptr>
auto setup_adj_arg(size_t n, size_t m) {
  return adj_arg_t<T>(nullptr, 0, 0);
}

template <typename T, bool SaveValue = true>
auto setup_adj_arg() {
  return adj_arg_t<T>();
}


/**
 * Stores the values for the functors passed to adj_jac_vari. Passing a var or
 *  or container of vars will store a pointer to the val. If a non-var
 *  type is passed this class will be empty.
 * @tparam T Type to store, for `var` stores a double and containers stores
 * a pointer to doubles.
 * @tparam SaveValues A boolean that the user can pass to override the default
 *  behavior of storing data for var types and having an empty class for
 *  non-var types.
 */
template <typename T, bool SaveValues = is_var<value_type_t<T>>::value,
          typename = void>
class adj_arg;

/**
 * For cases where `SaveValues` is true `adj_arg` allocates memory on the
 * stack allocator whose data can be accessed through an `Eigen::Map`
 */
template <typename T>
class adj_arg<T, true, require_container_t<T>> {
 public:
  double* mem_;  // values to store
  using FunctorReturnType = plain_type_t<decltype(value_of(std::declval<T>()))>;
  // If the type is an `std::vector` this stores an `Eigen::VectorXd`
  using RetType = std::conditional_t<is_std_vector<FunctorReturnType>::value,
                                     Eigen::Matrix<double, -1, 1>, FunctorReturnType>;
  // The Eigen map type used for storage
  using eigen_map = Eigen::Map<RetType>;
  eigen_map map_;
  static constexpr bool needs_adj{true};
  /**
   * Allocate unitialized memory of size `n` for a vector
   * @param n The size of memory to allocate on the stack.
   */
  explicit adj_arg(size_t n)
      : mem_(ChainableStack::instance_->memalloc_.alloc_array<double>(n)),
        map_(mem_, n) {}
  /**
   * Allocate unitialized memory of size `n * m` for a matrix
   * @param n Number of rows
   * @param m Number of columns
   */
  adj_arg(size_t n, size_t m)
      : mem_(ChainableStack::instance_->memalloc_.alloc_array<double>(n * m)),
        map_(mem_, n, m) {}
  /**
   * Allocate and initialize memory from a Eigen matrix of vars.
   * @tparam EigMat An Eigen type that holds a `var` scalar
   * @param x An Eigen object whose val portion of the `var` will be stored
   */
  template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
  explicit adj_arg(EigMat&& x)
      : mem_(
            ChainableStack::instance_->memalloc_.alloc_array<double>(x.size())),
        map_(eigen_map(mem_, x.rows(), x.cols()) = value_of(x)) {}

  /**
   * Allocate and initialize memory from a standard vector of vars.
   * @tparam StdVec a standard library vector with a `var` value type.
   * @param x Vector of `var` types whose `val_`'s will be stored
   */
  template <typename StdVec, require_std_vector_vt<is_var, StdVec>* = nullptr>
  explicit adj_arg(StdVec&& x)
      : mem_(
            ChainableStack::instance_->memalloc_.alloc_array<double>(x.size())),
        map_(eigen_map(mem_, x.size())
             = eigen_map(value_of(x).data(), x.size())) {}

  /**
   * Allocate and initialize memory from a Eigen matrix of doubles.
   * @tparam EigMat An Eigen type that holds an arithmetic scalar
   * @param x An Eigen object whose scalars are stored
   */
  template <typename EigMat,
            require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
  explicit adj_arg(EigMat&& x)
      : mem_(
            ChainableStack::instance_->memalloc_.alloc_array<double>(x.size())),
        map_(eigen_map(mem_, x.rows(), x.cols()) = x) {}

  /**
   * Allocate and initialize memory from a standard vector of doubles.
   * @tparam StdVec a standard library vector with an arithmetic value type
   * @param x Vector of arithmetic types to store
   */
  template <typename StdVec,
            require_std_vector_vt<std::is_arithmetic, StdVec>* = nullptr>
  explicit adj_arg(StdVec&& x)
      : mem_(
            ChainableStack::instance_->memalloc_.alloc_array<double>(x.size())),
        map_(eigen_map(mem_, x.size()) = eigen_map(x.data(), x.size())) {}

  /**
   * Query number of rows
   */
  inline auto rows() const { return map_.rows(); }

  /**
   * Query number of columns
   */
  inline auto cols() const { return map_.cols(); }

  /**
   * Query size of container
   */
  inline auto size() const { return map_.size(); }

  /**
   * Setter for map of the stack allocated data
   */
  inline auto& arg() { return map_; }

  /**
   * Getter for map of the stack allocated data
   */
  inline const auto& arg() const { return map_; }
  /**
   * Setter for the ith element from the map.
   * @param i Access the ith element.
   */
  inline double& operator()(size_t i) { return map_.coeffRef(i); }

  /**
   * Getter for the element in the ith row and jth from the map.
   * @param i Access the ith element.
   */
  inline const double& operator()(size_t i) const { return map_.coeffRef(i); }

  /**
   * Setter for the element in the ith row and jth column from the map.
   * @param i Row to access.
   * @param j Column to access.
   */
  inline double& operator()(size_t i, size_t j) { return map_.coeffRef(i, j); }

  /**
   * Getter for the element in the ith row and jth column from the map.
   * @param i Row to access.
   * @param j Column to access.
   */
  inline const double& operator()(size_t i, size_t j) const {
    return map_.coeffRef(i, j);
  }
};

/**
 * When `SaveValues` is `false` `adj_arg` holds an empty map of size `(0,0)` for
 *  containers.
 */
template <typename T>
class adj_arg<T, false, require_container_t<T>> {
 public:
  Eigen::Matrix<double, -1, -1> map_{0, 0};
  static constexpr bool needs_adj{false};
  explicit adj_arg(size_t n) {}
  adj_arg(size_t n, size_t m) {}
  template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
  explicit adj_arg(const EigMat& x) {}

  template <typename StdVec, require_std_vector_t<StdVec>* = nullptr>
  explicit adj_arg(const StdVec& x) {}

  inline static constexpr auto rows() { return 0; }

  inline static constexpr auto cols() { return 0; }

  inline static constexpr auto size() { return 0; }
  /**
   * Accessing data for the no-op class throws a domain error.
   * @throw domain error
   */
  inline auto& arg() {
    throw_domain_error("adj_arg", "", "Attempting to Access Empty adj_arg!",
                       "");
    return map_;
  }

  /**
   * Accessing data for the no-op class throws a domain error.
   * @throw domain error
   */
  inline const auto& arg() const { return map_; }
  inline double& operator()(size_t i) {
    throw_domain_error("adj_arg", "", "Attempting to Access Empty adj_arg!",
                       "");
    return map_(0, 0);
  }

  /**
   * Accessing data for the no-op class throws a domain error.
   * @throw domain error
   */
  inline const double& operator()(size_t i) const {
    throw_domain_error("adj_arg", "", "Attempting to Access Empty adj_arg!",
                       "");
    return map_(0);
  }

  /**
   * Accessing data for the no-op class throws a domain error.
   * @throw domain error
   */
  inline double& operator()(size_t i, size_t j) {
    throw_domain_error("adj_arg", "", "Attempting to Access Empty adj_arg!",
                       "");
    return map_(0, 0);
  }

  /**
   * Accessing data for the no-op class throws a domain error.
   * @throw domain error
   */
  inline const double& operator()(size_t i, size_t j) const {
    throw_domain_error("adj_arg", "", "Attempting to Access Empty adj_arg!",
                       "");
    return map_(0, 0);
  }
};

/**
 * When `SaveValues` is `true` `adj_arg` allocates a double on the stack
 * allocator for scalar types.
 */
template <typename T>
class adj_arg<T, true, require_stan_scalar_t<T>> {
 public:
  double* mem_{ChainableStack::instance_->memalloc_.alloc_array<double>(1)};
  static constexpr bool needs_adj{true};
  explicit adj_arg(size_t n) {}
  adj_arg(size_t n, size_t m) {}
  template <typename S, require_var_t<S>* = nullptr>
  explicit adj_arg(const S& x) {
    *mem_ = x.val();
  }
  template <typename S, require_floating_point_t<S>* = nullptr>
  explicit adj_arg(const S& x) {
    *mem_ = x;
  }

  inline static constexpr auto rows() { return 0; }

  inline static constexpr auto cols() { return 0; }

  inline static constexpr auto size() { return 1; }

  inline auto& arg() { return *mem_; }
  inline const auto& arg() const { return *mem_; }
};

/**
 * When `SaveValues` is `false` `adj_arg` holds a double with a value of zero.
 */
template <typename T>
class adj_arg<T, false, require_stan_scalar_t<T>> {
 public:
  double map_{0.0};
  static constexpr bool needs_adj{false};
  explicit adj_arg(size_t n) {}
  adj_arg(size_t n, size_t m) {}
  template <typename S, require_var_t<S>* = nullptr>
  explicit adj_arg(const S& x) {}
  template <typename S, require_floating_point_t<S>* = nullptr>
  explicit adj_arg(const S& x) {}

  inline static constexpr auto rows() { return 0; }

  inline static constexpr auto cols() { return 0; }

  inline static constexpr auto size() { return 0; }

  /**
   * Accessing data for the no-op class throws a domain error.
   * @throw domain error
   */
  inline auto& arg() {
    throw_domain_error("adj_arg", "", "Attempting to Access Empty adj_arg!",
                       "");
    return map_;
  }
  /**
   * Accessing data for the no-op class throws a domain error.
   * @throw domain error
   */
  inline const auto& arg() const {
    throw_domain_error("adj_arg", "", "Attempting to Access Empty adj_arg!",
                       "");
    return map_;
  }
};

}  // namespace math
}  // namespace stan
#endif
