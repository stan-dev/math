#ifndef STAN_MATH_REV_FUNCTOR_ADJ_OP_HPP
#define STAN_MATH_REV_FUNCTOR_ADJ_OP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <vector>

namespace stan {
namespace math {


/**
 * Stores the values for the functors passed to adj_jac_vari. Passing a var or
 *  or container of vars will store a pointer to the val. If a non-var
 *  type is passed this class will be empty.
 * @tparam T Type to store, for `var` stores a double and containers stores
 * a pointer to doubles.
 * @tparam IsVar A boolean that the user can pass to override the default
 *  behavior of storing data for var types and having an empty class for
 *  non-var types.
 */
template <typename T, bool IsVar = is_var<value_type_t<T>>::value,
          typename = void>
class adj_op;

/**
 * For cases where `IsVar` is true `adj_op` is allocates memory on the stack
 *  allocator and
 */
template <typename T>
class adj_op<T, true, require_container_t<T>> {
 public:
  double* mem_;  // values to store
  using FReturnType = plain_type_t<decltype(value_of(T()))>;
  // If the type is an `std::vector` this stores an `Eigen::VectorXd`
  using RetType = std::conditional_t<is_std_vector<FReturnType>::value,
                                     Eigen::Matrix<double, -1, 1>, FReturnType>;
  // The Eigen map type used for storage
  using eigen_map = Eigen::Map<RetType>;
  eigen_map map_;
  /**
   * Allocate unitialized memory of size `n` for a vector
   * @param n The size of memory to allocate on the stack.
   */
  explicit adj_op(size_t n)
      : mem_(ChainableStack::instance_->memalloc_.alloc_array<double>(n)),
        map_(mem_, n) {}
  /**
   * Allocate unitialized memory of size `n * m` for a matrix
   * @param n Number of rows
   * @param m Number of columns
   */
  explicit adj_op(size_t n, size_t m)
      : mem_(ChainableStack::instance_->memalloc_.alloc_array<double>(n * m)),
        map_(mem_, n, m) {}
  /**
   * Allocate and initialize memory from a Eigen matrix of vars.
   * @tparam EigMat An Eigen type that holds a `var` scalar
   * @param x An Eigen object whose val portion of the `var` will be stored
   */
  template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
  explicit adj_op(EigMat&& x)
      : mem_(
          ChainableStack::instance_->memalloc_.alloc_array<double>(x.size())),
        map_(eigen_map(mem_, x.rows(), x.cols()) = value_of(x)) {}

  /**
   * Allocate and initialize memory from a standard vector of vars.
   * @tparam StdVec a standard library vector with a `var` value type.
   * @param x Vector of `var` types whose `val_`'s will be stored
   */
  template <typename StdVec, require_std_vector_vt<is_var, StdVec>* = nullptr>
  explicit adj_op(StdVec&& x)
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
  explicit adj_op(EigMat&& x)
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
  explicit adj_op(StdVec&& x)
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
  inline auto& map() { return map_; }

  /**
   * Getter for map of the stack allocated data
   */
  inline const auto& map() const { return map_; }
  /**
   * Setter for the ith element from the map.
   * @param i Access the ith element.
   */
  inline double& operator()(size_t i) { return map_(i); }

  /**
   * Getter for the element in the ith row and jth from the map.
   * @param i Access the ith element.
   */
  inline const double& operator()(size_t i) const { return map_(i); }

  /**
   * Setter for the element in the ith row and jth column from the map.
   * @param i Row to access.
   * @param j Column to access.
   */
  inline double& operator()(size_t i, size_t j) { return map_(i, j); }

  /**
   * Getter for the element in the ith row and jth column from the map.
   * @param i Row to access.
   * @param j Column to access.
   */
  inline const double& operator()(size_t i, size_t j) const {
    return map_(i, j);
  }
};

/**
 * When `IsVar` is `false` `adj_op` holds an empty map of size `(0,0)` for
 *  containers.
 */
template <typename T>
class adj_op<T, false, require_container_t<T>> {
 public:
  Eigen::Matrix<double, -1, -1> map_{0, 0};
  explicit adj_op(size_t n) {}
  adj_op(size_t n, size_t m) {}
  template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
  explicit adj_op(const EigMat& x) {}

  template <typename StdVec, require_std_vector_t<StdVec>* = nullptr>
  explicit adj_op(const StdVec& x) {}

  inline auto rows() const { return 0; }

  inline auto cols() const { return 0; }

  inline auto size() const { return 0; }

  inline auto& map() {
    throw_domain_error("adj_op", "", "Attempting to Access Empty adj_op!", "");
    return map_;
  }
  inline const auto& map() const { return map_; }
  inline double& operator()(size_t i) {
    throw_domain_error("adj_op", "", "Attempting to Access Empty adj_op!", "");
    return map_(0, 0);
  }

  inline const double& operator()(size_t i) const {
    throw_domain_error("adj_op", "", "Attempting to Access Empty adj_op!", "");
    return map_(0);
  }

  inline double& operator()(size_t i, size_t j) {
    throw_domain_error("adj_op", "", "Attempting to Access Empty adj_op!", "");
    return map_(0, 0);
  }

  inline const double& operator()(size_t i, size_t j) const {
    throw_domain_error("adj_op", "", "Attempting to Access Empty adj_op!", "");
    return map_(0, 0);
  }
};

/**
 * When `IsVar` is `true` `adj_op` allocates a double on the stack allocator
 *  for scalar types.
 */
template <typename T>
class adj_op<T, true, require_stan_scalar_t<T>> {
 public:
  double* mem_{ChainableStack::instance_->memalloc_.alloc_array<double>(1)};
  explicit adj_op(size_t n) {}
  explicit adj_op(size_t n, size_t m) {}
  template <typename S, require_var_t<S>* = nullptr>
  explicit adj_op(const S& x) {
    *mem_ = x.val();
  }
  template <typename S, require_floating_point_t<S>* = nullptr>
  explicit adj_op(const S& x) {
    *mem_ = x;
  }

  inline auto rows() const { return 0; }

  inline auto cols() const { return 0; }

  inline auto size() const { return 1; }

  inline auto& map() { return *mem_; }
  inline const auto& map() const { return *mem_; }
};

/**
 * When `IsVar` is `false` `adj_op` holds a double with a value of zero.
 */
template <typename T>
class adj_op<T, false, require_stan_scalar_t<T>> {
 public:
  double map_{0.0};
  explicit adj_op(size_t n) {}
  adj_op(size_t n, size_t m) {}
  template <typename S, require_var_t<S>* = nullptr>
  explicit adj_op(const S& x) {}
  template <typename S, require_floating_point_t<S>* = nullptr>
  explicit adj_op(const S& x) {}

  inline auto rows() const { return 0; }

  inline auto cols() const { return 0; }

  inline auto size() const { return 0; }

  inline auto& map() {
    throw_domain_error("adj_op", "", "Attempting to Access Empty adj_op!", "");
    return map_;
  }
  inline const auto& map() const {
    throw_domain_error("adj_op", "", "Attempting to Access Empty adj_op!", "");
    return map_;
  }
};

}
}
#endif
