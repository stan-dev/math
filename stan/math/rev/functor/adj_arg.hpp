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

/**
 * Deduces the type of data members needed in functors for `adj_jac_apply`.
 * @tparam Type used to deduce this traits `type`. For a scalar this will be
 * a double and for containers this will be an `Eigen::Map` of either the
 * plain type for Eigen types or an `Eigen::VectorXd` for `std::vector`.
 */
template <typename T, typename = void>
struct adj_arg;

/**
 * Specialization for container types of data members needed in functors
 *  for `adj_jac_apply`.
 * @tparam Type used to deduce this traits `type`. This will be an `Eigen::Map`
 *  of either the plain type for Eigen types or an `Eigen::VectorXd` for
 *  `std::vector`.
 */
template <typename T>
struct adj_arg<T, require_container_t<T>> {
  using FunctorReturnType = plain_type_t<decltype(value_of(std::declval<T>()))>;
  // If the type is an `std::vector` this stores an `Eigen::VectorXd`
  using RetType
      = std::conditional_t<is_std_vector<FunctorReturnType>::value,
                           Eigen::Matrix<double, -1, 1>, FunctorReturnType>;
  // The Eigen map type used for storage
  using type = Eigen::Map<RetType>;
  using const_type = Eigen::Map<const RetType>;
};

/**
 * Specialization for scalar types of data members needed in functors
 *  for `adj_jac_apply`.
 * @tparam Type used to deduce this traits `type`. For a scalar this will be
 * a double.
 */
template <typename T>
struct adj_arg<T, require_stan_scalar_t<T>> {
  using type = double;
  using const_type = const double;
};

/**
 * Deduces the type of data members needed in functors for `adj_jac_apply`.
 * @tparam Type used to deduce this traits `type`. For a scalar this will be
 * a double and for containers this will be an `Eigen::Map` of either the
 * plain type for Eigen types or an `Eigen::VectorXd` for `std::vector`.
 */
template <typename T>
using adj_arg_t = typename adj_arg<std::decay_t<T>>::type;

/**
 * Deduces type for writing to data members for functors used by `adj_jac_apply`
 * @tparam Type used to deduce this traits `type`. For a scalar this will be
 * a `const double` and for containers this will be an `Eigen::Map` of either the
 * const plain type for Eigen types or an `const Eigen::VectorXd` for
 * `std::vector`.
 */
template <typename T>
using adj_arg_reader_t = typename adj_arg<std::decay_t<T>>::const_type;

/**
 * Creates memory on Stan's stack allocator for an associated data member
 * in functors passed to `adj_jac_apply`
 * @tparam T User provided type for deducing the return from `adj_arg_t`.
 * @tparam SaveValue An optional boolean for deducing whether the data member
 *  needs to allocate memory. When true will allocate memory and when false will
 *  pass either a zero size `Eigen::Map` or a default initilized scalar.
 * @param n Size of allocation
 */
template <typename T, bool SaveValue = is_var<scalar_type_t<T>>::value,
 require_container_t<T>* = nullptr,
 require_t<bool_constant<SaveValue>>* = nullptr>
inline adj_arg_t<T> setup_adj_arg(Eigen::Index n) {
  double* mem = ChainableStack::instance_->memalloc_.alloc_array<double>(n);
  return {mem, n};
}

/**
 * Overload of `setup_adj_arg` that does not allocate memory for containers.
 * @tparam T User provided type for deducing the return from `adj_arg_t`.
 * @tparam SaveValue An optional boolean for deducing whether the data member
 *  needs to allocate memory. When true will allocate memory and when false will
 *  pass either a zero size `Eigen::Map` or a default initilized scalar.
 */
template <typename T, bool SaveValue,
  require_container_t<T>* = nullptr,
  require_t<bool_constant<!SaveValue>>* = nullptr>
inline adj_arg_t<T> setup_adj_arg(Eigen::Index /* n */) {
  return {nullptr, 0};
}

/**
 * Overload for scalars
 * @tparam T User provided type for deducing the return from `adj_arg_t`.
 * @tparam SaveValue An optional boolean for deducing whether the data member
 *  needs to allocate memory. When true will allocate memory and when false will
 *  pass either a zero size `Eigen::Map` or a default initilized scalar.
 */
template <typename T, bool SaveValue = is_var<scalar_type_t<T>>::value,
 require_stan_scalar_t<T>* = nullptr>
inline adj_arg_t<T> setup_adj_arg(Eigen::Index /* n */) {
  return {};
}

/**
 * Creates memory on Stan's stack allocator for an associated data member
 * in functors passed to `adj_jac_apply`
 * @tparam T User provided type for deducing the return from `adj_arg_t`.
 * @tparam SaveValue An optional boolean for deducing whether the data member
 *  needs to allocate memory. When true will allocate memory and when false will
 *  pass either a zero size `Eigen::Map` or a default initilized scalar.
 * @param n Size of allocation
 */
template <typename T, bool SaveValue = is_var<scalar_type_t<T>>::value,
  require_container_t<T>* = nullptr,
  require_t<bool_constant<SaveValue>>* = nullptr>
inline adj_arg_t<T> setup_adj_arg(Eigen::Index n, Eigen::Index m) {
  double* mem = ChainableStack::instance_->memalloc_.alloc_array<double>(n * m);
  return {mem, n, m};
}

/**
 * Overload of `setup_adj_arg` that does not allocate memory.
 * @tparam T User provided type for deducing the return from `adj_arg_t`.
 * @tparam SaveValue An optional boolean for deducing whether the data member
 *  needs to allocate memory. When true will allocate memory and when false will
 *  pass either a zero size `Eigen::Map` or a default initilized scalar.
 */
template <typename T, bool SaveValue = is_var<scalar_type_t<T>>::value,
          require_container_t<T>* = nullptr,
          require_t<bool_constant<!SaveValue>>* = nullptr>
inline adj_arg_t<T> setup_adj_arg(Eigen::Index /* n */, Eigen::Index /* m */) {
  return {nullptr, 0, 0};
}


/**
 * Overload of `setup_adj_arg` for scalars
 * @tparam T User provided type for deducing the return from `adj_arg_t`.
 * @tparam SaveValue An optional boolean for deducing whether the data member
 *  needs to allocate memory. When true will allocate memory and when false will
 *  pass either a zero size `Eigen::Map` or a default initilized scalar.
 */
template <typename T, bool SaveValue = is_var<scalar_type_t<T>>::value,
  require_stan_scalar_t<T>* = nullptr>
inline adj_arg_t<T> setup_adj_arg(Eigen::Index /* n */, Eigen::Index /* m */) {
  return {};
}

/**
 * Overload of `setup_adj_arg` for that does not allocate memory.
 * @tparam T User provided type for deducing the return from `adj_arg_t`.
 * @tparam SaveValue An optional boolean for deducing whether the data member
 *  needs to allocate memory. When true will allocate memory and when false will
 *  pass either a zero size `Eigen::Map` or a default initilized scalar.
 */
template <typename T, bool SaveValue = true>
inline adj_arg_t<T> setup_adj_arg() {
  return {};
}

}  // namespace math
}  // namespace stan
#endif
