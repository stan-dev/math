#ifndef STAN_MATH_REV_FUN_ADJOINT_OF_HPP
#define STAN_MATH_REV_FUN_ADJOINT_OF_HPP

#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

namespace internal {
struct nonexisting_adjoint {
  template <typename T>
  nonexisting_adjoint operator+(const T&) {
    return *this;
  }
  template <typename T>
  nonexisting_adjoint operator+=(T) const {
    throw std::runtime_error(
        "internal::nonexisting_adjoint::operator+= should never be called! "
        "Please file a bug report.");
  }
  template <typename T>
  nonexisting_adjoint operator-=(T) const {
    throw std::runtime_error(
        "internal::nonexisting_adjoint::operator-= should never be called! "
        "Please file a bug report.");
  }

  static inline nonexisting_adjoint array() {
    throw std::runtime_error(
        "internal::nonexisting_adjoint.array() should never be called! "
        "Please file a bug report.");
  }
};
}  // namespace internal

/**
 * Returns a reference to a variable's adjoint.
 *
 * @param x a var
 * @return reference to `x`'s adjoint
 */
template <typename T, require_var_t<T>* = nullptr>
auto& adjoint_of(const T& x) {
  return x.adj();
}

template <typename T, require_var_t<T>* = nullptr>
auto& get_adj(const T& x) {
  return x.adj();
}

template <typename T, require_eigen_vt<is_var, T>* = nullptr>
auto get_adj(const T& x) {
  return x.adj();
}

template <typename T, require_st_var<T>* = nullptr,
          require_std_vector_t<T>* = nullptr>
auto get_adj(const T& x) {
  std::vector<promote_scalar_t<double, value_type_t<T>>> res(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    res[i] = get_adj(x[i]);
  }
  return res;
}

/**
 * Returns a reference to a variable's adjoint. If the input object is not var,
 * it does not have an adjoint and this returns a dummy object. It defines
 * operators += and -=, but they should not actually be called.
 *
 * @param x any non-var object
 * @return a dummy adjoint
 */
template <typename T, require_not_var_t<T>* = nullptr>
internal::nonexisting_adjoint adjoint_of(const T& x) {
  return {};
}

}  // namespace math
}  // namespace stan

#endif  // ADJOINT_OF_HPP
