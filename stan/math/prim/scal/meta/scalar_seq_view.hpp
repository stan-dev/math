#ifndef STAN_MATH_PRIM_SCAL_META_SCALAR_SEQ_VIEW_HPP
#define STAN_MATH_PRIM_SCAL_META_SCALAR_SEQ_VIEW_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <type_traits>
#include <utility>

namespace stan {
/**
 * scalar_seq_view provides a uniform sequence-like wrapper around either a
 * scalar or a sequence of scalars.
 *
 * @tparam C the container type; will be the scalar type if wrapping a scalar
 * @tparam T the scalar type
 */
template <typename C, typename T = typename scalar_type<C>::type>
class scalar_seq_view {
 public:
  // Alias to check perf forwarded param and declared type are the same.
  template <typename Type1, typename Type2>
  using enable_if_same_container = std::enable_if_t<
      std::is_same<std::decay_t<Type1>, std::decay_t<Type2>>::value>;

  template <typename K, typename = enable_if_same_container<C, K>>
  explicit scalar_seq_view(K&& c) : c_(std::forward<K>(c)) {}

  /**
   * Segfaults if out of bounds.
   * @param i index
   * @return the element at the specified position in the container
   */
  const auto&& operator[](int i) const { return c_[i]; }
  auto&& operator[](int i) { return c_[i]; }

  const int size() const { return c_.size(); }
  int size() { return c_.size(); }

 private:
  C c_;
};

/**
 * This specialization handles wrapping a scalar as if it were a sequence.
 *
 * @tparam T the scalar type
 */
template <typename T>
class scalar_seq_view<T, T> {
 public:
  explicit scalar_seq_view(const T& t) : t_(t) {}

  const T& operator[](int /* i */) const { return t_; }
  T& operator[](int /* i */) { return t_; }

  int size() const { return 1; }

 private:
  T t_;
};
}  // namespace stan
#endif
