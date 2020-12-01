#ifndef STAN_MATH_PRIM_META_SCALAR_SEQ_VIEW_HPP
#define STAN_MATH_PRIM_META_SCALAR_SEQ_VIEW_HPP

#include <stan/math/prim/meta/plain_type.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/ref_type.hpp>
#include <stan/math/prim/meta/is_autodiff.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_var_matrix.hpp>
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <type_traits>
#include <utility>

namespace stan {
/** \ingroup type_trait
 * scalar_seq_view provides a uniform sequence-like wrapper around either a
 * scalar or a sequence of scalars.
 *
 * @tparam C the container type; will be the scalar type if wrapping a scalar
 * @tparam T the scalar type
 */
template <typename C, typename = void>
class scalar_seq_view;

template <typename C>
class scalar_seq_view<C, require_eigen_vector_t<C>> {
 public:
  template <typename T,
            typename = require_same_t<plain_type_t<T>, plain_type_t<C>>>
  explicit scalar_seq_view(T&& c) : c_(std::forward<T>(c)) {}

  /** \ingroup type_trait
   * Segfaults if out of bounds.
   * @param i index
   * @return the element at the specified position in the container
   */
  auto operator[](int i) const { return c_.coeffRef(i); }
  auto& operator[](int i) { return c_.coeffRef(i); }

  int size() const { return c_.size(); }

  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  auto val(int i) const {
    return c_.coeffRef(i);
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  auto val(int i) const {
    return c_.coeffRef(i).val();
  }

 private:
  ref_type_t<C> c_;
};

template <typename C>
class scalar_seq_view<C, require_var_matrix_t<C>> {
 public:
  template <typename T,
            typename = require_same_t<plain_type_t<T>, plain_type_t<C>>>
  explicit scalar_seq_view(T&& c) : c_(std::forward<T>(c)) {}

  /** \ingroup type_trait
   * Segfaults if out of bounds.
   * @param i index
   * @return the element at the specified position in the container
   */
  auto operator[](int i) const { return c_.coeffRef(i); }
  auto& operator[](int i) { return c_.coeffRef(i); }

  int size() const { return c_.size(); }

  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  auto val(int i) const {
    return c_.val().coeffRef(i);
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  auto val(int i) const {
    return c_.val().coeffRef(i);
  }

 private:
  std::decay_t<C> c_;
};

template <typename C>
class scalar_seq_view<C, require_std_vector_t<C>> {
 public:
  template <typename T,
            typename = require_same_t<plain_type_t<T>, plain_type_t<C>>>
  explicit scalar_seq_view(T&& c) : c_(std::forward<T>(c)) {}

  /** \ingroup type_trait
   * Segfaults if out of bounds.
   * @param i index
   * @return the element at the specified position in the container
   */
  decltype(auto) operator[](int i) const { return c_[i]; }
  decltype(auto) operator[](int i) { return c_[i]; }
  int size() const { return c_.size(); }

  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  auto val(int i) const {
    return c_[i];
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  auto val(int i) const {
    return c_[i].val();
  }

 private:
  ref_type_t<C> c_;
};

/** \ingroup type_trait
 * This specialization handles wrapping a scalar as if it were a sequence.
 *
 * @tparam T the scalar type
 */
template <typename C>
class scalar_seq_view<C, require_stan_scalar_t<C>> {
 public:
  explicit scalar_seq_view(const C& t) : t_(t) {}

  auto& operator[](int /* i */) const { return t_; }
  auto& operator[](int /* i */) { return t_; }
  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  auto val(int /* i */) const {
    return t_;
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  auto val(int /* i */) const {
    return t_.val();
  }

  int size() const { return 1; }

 private:
  std::decay_t<C> t_;
};
}  // namespace stan
#endif
