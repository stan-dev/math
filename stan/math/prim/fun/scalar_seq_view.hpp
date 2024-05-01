#ifndef STAN_MATH_PRIM_FUN_SCALAR_SEQ_VIEW_HPP
#define STAN_MATH_PRIM_FUN_SCALAR_SEQ_VIEW_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/meta/is_tuple.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/cumulative_sum.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/apply_at.hpp>
#include <stan/math/prim/functor/for_each.hpp>

namespace stan {
namespace internal {

/**
 * Utility function to help with indexing nested containers. Given an index
 * and a vector containing the cumulative sum of the container sizes, returns
 * the index of the desired container and the index of the desired value
 * within that container.
 *
 * @param idx The index passed to operator[]
 * @param cumulative_sizes The cumulative sum of container sizes
 * @return Two-element array containing index of container and index of value
 *            within container
 */
std::array<size_t, 2> lookup_index(
    size_t idx, const std::vector<size_t>& cumulative_sizes) {
  size_t element
      = std::find_if(cumulative_sizes.cbegin(), cumulative_sizes.cend(),
                     [idx](size_t i) { return (idx + 1) <= i; })
        - cumulative_sizes.cbegin();
  size_t offset = element == 0 ? 0 : cumulative_sizes[element - 1];
  return {element, idx - offset};
}
}  // namespace internal

/**
 * scalar_seq_view provides a uniform sequence-like wrapper around either a
 * scalar or a sequence of scalars.
 *
 * @tparam C the container type; will be the scalar type if wrapping a scalar
 * @tparam T the scalar type
 */
template <typename C, typename = void>
class scalar_seq_view;

template <typename C>
class scalar_seq_view<C, require_eigen_t<C>> {
 public:
  template <typename T,
            typename = require_same_t<plain_type_t<T>, plain_type_t<C>>>
  explicit scalar_seq_view(T&& c) : c_(std::forward<T>(c)) {}

  /**
   * Segfaults if out of bounds.
   * @param i index
   * @return the element at the specified position in the container
   */
  inline auto operator[](size_t i) const { return c_.coeff(i); }
  inline auto& operator[](size_t i) { return c_.coeffRef(i); }

  inline auto size() const noexcept { return c_.size(); }

  inline const value_type_t<C>* data() const noexcept { return c_.data(); }
  inline value_type_t<C>* data() noexcept { return c_.data(); }

  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return c_.coeffRef(i);
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return c_.coeffRef(i).val();
  }

 private:
  plain_type_t<C> c_;
};

template <typename C>
class scalar_seq_view<C, require_var_matrix_t<C>> {
 public:
  template <typename T,
            typename = require_same_t<plain_type_t<T>, plain_type_t<C>>>
  explicit scalar_seq_view(T&& c) : c_(std::forward<T>(c)) {}

  /**
   * Segfaults if out of bounds.
   * @param i index
   * @return the element at the specified position in the container
   */
  inline auto operator[](size_t i) const { return c_.coeff(i); }

  inline const auto* data() const noexcept { return c_.vi_; }
  inline auto* data() noexcept { return c_.vi_; }

  inline auto size() const noexcept { return c_.size(); }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  inline auto val(size_t i) const {
    return c_.val().coeff(i);
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  inline auto& val(size_t i) {
    return c_.val().coeffRef(i);
  }

 private:
  std::decay_t<C> c_;
};

template <typename C>
class scalar_seq_view<C, require_std_vector_vt<is_stan_scalar, C>> {
 public:
  template <typename T,
            typename = require_same_t<plain_type_t<T>, plain_type_t<C>>>
  explicit scalar_seq_view(T&& c) : c_(std::forward<T>(c)) {}

  /**
   * Segfaults if out of bounds.
   * @param i index
   * @return the element at the specified position in the container
   */
  inline auto operator[](size_t i) const { return c_[i]; }
  inline auto& operator[](size_t i) { return c_[i]; }
  inline auto size() const noexcept { return c_.size(); }
  inline const auto* data() const noexcept { return c_.data(); }

  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return c_[i];
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return c_[i].val();
  }

 private:
  std::decay_t<C> c_;
};

template <typename C>
class scalar_seq_view<C, require_std_vector_vt<is_container, C>> {
 public:
  template <typename T,
            typename = require_same_t<plain_type_t<T>, plain_type_t<C>>>
  explicit scalar_seq_view(T&& c) {
    std::vector<size_t> sizes_;
    sizes_.reserve(c.size());
    c_.reserve(c.size());

    for (auto&& c_val : c) {
      c_.push_back(scalar_seq_view<value_type_t<C>>(c_val));
      sizes_.push_back(c_.back().size());
    }
    cumulative_sizes_ = math::cumulative_sum(sizes_);
  }

  inline auto size() const noexcept { return cumulative_sizes_.back(); }

  inline auto operator[](size_t i) const {
    std::array<size_t, 2> idxs = internal::lookup_index(i, cumulative_sizes_);
    return c_[idxs[0]][idxs[1]];
  }

  inline auto& operator[](size_t i) {
    std::array<size_t, 2> idxs = internal::lookup_index(i, cumulative_sizes_);
    return c_[idxs[0]][idxs[1]];
  }

  inline const auto* data() const noexcept { return c_.data(); }

  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return this[i];
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return this[i].val();
  }

 private:
  math::modify_nested_value_type_t<scalar_seq_view, C> c_;
  std::vector<size_t> cumulative_sizes_;
};

template <typename C>
class scalar_seq_view<C, math::require_tuple_t<C>> {
 public:
  template <typename T>
  explicit scalar_seq_view(T&& c)
      : c_(math::apply(
          [](auto&&... args) {
            return std::make_tuple(
                scalar_seq_view<std::decay_t<decltype(args)>>(args)...);
          },
          std::forward<decltype(c)>(c))) {
    std::vector<size_t> sizes_;
    math::for_each([&sizes_](auto&& elem) { sizes_.push_back(elem.size()); },
                   std::forward<decltype(c_)>(c_));
    cumulative_sizes_ = math::cumulative_sum(sizes_);
  }

  inline const auto operator[](size_t i) const {
    std::array<size_t, 2> idxs = internal::lookup_index(i, cumulative_sizes_);
    auto index_func = [idxs](auto&& tuple_elem) { return tuple_elem[idxs[1]]; };
    return math::apply_at(index_func, idxs[0], std::forward<decltype(c_)>(c_));
  }

  inline auto& operator[](size_t i) {
    std::array<size_t, 2> idxs = internal::lookup_index(i, cumulative_sizes_);
    auto index_func = [idxs](auto&& tuple_elem) -> decltype(auto) {
      return tuple_elem[idxs[1]];
    };
    return math::apply_at(index_func, idxs[0], std::forward<decltype(c_)>(c_));
  }

  inline auto size() const noexcept { return cumulative_sizes_.back(); }

  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return this[i];
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return this[i].val();
  }

 private:
  math::modify_nested_value_type_t<scalar_seq_view, C> c_;
  std::vector<size_t> cumulative_sizes_;
};

template <typename C>
class scalar_seq_view<C, require_t<std::is_pointer<C>>> {
 public:
  template <typename T,
            typename = require_same_t<plain_type_t<T>, plain_type_t<C>>>
  explicit scalar_seq_view(const T& c) : c_(c) {}

  /**
   * Segfaults if out of bounds.
   * @param i index
   * @return the element at the specified position in the container
   */
  inline auto operator[](size_t i) const { return c_[i]; }
  inline auto size() const noexcept {
    static_assert(1, "Cannot Return Size of scalar_seq_view with pointer type");
  }
  inline const auto* data() const noexcept { return &c_[0]; }

  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return c_[i];
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return c_[i].val();
  }

 private:
  const C& c_;
};

/**
 * This specialization handles wrapping a scalar as if it were a sequence.
 *
 * @tparam T the scalar type
 */
template <typename C>
class scalar_seq_view<C, require_stan_scalar_t<C>> {
 public:
  explicit scalar_seq_view(const C& t) noexcept : t_(t) {}

  inline auto operator[](size_t /* i */) const { return t_; }
  inline auto& operator[](size_t /* i */) { return t_; }

  template <typename T = C, require_st_arithmetic<T>* = nullptr>
  inline decltype(auto) val(size_t /* i */) const noexcept {
    return t_;
  }

  template <typename T = C, require_st_autodiff<T>* = nullptr>
  inline decltype(auto) val(size_t /* i */) const noexcept {
    return t_.val();
  }

  static constexpr auto size() { return 1; }
  inline const auto* data() const noexcept { return &t_; }
  inline auto* data() noexcept { return &t_; }

 private:
  std::decay_t<C> t_;
};
}  // namespace stan
#endif
