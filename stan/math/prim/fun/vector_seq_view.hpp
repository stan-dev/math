#ifndef STAN_MATH_PRIM_FUN_VECTOR_SEQ_VIEW_HPP
#define STAN_MATH_PRIM_FUN_VECTOR_SEQ_VIEW_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {

/**
 * This class provides a low-cost wrapper for situations where you either need
 * an Eigen Vector or RowVector or a std::vector of them and you want to be
 * agnostic between those two options. This is similar to scalar_seq_view but
 * instead of being a sequence-like view over a scalar or seq of scalars, it's
 * a sequence-like view over a Vector or seq of Vectors. Notably this version
 * only allows std::vectors as the outer container type, since we would have
 * difficulty figuring out which contained type was the container otherwise.
 *
 * @tparam T the wrapped type, either a Vector or std::vector of them.
 */
template <typename T, typename = void>
class vector_seq_view;

/**
 * This class provides a low-cost wrapper for situations where you either need
 * an Eigen Vector or RowVector or a std::vector of them and you want to be
 * agnostic between those two options. This is similar to scalar_seq_view but
 * instead of being a sequence-like view over a scalar or seq of scalars, it's
 * a sequence-like view over a Vector or seq of Vectors. Notably this version
 * only allows std::vectors as the outer container type, since we would have
 * difficulty figuring out which contained type was the container otherwise.
 *
 * @tparam T the type of the underlying Vector
 */
template <typename T>
class vector_seq_view<T, require_matrix_t<T>> {
 public:
  explicit vector_seq_view(const T& m) : m_(m) {}
  static constexpr auto size() { return 1; }
  inline const auto& operator[](size_t /* i */) const noexcept { return m_; }

  template <typename C = T, require_st_arithmetic<C>* = nullptr>
  inline const auto& val(size_t /* i */) const noexcept {
    return m_;
  }

  template <typename C = T, require_st_autodiff<C>* = nullptr>
  inline auto val(size_t /* i */) const noexcept {
    return m_.val();
  }

 private:
  const ref_type_t<T> m_;
};

namespace internal {
template <typename T>
using is_matrix_or_std_vector
    = math::disjunction<is_matrix<T>, is_std_vector<T>>;
}

/**
 * This class provides a low-cost wrapper for situations where you either need
 * an Eigen Vector or RowVector or a std::vector of them and you want to be
 * agnostic between those two options. This is similar to scalar_seq_view but
 * instead of being a sequence-like view over a scalar or seq of scalars, it's
 * a sequence-like view over a Vector or seq of Vectors. Notably this version
 * only allows std::vectors as the outer container type, since we would have
 * difficulty figuring out which contained type was the container otherwise.
 *
 * @tparam S the type inside of the std::vector
 */
template <typename T>
class vector_seq_view<
    T, require_std_vector_vt<internal::is_matrix_or_std_vector, T>> {
 public:
  explicit vector_seq_view(const T& v) noexcept : v_(v) {}
  inline auto size() const noexcept { return v_.size(); }

  inline decltype(auto) operator[](size_t i) const { return v_[i]; }

  template <typename C = T, require_st_arithmetic<C>* = nullptr>
  inline decltype(auto) val(size_t i) const {
    return v_[i];
  }

  template <typename C = T, require_st_autodiff<C>* = nullptr>
  inline auto val(size_t i) const {
    return value_of(v_[i]);
  }

 private:
  const T& v_;
};

}  // namespace stan

#endif
