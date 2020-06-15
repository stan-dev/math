#ifndef STAN_MATH_PRIM_META_vector_SEQ_VIEW_HPP
#define STAN_MATH_PRIM_META_vector_SEQ_VIEW_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {

/** \ingroup type_trait
 * This class provides a low-cost wrapper for situations where you either need
 * an Eigen Vector or RowVector or a std::vector of them and you want to be
 * agnostic between those two options. This is similar to scalar_seq_view but
 * instead of being a sequence-like view over a scalar or seq of scalars, it's
 * a sequence-like view over a Vector or seq of Vectors. Notably this version
 * only allows std::vectors as the container type, since we would have
 * difficulty figuring out which contained type was the container otherwise.
 *
 * @tparam T the wrapped type, either a Vector or std::vector of them.
 */
template <typename T, typename = void>
class vector_seq_view;

/** \ingroup type_trait
 * This class provides a low-cost wrapper for situations where you either need
 * an Eigen Vector or RowVector or a std::vector of them and you want to be
 * agnostic between those two options. This is similar to scalar_seq_view but
 * instead of being a sequence-like view over a scalar or seq of scalars, it's
 * a sequence-like view over a Vector or seq of Vectors. Notably this version
 * only allows std::vectors as the container type, since we would have
 * difficulty figuring out which contained type was the container otherwise.
 *
 * @tparam T the underlying Vector or expression type
 */
template <typename T>
class vector_seq_view<T, require_eigen_t<T>> {
 public:
  explicit vector_seq_view(T& m) : m_(m) {}
  int size() const { return 1; }
  T& operator[](int /* i */) const { return m_; }

 private:
  T& m_;
};

/** \ingroup type_trait
 * This class provides a low-cost wrapper for situations where you either need
 * an Eigen Vector or RowVector or a std::vector of them and you want to be
 * agnostic between those two options. This is similar to scalar_seq_view but
 * instead of being a sequence-like view over a scalar or seq of scalars, it's
 * a sequence-like view over a Vector or seq of Vectors. Notably this version
 * only allows std::vectors as the container type, since we would have
 * difficulty figuring out which contained type was the container otherwise.
 *
 * @tparam T std vector of eigen types
 */
template <typename T>
class vector_seq_view<T, require_std_vector_vt<is_eigen, T>> {
  using inner = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                                   const value_type_t<T>, value_type_t<T>>;

 public:
  explicit vector_seq_view(T& v) : v_(v) {}
  int size() const { return v_.size(); }
  inner& operator[](int i) const { return v_[i]; }

 private:
  T& v_;
};

}  // namespace stan

#endif
