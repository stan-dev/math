#ifndef STAN_MATH_REV_FUN_DOT_SELF_HPP
#define STAN_MATH_REV_FUN_DOT_SELF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <vector>

namespace stan {
namespace math {

namespace internal {
class dot_self_vari : public vari {
 protected:
  vari** v_;
  size_t size_;

 public:
  dot_self_vari(vari** v, size_t size)
      : vari(Eigen::Map<vector_vi>(v, size).val().squaredNorm()),
        v_(v),
        size_(size) {}
  template <typename T, require_eigen_t<T>* = nullptr>
  explicit dot_self_vari(const T& v)
      : vari(v.val().squaredNorm()), size_(v.size()) {
    v_ = reinterpret_cast<vari**>(
        ChainableStack::instance_->memalloc_.alloc(size_ * sizeof(vari*)));
    Eigen::Map<matrix_vi>(v_, v.rows(), v.cols()) = v.vi();
  }
  virtual void chain() {
    Eigen::Map<vector_vi> v_map(v_, size_);
    v_map.adj() += adj_ * 2.0 * v_map.val();
  }
};
}  // namespace internal

/**
 * Returns the dot product of a vector of var with itself.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile time dimension equal to 1)
 * @param[in] v Vector.
 * @return Dot product of the vector with itself.
 */
template <typename T, require_eigen_vector_vt<is_var, T>* = nullptr>
inline var dot_self(const T& v) {
  const Eigen::Ref<const plain_type_t<T>>& v_ref = v;
  return {new internal::dot_self_vari(v_ref)};
}

}  // namespace math
}  // namespace stan
#endif
