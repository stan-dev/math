#ifndef STAN_MATH_REV_FUN_SUM_HPP
#define STAN_MATH_REV_FUN_SUM_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Class for sums of variables constructed with standard vectors.
 * There's an extension for Eigen matrices.
 */
class sum_v_vari : public vari {
 protected:
  vari** v_;
  size_t length_;

  inline static double sum_of_val(const std::vector<var>& v) {
    double result = 0;
    for (auto x : v) {
      result += x.val();
    }
    return result;
  }

 public:
  explicit sum_v_vari(double value, vari** v, size_t length)
      : vari(value), v_(v), length_(length) {}

  explicit sum_v_vari(const std::vector<var>& v1)
      : vari(sum_of_val(v1)),
        v_(reinterpret_cast<vari**>(ChainableStack::instance_->memalloc_.alloc(
            v1.size() * sizeof(vari*)))),
        length_(v1.size()) {
    for (size_t i = 0; i < length_; i++) {
      v_[i] = v1[i].vi_;
    }
  }

  virtual void chain() {
    for (size_t i = 0; i < length_; i++) {
      v_[i]->adj_ += adj_;
    }
  }
};

/**
 * Returns the sum of the entries of the specified vector.
 *
 * @param m Vector.
 * @return Sum of vector entries.
 */
inline var sum(const std::vector<var>& m) {
  if (m.empty()) {
    return 0.0;
  }
  return var(new sum_v_vari(m));
}

/**
 * Class for representing sums with constructors for Eigen.
 * The <code>chain()</code> method and member variables are
 * managed by the superclass <code>sum_v_vari</code>.
 */
class sum_eigen_v_vari : public sum_v_vari {
 public:
  template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
  explicit sum_eigen_v_vari(const EigMat& v1)
      : sum_v_vari(
            v1.val().sum(),
            reinterpret_cast<vari**>(ChainableStack::instance_->memalloc_.alloc(
                v1.size() * sizeof(vari*))),
            v1.size()) {
    Eigen::Map<matrix_vi>(v_, v1.rows(), v1.cols()) = v1.vi();
  }
};

/**
 * Returns the sum of the coefficients of the specified
 * matrix, column vector or row vector.
 *
 * @tparam T type of the matrix of vector (Must be derived from \c
 * Eigen::MatrixBase and contain \c var scalars)
 * @param m Specified matrix or vector.
 * @return Sum of coefficients of matrix.
 */
template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
inline var sum(const EigMat& m) {
  if (m.size() == 0) {
    return 0.0;
  }
  const Eigen::Ref<const plain_type_t<EigMat>>& m_ref = m;
  return var(new sum_eigen_v_vari(m_ref));
}

}  // namespace math
}  // namespace stan
#endif
