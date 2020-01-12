#ifndef STAN_MATH_REV_FUN_SUM_HPP
#define STAN_MATH_REV_FUN_SUM_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/sum.hpp>
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
 protected:
  template <typename Derived>
  inline static double sum_of_val(const Eigen::DenseBase<Derived>& v) {
    return Eigen::Ref<const matrix_v>(v).val().sum();
  }

 public:
  template <int R1, int C1>
  explicit sum_eigen_v_vari(const Eigen::Matrix<var, R1, C1>& v1)
      : sum_v_vari(
            sum_of_val(v1),
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
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param m Specified matrix or vector.
 * @return Sum of coefficients of matrix.
 */
template <int R, int C>
inline var sum(const Eigen::Matrix<var, R, C>& m) {
  if (m.size() == 0) {
    return 0.0;
  }
  return var(new sum_eigen_v_vari(m));
}

}  // namespace math
}  // namespace stan
#endif
