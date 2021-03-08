#ifndef STAN_MATH_REV_FUN_SUM_HPP
#define STAN_MATH_REV_FUN_SUM_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/typedefs.hpp>
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
 * Returns the sum of the coefficients of the specified
 * matrix.
 *
 * @tparam T type of the matrix of vector. Can be either a var matrix or
 *  matrix of vars.
 * @param x Specified var_value containing a matrix or vector.
 * @return Sum of coefficients of matrix.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline var sum(const T& x) {
  arena_t<T> x_arena = x;
  return make_callback_var(sum(x_arena.val()), [x_arena](auto& vi) mutable {
    x_arena.adj().array() += vi.adj();
  });
}

}  // namespace math
}  // namespace stan
#endif
