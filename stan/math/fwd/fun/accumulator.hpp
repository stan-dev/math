#ifndef STAN_MATH_FWD_FUN_ACCUMULATOR_HPP
#define STAN_MATH_FWD_FUN_ACCUMULATOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/sum.hpp>
#include <stan/math/prim/fun/accumulator.hpp>

#include <vector>
#include <type_traits>

namespace stan {
namespace math {
template <typename T, typename>
class accumulator;
/**
 * Class to accumulate values and eventually return their sum.  If
 * no values are ever added, the return value is 0.
 *
 * This class is useful for speeding up autodiff of long sums
 * because it uses the <code>sum()</code> operation (either from
 * <code>stan::math</code> or one defined by argument-dependent lookup.
 *
 * @tparam T Type of scalar added
 */
template <typename T>
class accumulator<T, require_fvar_t<T>> {
 private:
  std::vector<T> buf_;

 public:
  /**
   * Add the specified arithmetic type value to the buffer after
   * static casting it to the class type <code>T</code>.
   *
   * <p>See the std library doc for <code>std::is_arithmetic</code>
   * for information on what counts as an arithmetic type.
   *
   * @tparam S Type of argument
   * @param x Value to add
   */
  template <typename S, typename = require_stan_scalar_t<S>>
  inline void add(S x) {
    buf_.push_back(x);
  }

  /**
   * Add each entry in the specified matrix, vector, or row vector
   * of values to the buffer.
   *
   * @tparam S type of the matrix
   * @param m Matrix of values to add
   */
  template <typename S, require_matrix_t<S>* = nullptr>
  inline void add(const S& m) {
    buf_.push_back(stan::math::sum(m));
  }

  /**
   * Recursively add each entry in the specified standard vector
   * to the buffer.  This will allow vectors of primitives,
   * autodiff variables to be added; if the vector entries
   * are collections, their elements are recursively added.
   *
   * @tparam S Type of value to recursively add.
   * @param xs Vector of entries to add
   */
  template <typename S>
  inline void add(const std::vector<S>& xs) {
    for (size_t i = 0; i < xs.size(); ++i) {
      this->add(xs[i]);
    }
  }

#ifdef STAN_OPENCL

  /**
   * Sum each entry and then push to the buffer.
   * @tparam S A Type inheriting from `matrix_cl_base`
   * @param xs An OpenCL matrix
   */
  template <typename S,
            require_all_kernel_expressions_and_none_scalar_t<S>* = nullptr>
  inline void add(const S& xs) {
    buf_.push_back(stan::math::sum(xs));
  }

#endif

  /**
   * Return the sum of the accumulated values.
   *
   * @return Sum of accumulated values.
   */
  inline T sum() const { return stan::math::sum(buf_); }
};

}  // namespace math
}  // namespace stan

#endif
