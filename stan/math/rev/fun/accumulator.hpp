#ifndef STAN_MATH_REV_FUN_ACCUMULATOR_HPP
#define STAN_MATH_REV_FUN_ACCUMULATOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/prim/fun/accumulator.hpp>
#include <vector>
#include <type_traits>

namespace stan {
namespace math {

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
class accumulator<T, require_var_t<T>> {
 private:
  static const int max_size_ = 128;
  std::vector<var, arena_allocator<var>> buf_;

  /**
   * Checks if the internal buffer is full and if so reduces it to 1 element.
   */
  void check_size() {
    if (buf_.size() == max_size_) {
      var tmp = stan::math::sum(buf_);
      buf_.resize(1);
      buf_[0] = tmp;
    }
  }

 public:
  /**
   * Constructor. Reserves space.
   */
  accumulator() : buf_() { buf_.reserve(max_size_); }

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
    check_size();
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
    check_size();
    buf_.push_back(stan::math::sum(m));
  }

  /**
   * Add each entry in the specified std vector of values to the buffer.
   *
   * @tparam S type of the matrix
   * @param m Matrix of values to add
   */
  template <typename S, require_std_vector_vt<is_stan_scalar, S>* = nullptr>
  inline void add(const S& xs) {
    check_size();
    this->add(stan::math::sum(xs));
  }

  /**
   * Recursively add each entry in the specified standard vector containint
   * containers to the buffer.
   *
   * @tparam S Type of value to recursively add.
   * @param xs Vector of entries to add
   */
  template <typename S,
            require_std_vector_vt<is_container_or_var_matrix, S>* = nullptr>
  inline void add(const S& xs) {
    check_size();
    for (const auto& i : xs) {
      this->add(i);
    }
  }

#ifdef STAN_OPENCL
  /**
   * Sum each entry and then push to the buffer.
   * @tparam S A Type inheriting from `matrix_cl_base`
   * @param x An OpenCL matrix
   */
  template <typename S,
            require_all_kernel_expressions_and_none_scalar_t<S>* = nullptr>
  inline void add(const var_value<S>& xs) {
    check_size();
    buf_.push_back(stan::math::sum(xs));
  }

  /**
   * Sum each entry and then push to the buffer.
   * @tparam S A Type inheriting from `matrix_cl_base`
   * @param x An OpenCL matrix
   */
  template <typename S,
            require_all_kernel_expressions_and_none_scalar_t<S>* = nullptr>
  inline void add(const S& xs) {
    check_size();
    buf_.push_back(stan::math::sum(xs));
  }

#endif
  /**
   * Return the sum of the accumulated values.
   *
   * @return Sum of accumulated values.
   */
  inline var sum() const { return stan::math::sum(buf_); }
};

}  // namespace math
}  // namespace stan

#endif
