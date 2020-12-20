#ifndef STAN_MATH_PRIM_FUN_LDLT_FACTOR_HPP
#define STAN_MATH_PRIM_FUN_LDLT_FACTOR_HPP

#include <stan/math/prim/err/check_square.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <memory>

namespace stan {
namespace math {

/**
 * LDLT_factor is a structure that holds a matrix of type T and the
 * LDLT of its values.
 *
 * If the template argument `alloc_in_arena` is true then the object
 * allocates all necessary memory in the memory arena.
 *
 * LDLT_factors are cheap to copy (they copy references to the
 * underlying storage).
 *
 * @tparam T type of elements in the matrix
 * @tparam alloc_in_arena allocate members of LDLT_factor in the arena
 */
template <typename T, bool alloc_in_arena>
class LDLT_factor;

/**
 * An LDLT_factor with `alloc_in_arena = False` is a structure that
 * holds a matrix of type T and the LDLT of its values.
 *
 * @tparam T type of elements in the matrix
 */
template <typename T>
class LDLT_factor<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, false> {
 private:
  using ldlt_type
      = Eigen::LDLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
  std::shared_ptr<ldlt_type> ldlt_ptr_;

 public:
  template <typename S, require_eigen_t<S>* = nullptr>
  explicit LDLT_factor(const S& matrix) : ldlt_ptr_(new ldlt_type(matrix)) {}

  /**
   * Return a const reference to the underlying matrix
   */
  const auto& matrix() const { return ldlt_ptr_->matrixLDLT(); }

  /**
   * Return a const reference to the LDLT factor of the matrix
   */
  const auto& ldlt() const { return *ldlt_ptr_; }
};

/**
 * Make an LDLT_factor with matrix type `T` and allocate any necessary
 * memory in the arena if any of the types `T` or `Args...` contain
 * `var_value<T>` types.
 *
 * @tparam T Type of matrix to take the LDLT of
 * @tparam Args Other arguments to decide where to allocate LDLT
 * @param A Matrix to take the LDLT of
 * @return LDLT_factor of A
 */
template <typename T, typename... Args, require_matrix_t<T>* = nullptr>
inline auto make_ldlt_factor(const T& A) {
  return LDLT_factor<
      plain_type_t<T>,
      disjunction<is_var<scalar_type_t<return_type_t<T, Args...>>>,
                  is_var<partials_type_t<return_type_t<Args...>>>,
                  is_var<partials_type_t<
                      partials_type_t<return_type_t<Args...>>>>>::value>(A);
}

}  // namespace math
}  // namespace stan

#endif
