#ifndef STAN_MATH_PRIM_ERR_CHECK_VECTOR_HPP
#define STAN_MATH_PRIM_ERR_CHECK_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/invalid_argument.hpp>
#include <sstream>
#include <string>
#include <typeinfo>

namespace stan {
namespace math {

/**
 * Check if the matrix is either a row vector or column vector.
 * This function checks the runtime size of the matrix to check
 * whether it is a row or column vector.
 * @tparam EigMat A type derived from `EigenBase` with dynamic rows and columns
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param x Matrix
 * @throw <code>std::invalid_argument</code> if x is not a row or column
 *   vector.
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline void check_vector(const char* function, const char* name,
                         const EigMat& x) {
  if (x.rows() == 1 || x.cols() == 1) {
    return;
  } else {
    std::ostringstream msg;
    msg << ") has " << x.rows() << " rows and " << x.cols()
        << " columns but it should be a vector so it should "
        << "either have 1 row or 1 column";
    std::string msg_str(msg.str());
    invalid_argument(function, name, value_type_t<EigMat>(), "(",
                     msg_str.c_str());
  }
}

/**
 * Overload for check_vector that returns immedietly at compile time for
 * Eigen types with compile time rows or columns equal to 1.
 * whether it is a row or column vector.
 * @tparam EigVec A type derived from `EigenBase` with 1 row or column
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param x Matrix
 * @throw <code>std::invalid_argument</code> if x is not a row or column
 *   vector.
 */
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
constexpr inline void check_vector(const char* function, const char* name,
                                   const EigVec& x) {
  return;
}
}  // namespace math
}  // namespace stan
#endif
