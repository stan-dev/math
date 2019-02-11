#ifndef STAN_MATH_PRIM_MAT_ERR_IS_SQUARE_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_SQUARE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/is_scal_size_match.hpp>

namespace stan {
namespace math {

/**
 * Check if the specified matrix is square.
 *
 * This check allows 0x0 matrices.
 *
 * @tparam T type of scalar
 *
 * @param y Matrix to test
 * 
 * @return <code>true</code> if matrix is square
 */
template<typename T_y>
inline bool is_square(
	    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  is_scal_size_match(y.rows(), y.cols());
}

} // namespace math
} // namepsace stan
#endif
