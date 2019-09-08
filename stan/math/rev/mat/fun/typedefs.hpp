#ifndef STAN_MATH_REV_MAT_FUN_TYPEDEFS_HPP
#define STAN_MATH_REV_MAT_FUN_TYPEDEFS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

using size_type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Index;

/**
 * The type of a matrix holding <code>var</code>
 * values.
 */
using matrix_v = Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * The type of a (column) vector holding <code>var</code>
 * values.
 */
using vector_v = Eigen::Matrix<var, Eigen::Dynamic, 1>;

/**
 * The type of a row vector holding <code>var</code>
 * values.
 */
using row_vector_v = Eigen::Matrix<var, 1, Eigen::Dynamic>;

/**
 * The type of a matrix holding <code>vari*</code>
 * values.
 */
using matrix_vi = Eigen::Matrix<vari*, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * The type of a (column) vector holding <code>vari*</code>
 * values.
 */
using vector_vi = Eigen::Matrix<vari*, Eigen::Dynamic, 1>;

/**
 * The type of a row vector holding <code>vari*</code>
 * values.
 */
using row_vector_vi = Eigen::Matrix<vari*, 1, Eigen::Dynamic>;

}  // namespace math
}  // namespace stan
#endif
