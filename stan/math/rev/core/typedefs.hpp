#ifndef STAN_MATH_REV_CORE_TYPEDEFS_HPP
#define STAN_MATH_REV_CORE_TYPEDEFS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/Eigen_NumTraits.hpp>

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

/**
 * Types used for profiling.
 */
struct profile_meta {
  bool fwd_pass_active;
  bool rev_pass_active;
  std::chrono::time_point<std::chrono::steady_clock> fwd_pass_start;
  std::chrono::time_point<std::chrono::steady_clock> rev_pass_start;
};

struct profile_info {
  double fwd_pass;
  double rev_pass;
  int count_rev;
  profile_meta meta;
};

using profile_map = std::map<std::string, profile_info>;

}  // namespace math
}  // namespace stan
#endif
