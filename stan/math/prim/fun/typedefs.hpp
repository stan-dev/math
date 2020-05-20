#ifndef STAN_MATH_PRIM_FUN_TYPEDEFS_HPP
#define STAN_MATH_PRIM_FUN_TYPEDEFS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Type for sizes and indexes in an Eigen matrix with double elements.
 */
using size_type
    = index_type_t<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>;

/**
 * Type for matrix of double values.
 */
using matrix_d = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * Type for (column) vector of double values.
 */
using vector_d = Eigen::Matrix<double, Eigen::Dynamic, 1>;

/**
 * Type for (row) vector of double values.
 */
using row_vector_d = Eigen::Matrix<double, 1, Eigen::Dynamic>;

typedef struct {
    std::chrono::time_point<std::chrono::steady_clock> fwd_pass_time_start;
    std::chrono::time_point<std::chrono::steady_clock> fwd_pass_time_stop;
    std::chrono::time_point<std::chrono::steady_clock> bkcwd_pass_time_start;
    std::chrono::time_point<std::chrono::steady_clock> bkcwd_pass_time_stop;
} profile;

using profiles = std::map<int, profile>;

}  // namespace math
}  // namespace stan

#endif
