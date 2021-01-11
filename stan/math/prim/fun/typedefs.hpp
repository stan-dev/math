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

/**
 * Types used for profiling.
 */
struct profile_meta {
  bool fwd_pass_active;
  bool rev_pass_active;
  std::chrono::time_point<std::chrono::steady_clock> fwd_pass_start_tp;
  std::chrono::time_point<std::chrono::steady_clock> rev_pass_start_tp;
  size_t start_chain_stack_size;
  size_t start_nochain_stack_size;
};

struct profile_info {
  double fwd_pass_time;
  double rev_pass_time;
  size_t chain_stack_size_sum;
  size_t nochain_stack_size_sum;
  size_t n_fwd_pass;
  size_t n_rev_pass;
  profile_meta meta;
};

using profile_key = std::pair<const char*, std::thread::id>;

using profile_map = std::map<profile_key, profile_info>;

}  // namespace math
}  // namespace stan

#endif
