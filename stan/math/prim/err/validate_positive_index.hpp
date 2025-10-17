#ifndef STAN_MATH_PRIM_ERR_VALIDATE_POSITIVE_INDEX_HPP
#define STAN_MATH_PRIM_ERR_VALIDATE_POSITIVE_INDEX_HPP

#include <stan/math/prim/meta.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

namespace stan {
namespace math {

/**
 * Check that size is at least 1. Used for simplexes and
 * other constraints that do a (size - 1) operation.
 *
 * @param var_name Name of variable
 * @param expr Expression in which variable is declared
 * @param val Size to check
 * @throw std::invalid_argument if size is less than 1
 */
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << "Found dimension size less than one in constrained type "
             "declaration (simplex, sum_to_zero_vector, etc.)"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
      std::string msg_str(msg.str());
      throw std::invalid_argument(msg_str.c_str());
    }();
  }
}

}  // namespace math
}  // namespace stan
#endif
