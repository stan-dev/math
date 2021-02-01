#ifndef STAN_MATH_PRIM_ERR_VALIDATE_UNIT_VECTOR_INDEX_HPP
#define STAN_MATH_PRIM_ERR_VALIDATE_UNIT_VECTOR_INDEX_HPP

#include <stan/math/prim/meta.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

namespace stan {
namespace math {

/**
 * Check that unit vector is at least size 2
 *
 * @param var_name Name of unit vector variable
 * @param expr Expression in which variable is declared
 * @param val Size of unit vector
 * @throw std::invalid_argument if simplex size is less than 2
 */
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      if (val == 1) {
        msg << "Found dimension size one in unit vector declaration."
            << " One-dimensional unit vector is discrete"
            << " but the target distribution must be continuous."
            << " variable=" << var_name
            << "; dimension size expression=" << expr;
      } else {
        msg << "Found dimension size less than one in unit vector declaration"
            << "; variable=" << var_name
            << "; dimension size expression=" << expr
            << "; expression value=" << val;
      }
      std::string msg_str(msg.str());
      throw std::invalid_argument(msg_str.c_str());
    }();
  }
}

}  // namespace math
}  // namespace stan
#endif
