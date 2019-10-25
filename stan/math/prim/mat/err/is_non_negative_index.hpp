#ifndef STAN_MATH_PRIM_MAT_ERR_IS_NON_NEGATIVE_INDEX_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_NON_NEGATIVE_INDEX_HPP

#include <string>

namespace stan {
namespace math {

inline bool is_validate_non_negative_index(int val) { return val >= 0; }

}  // namespace math
}  // namespace stan
#endif
