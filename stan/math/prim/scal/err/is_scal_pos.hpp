#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_SCAL_POS
#define STAN_MATH_PRIM_SCAL_ERR_IS_SCAL_POS

#include <stan/math/prim/scal/meta/value_type.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/get.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <boost/type_traits/is_unsigned.hpp>

namespace stan {
namespace math {

template <typename T_y>
struct positive<T_y, true> {
  static bool check(const T_y& y) {
    for (size_t n = 0; n < length(y); n++) {
      if (!boost::is_unsigned<typename value_type<T_y>::type>::value
          && !(stan::get(y, n) > 0))
        return false;
    }
    return true;
  }
};

/**
 * Check if <code>y</code> is positive.
 *
 * This function is vectorized and will check each element of
 * <code>y</code>.
 *
 * @tparam T_y Type of y
 *
 * @param y Variable to check
 *
 * @return <code>true</code> is vector contains only positive elements
 */
template <typename T_y>
inline bool is_scal_pos(const T_y& y) {
  positive<T_y>::check(y);
}

}  // namespace math
}  // namespace stan
#endif
