#ifndef STAN_MATH_PRIM_SCAL_FUN_DISTANCE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_DISTANCE_HPP

#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/abs.hpp>

namespace stan {
  namespace math {

    /**
     * Returns the squared distance.
     *
     * @param x1 First vector.
     * @param x2 Second vector.
     * @return Distance between two scalars
     * @throw std::domain_error If the vectors are not the same
     * size or if they are both not vector dimensioned.
     */
    template<typename T1, typename T2>
    inline typename boost::math::tools::promote_args<T1, T2>::type
    distance(const T1& x1,
                     const T2& x2) {
      check_finite("distance", "x1", x1);
      check_finite("distance", "x2", x2);
      return abs(x1 - x2);
    }
  }
}
#endif
