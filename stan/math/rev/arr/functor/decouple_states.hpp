#ifndef STAN_MATH_REV_ARR_FUNCTOR_DECOUPLE_STATES_HPP
#define STAN_MATH_REV_ARR_FUNCTOR_DECOUPLE_STATES_HPP

#include <stan/math/rev/core.hpp>

#include <boost/type_traits/is_same.hpp>
#include <vector>

namespace stan {

  namespace math {

    // special case if both (initials and parameters) are known, then
    // this function just returns its input.
    template <>
    inline
    std::vector<std::vector<double> >
    decouple_states(const std::vector<std::vector<double> >& y,
                    const std::vector<double>& y0,
                    const std::vector<double>& theta) {
      return y;
    }

  }  // ns math
}  // ns stan

#endif
