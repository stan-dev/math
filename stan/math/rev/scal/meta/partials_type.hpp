#ifndef STAN_MATH_REV_SCAL_META_PARTIALS_TYPE_HPP
#define STAN_MATH_REV_SCAL_META_PARTIALS_TYPE_HPP

#include <stan/math/prim/scal/meta/partials_type.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {

template <>
struct partials_type<stan::math::var> {
  typedef double type;
};

}  // namespace stan
#endif
