#ifndef STAN_MATH_REV_META_VAR_MATRIX_CONVERTER_HPP
#define STAN_MATH_REV_META_VAR_MATRIX_CONVERTER_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>

namespace stan {

template <typename T>
using var_matrix_converter_t = math::var_value<math::promote_scalar_t<double, T>>;

}

#endif
