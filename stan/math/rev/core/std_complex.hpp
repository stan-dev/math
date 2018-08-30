#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/std_numeric_limits.hpp>
#include <complex>

namespace std {

/**
 * Specialization of complex for var objects.
 *
 * This implementation of std::numeric_limits<stan::math::var>
 * is used to treat var objects like doubles.
 */
/*
template <>
struct complex<stan::math::var> {

};*/

}  // namespace std
#endif
