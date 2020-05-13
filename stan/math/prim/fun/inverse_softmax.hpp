#ifndef STAN_MATH_ARR_SCAL_FUN_INVERSE_SOFTMAX_HPP
#define STAN_MATH_ARR_SCAL_FUN_INVERSE_SOFTMAX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Writes the inverse softmax of the simplex argument into the second
 * argument.  See <code>softmax</code> for the inverse
 * function and a definition of the relation.
 *
 * The inverse softmax function is defined by
 *
 * \f$\mbox{inverse\_softmax}(x)[i] = \log x[i]\f$.
 *
 * This function defines the inverse of <code>softmax</code>
 * up to a scaling factor.
 *
 * Because of the definition, values of 0.0 in the simplex
 * are converted to negative infinity, and values of 1.0
 * are converted to 0.0.
 *
 * There is no check that the input vector is a valid simplex vector.
 *
 * @tparam Vector type of the simplex vector
 * @param simplex Simplex vector input.
 * @param y Vector into which the inverse softmax is written.
 * @throw std::invalid_argument if size of the input and
 *    output vectors differ.
 */
template <typename Vector, require_vector_t<Vector>* = nullptr>
void inverse_softmax(const Vector& simplex, Vector& y) {
  check_matching_sizes("inverse_softmax", "simplex", simplex, "y", y);
  y = log(simplex);
}

}  // namespace math
}  // namespace stan
#endif
