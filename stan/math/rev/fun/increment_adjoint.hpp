#ifndef STAN_MATH_REV_FUN_INCREMENT_ADJOINT_HPP
#define STAN_MATH_REV_FUN_INCREMENT_ADJOINT_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/sum.hpp>

namespace stan {
namespace math {

/**
 * Incrementing adjoint of a prim type is a no-op.
 * @tparam T1 a non-var type
 * @tparam T2 type of increment
 */
template<typename T1, typename T2, require_not_st_var<T1>* = nullptr>
void increment_adjoint(T1&, const T2&){}

/**
 * Increments adjoint of given scalar var by given value.
 * @tparam T2 type of increment
 * @param scal scalar var, adjoint of which is to be incremented
 * @param increment increment to add to the adjoint
 */
template<typename T2>
void increment_adjoint(var& scal, const T2& increment){
  scal.adj() += sum(increment);
}

/**
 * Increments adjoints of given rev matrix by given values.
 * @tparam T1 a rev matrix type
 * @tparam T2 type of increment
 * @param mat matrix, adjoint of which are to be incremented
 * @param increment increment to add to the adjoints
 */
template<typename T1, typename T2, require_rev_matrix_t<T1>* = nullptr,
         require_eigen_t<T2>* = nullptr>
void increment_adjoint(T1& mat, const T2& increment){
  mat.adj() += increment;
}

/**
 * Increments adjoints of given rev matrix by given scalar value.
 * @tparam T1 a rev matrix type
 * @tparam T2 type of increment
 * @param mat matrix, adjoint of which are to be incremented
 * @param increment increment to add to the adjoints
 */
template<typename T1, typename T2, require_rev_matrix_t<T1>* = nullptr,
         require_arithmetic_t<T2>* = nullptr>
void increment_adjoint(T1& mat, const T2& increment){
  mat.adj().array() += increment;
}

}
}

#endif // INCREMENT_ADJOINT_HPP
