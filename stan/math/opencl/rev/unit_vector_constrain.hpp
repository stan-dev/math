#ifndef STAN_MATH_OPENCL_REV_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_REV_UNIT_VECTOR_CONSTRAIN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return the unit length vector corresponding to the free vector y.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 * @tparam T type of the value of input matrix
 * @param A vector of K unrestricted variables
 * @return Unit length vector
 **/
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> unit_vector_constrain(
    const var_value<T>& A) {
  using std::sqrt;
  const double r = sqrt(dot_self(A.val()));
  return make_callback_var(
      elt_divide(A.val(), r),
      [A, r](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_divide(res.adj(), r)
                   - A.val() * (dot_product(A.val(), res.adj()) / (r * r * r));
      });
}

/**
 * Return the unit length vector corresponding to the free vector y.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 * @tparam T type of the value of input matrix
 * @param A vector of K unrestricted variables
 * @param lp Log probability reference to increment.
 * @return Unit length vector
 **/
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> unit_vector_constrain(const var_value<T>& A,
                                                          var& lp) {
  using std::sqrt;
  double r = dot_self(A.val());
  lp -= 0.5 * r;
  r = sqrt(r);
  return make_callback_var(
      elt_divide(A.val(), r),
      [A, r, lp](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_divide(res.adj(), r)
                   - A.val() * (dot_product(A.val(), res.adj()) / (r * r * r))
                   - A.val() * lp.adj();
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
