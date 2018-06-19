#ifndef STAN_MATH_PRIM_MAT_META_EIGEN_SCALARBINARYOPTRAITS_HPP
#define STAN_MATH_PRIM_MAT_META_EIGEN_SCALARBINARYOPTRAITS_HPP

#include <stan/math/prim/scal/meta/is_arith_like.hpp>
#include <stan/math/prim/scal/meta/is_complex.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/rm_complex.hpp>
#include <stan/math/prim/scal/meta/rm_zeroing.hpp>
#include <Eigen/Dense>
#include <complex>
#include <type_traits>

/**
 * Eigen scalar op traits specialization for complex variables
 */
namespace Eigen {

template <class T1, class T2, template <class, class> class OP>
struct ScalarBinaryOpTraits<
    T1,
    std::enable_if_t<
        // !VectorBlock !is_eigen
        (stan::is_complex<T1>::value || stan::is_arith_like<T1>::value)
            && (stan::is_complex<T2>::value || stan::is_arith_like<T2>::value)
            && !std::is_same<T1, T2>::value &&  // avoid Eigen's template
            ((stan::is_complex<T1>::value &&    // next boolean avoids Eigen
              !std::is_same<stan::rm_complex_t<T1>,
                            stan::rm_zeroing_t<T2>>::value)
             || (stan::is_complex<T2>::value &&  // next bool avoids Eigen
                 !std::is_same<stan::rm_zeroing_t<T1>,
                               stan::rm_complex_t<T2>>::value)),
        T2>,
    OP<T1, T2>> {
  typedef std::complex<typename stan::return_type<stan::rm_complex_t<T1>,
                                                  stan::rm_complex_t<T2>>::type>
      ReturnType;
};

}  // namespace Eigen
#endif
