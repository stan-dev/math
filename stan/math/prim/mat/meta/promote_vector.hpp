#ifndef STAN_MATH_PRIM_MAT_META_PROMOTE_VECTOR_HPP
#define STAN_MATH_PRIM_MAT_META_PROMOTE_VECTOR_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>

namespace stan {
  /**
   * promote_vector is a metaprogram that returns the type of the first
   * vector-like type in the arguments list.
   * If there are no vector-like types, it returns the first scalar type in its
   * argument list promoted according to the boost type promotion rules.
   *
   * @tparam T1 First argument
   * @tparam T2 Second argument
   * @tparam Enable (used internally, don't specify)
   */

  template <typename T1, typename T2, typename Enable = void>
  struct promote_vector {
    typedef typename boost::math::tools::promote_args<T1>::type type;
  };

  template <typename T1, typename T2>
  struct promote_vector<T1, T2, typename boost::enable_if_c<
                                  is_vector<T1>::value &&
                                  !is_vector<T2>::value
                                  >::type > {
    typedef T1 type;
  };

  template <typename T1, typename T2>
  struct promote_vector<T1, T2, typename boost::enable_if_c<
                                  !is_vector<T1>::value &&
                                  is_vector<T2>::value
                                  >::type > {
    typedef T2 type;
  };
}
#endif

