#ifndef STAN_MATH_MANUAL_FORWARD_DECLS_HPP
#define STAN_MATH_MANUAL_FORWARD_DECLS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/fwd_decl.hpp>
#include <stan/math/rev/core/var_value_fwd_declare.hpp>
#include <stan/math/prim/meta.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#ifndef TBB_INTERFACE_NEW
#include <tbb/tbb_stddef.h>

#if TBB_VERSION_MAJOR >= 2020
#define TBB_INTERFACE_NEW
#endif
#endif

#ifdef TBB_INTERFACE_NEW
#include <tbb/global_control.h>
#include <tbb/task_arena.h>
#else
#include <tbb/task_scheduler_init.h>
#endif

#include <cstdlib>
#include <thread>
#include <complex>
#include <vector>
#include <type_traits>

namespace stan {
namespace math {
template <typename Ta1, typename Ta2, typename Tb, typename Tz,
          require_all_stan_scalar_t<Ta1, Ta2, Tb, Tz>* = nullptr,
          require_any_var_t<Ta1, Ta2, Tb, Tz>* = nullptr>
inline return_type_t<Ta1, Ta1, Tb, Tz> hypergeometric_2F1(const Ta1& a1,
                                                          const Ta2& a2,
                                                          const Tb& b,
                                                          const Tz& z);
}
}
namespace std {

/**
 * Specialization of the standard library complex number type for
 * reverse-mode autodiff type `stan::math::var`.
 */
template <>
class complex<stan::math::var>;
}

#endif
