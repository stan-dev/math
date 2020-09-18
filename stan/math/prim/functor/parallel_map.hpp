#ifndef STAN_MATH_PRIM_FUNCTOR_PARALLEL_MAP_HPP
#define STAN_MATH_PRIM_FUNCTOR_PARALLEL_MAP_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace stan {
namespace math {

template <typename ApplyFunction, typename IndexFunction,
          typename Res, typename ArgsTuple,
          require_st_arithmetic<Res>* = nullptr>
inline decltype(auto) parallel_map(const ApplyFunction& app_fun,
                                   const IndexFunction& index_fun,
                                   Res&& result, ArgsTuple&& x) {
  tbb::parallel_for(
    tbb::blocked_range<size_t>(0, result.size()), 
    [&x,&app_fun,&index_fun,&result](
     const tbb::blocked_range<size_t>& r) {
      for (size_t i = r.begin(); i < r.end(); ++i) {
        // Apply specified function to arguments at current iteration
        result(i) = apply(
            [&](auto&&... args) {
              return index_fun(i, app_fun, args...);
            }, std::forward<ArgsTuple>(x));
      }
    });
  return std::forward<Res>(result);
}

}  // namespace math
}  // namespace stan
#endif
