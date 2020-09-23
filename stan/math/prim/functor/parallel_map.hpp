#ifndef STAN_MATH_PRIM_FUNCTOR_PARALLEL_MAP_HPP
#define STAN_MATH_PRIM_FUNCTOR_PARALLEL_MAP_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace stan {
namespace math {

template <bool Ranged, typename Res, typename ApplyFunction, typename IndexFunction,
          require_arithmetic_t<return_type_t<Res>>* = nullptr,
          std::enable_if_t<!Ranged>* = nullptr,
          typename... Args>
inline auto parallel_map(ApplyFunction&& app_fun,
                                   IndexFunction&& index_fun, int grainsize,
                                   Args&&... args) {
  std::decay_t<decltype(app_fun(args...))> result(max_size(args...));
  tbb::parallel_for(
    tbb::blocked_range<size_t>(0, result.size(), grainsize),
    [&](const tbb::blocked_range<size_t>& r) {
      for (size_t i = r.begin(); i < r.end(); ++i) {
        // Apply specified function to arguments at current iteration
        result(i) =  index_fun(i, app_fun, args...);
      }
    });
    return result;
}

template <bool Ranged, typename Res, typename ApplyFunction, typename IndexFunction,
          require_arithmetic_t<return_type_t<Res>>* = nullptr,
          std::enable_if_t<Ranged>* = nullptr,
          typename... Args>
inline auto parallel_map(ApplyFunction&& app_fun,
                                   IndexFunction&& index_fun, int grainsize,
                                   Args&&... args) {
  std::decay_t<decltype(app_fun(args...))> result(max_size(args...));
  tbb::parallel_for(
    tbb::blocked_range<size_t>(0, result.size(), grainsize),
    [&](const tbb::blocked_range<size_t>& r) {
      result.segment(r.begin(), r.size()) =  index_fun(r.begin(), r.size(), app_fun, args...);
    });
    return result;
}


template <typename ApplyFunction, typename IndexFunction,
          typename... Args>
inline auto parallel_map(ApplyFunction&& app_fun,
                                   IndexFunction&& index_fun, int grainsize,
                                   Args&&... x) {
using ret_type = std::decay_t<decltype(app_fun(x...))>;
return parallel_map<true, ret_type>(std::forward<ApplyFunction>(app_fun), std::forward<IndexFunction>(index_fun),
                    grainsize, std::forward<Args>(x)...);
}
}  // namespace math
}  // namespace stan
#endif
