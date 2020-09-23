#ifndef STAN_MATH_REV_FUNCTOR_PARALLEL_MAP_HPP
#define STAN_MATH_REV_FUNCTOR_PARALLEL_MAP_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace stan {
namespace math {

template <bool Ranged, typename Res, typename ApplyFunction, typename IndexFunction,
          require_var_t<return_type_t<Res>>* = nullptr,
          typename... Args>
inline auto parallel_map(ApplyFunction&& app_fun, IndexFunction&& index_fun,
   int grainsize, Args&&... args) {
  // Functors for manipulating vars at a given iteration of the loop
  auto var_counter = [&](auto&... xargs) {
   return count_vars(xargs...);
  };
  auto var_copier = [&](auto&&... xargs) {
   return std::tuple<decltype(deep_copy_vars(xargs))...>(
     deep_copy_vars(xargs)...);
  };
  auto vari_saver = [&](int i, int nvars, vari** varis) {
   return [=](const auto&... xargs) {
     save_varis(varis + nvars*i, xargs...);
   };
  };
  std::decay_t<decltype(app_fun(args...))> result(max_size(args...));
  int S = result.size();
  // Assuming that the number of the vars at each iteration of the loop is
  // the same (as the operations at each iteration should be the same), we can
  // just count vars at the first iteration.
  int nvars = index_fun(0, var_counter, args...);

  vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(
   S * nvars);
  double* values = ChainableStack::instance_->memalloc_.alloc_array<double>(
   S);
  double* partials = ChainableStack::instance_->memalloc_.alloc_array<double>(
   S * nvars);
  Eigen::Map<Eigen::VectorXd>(partials, S * nvars).setZero();

  tbb::parallel_for(
   tbb::blocked_range<size_t>(0, S, grainsize),
   [&](
    const tbb::blocked_range<size_t>& r) {
     // Run nested autodiff in this scope

     for (size_t i = r.begin(); i < r.end(); ++i) {
       nested_rev_autodiff nested;
       // Save varis from arguments at current iteration
       index_fun(i, vari_saver(i, nvars, varis), args...);
       // Create nested autodiff copies of all arguments at current
       // iteration that do not point back to main autodiff stack
       auto args_tuple_local_copy = index_fun(i, var_copier, args...);

       // Apply specified function to arguments at current iteration
       var out = apply(
           [&](auto&&... args) {
             return app_fun(args...);
           }, args_tuple_local_copy);

       out.grad();

       // Extract value and adjoints to be put into vars on main
       // autodiff stack
       values[i] = out.vi_->val_;
       apply([&](auto&&... args) {
         accumulate_adjoints(partials + nvars*i,
                             std::forward<decltype(args)>(args)...); },
         args_tuple_local_copy);
     }
   });
  // Pack values and adjoints into new vars on main autodiff stack
  for(int i = 0; i < S; ++i) {
  result.coeffRef(i) = var(new precomputed_gradients_vari(
   values[i],
   nvars,
   varis + nvars*i,
   partials + nvars*i));
  }
  return result;
  }

  }  // namespace math
  }  // namespace stan
  #endif
