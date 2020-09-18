#ifndef STAN_MATH_REV_CORE_PARALLEL_MAP_HPP
#define STAN_MATH_REV_CORE_PARALLEL_MAP_HPP

#include <stan/math/rev/fun/typedefs.hpp>
#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace stan {
namespace math {

template <typename ApplyFunction, typename IndexFunction,
          typename Res, typename ArgsTuple>
inline decltype(auto) parallel_map(const ApplyFunction& app_fun,
                                   const IndexFunction& index_fun,
                                   Res&& result, ArgsTuple&& x) {
    // Functors for manipulating vars at a given iteration of the loop
    auto var_counter = [&](auto&... xargs) {
      return count_vars(xargs...);
    };
    auto make_tuple_refs = [&](const auto... xargs){
      return std::make_tuple(xargs...);
    };
    auto var_copier = [&](const auto&... xargs) {
      return deep_copy_vars(xargs...);
    };


    int S = result.size();
    // Assuming that the number of the vars at each iteration of the loop is
    // the same, as the operations at each iteration should be the same, can
    // just count vars at the first iteration.
    int nvars = apply(
      [&](auto&&... args) { return index_fun(0, var_counter, args...); }, x);

    vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(
      S * nvars);
    double* values = ChainableStack::instance_->memalloc_.alloc_array<double>(
      S);
    double* partials = ChainableStack::instance_->memalloc_.alloc_array<double>(
      S * nvars);
    Eigen::Map<Eigen::VectorXd>(partials, S * nvars).setZero();

    tbb::parallel_for(
      tbb::blocked_range<size_t>(0, S), 
      [&x,&partials,&values,&app_fun,
       &var_copier,&index_fun,&varis,&nvars,&make_tuple_refs](
       const tbb::blocked_range<size_t>& r) {
        // Run nested autodiff in this scope
        nested_rev_autodiff nested;

        for (size_t i = r.begin(); i < r.end(); ++i) {
          // Create tuple of args at current iteration
          auto local_vars = apply(
            [&](auto&&... args) {
              return index_fun(i, make_tuple_refs, args...);
            }, x);

          // Save varis for vars at current iteration
          apply([&](auto&&... args) { save_varis(varis + nvars*i, args...); },
                local_vars);

          // Create nested autodiff copies of all arguments at current iteration
          // that do not point back to main autodiff stack
          auto args_tuple_local_copy = apply(
          [&](auto&&... args) {
            return std::tuple<decltype(deep_copy_vars(args))...>(
                deep_copy_vars(args)...);
          }, local_vars);

          // Apply specified function to arguments at current iteration
          var out = apply(
              [&](auto&&... args) {
                return app_fun(args...);
              }, args_tuple_local_copy);

          out.grad();

          // Extract values and adjoints to be put into vars on main
          // autodiff stack
          values[i] = out.val();
          apply([&](auto&&... args) {
            accumulate_adjoints(partials + nvars*i,
                          std::forward<decltype(args)>(args)...); },
            std::move(args_tuple_local_copy));

        }
      });

  for(int i = 0; i < S; ++i) {
    result.coeffRef(i) = var(new precomputed_gradients_vari(
      values[i],
      nvars,
      varis + nvars*i,
      partials + nvars*i));
  }

  return std::forward<Res>(result);
}

}  // namespace math
}  // namespace stan
#endif
