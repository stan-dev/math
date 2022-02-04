#ifndef STAN_MATH_PRIM_FUNCTOR_FILTER_MAP_HPP
#define STAN_MATH_PRIM_FUNCTOR_FILTER_MAP_HPP

#include <functional>
#include <tuple>
#include <utility>
#include <stan/math/prim/functor/apply.hpp>

namespace stan {
namespace math {


  namespace internal {
  /**
   * Evaluator for true and false cases of filtering when calling `filter_map`
   * @tparam IncludeElem
   */
  template <bool IncludeElem>
  struct filter_map_evaluator;

  template <>
  struct filter_map_evaluator<true> {
      /**
       * Base case with no more arguments left to check so the accumulated tuple
       * is returned.
       */
      template <template <typename> class Filter, typename F,
                typename... AccArgs, typename PassArg>
      inline static auto run(F&& f, std::tuple<AccArgs...>&& acc_tuple,
                             PassArg&& arg) {
          return std::tuple_cat(std::move(acc_tuple),
                                std::make_tuple(f(std::forward<PassArg>(arg))));
      }

      /*
       * `true` path where the `f` is applied to the element and it's result is
       * added to the result tuple
       *
       * @tparam Filter a struct that accepts one template parameter and has a
       * static constexpr bool member named value
       * @tparam F Type of functor
       * @tparam AccArgs Parameter pack of arguments that will be returned in the
       * final tuple
       * @tparam PassArg First template parameter peeled off from `OtherArgs` that
       * passed the check
       * @tparam CheckArg Second template parameter peeled off from `OtherArgs` to
       * check
       * @tparam OtherArgs Rest of the parameter pack of arguments
       * @param f functor callable
       * @param acc_tuple Tuple of arguments that passed filter and had `f`
       * applied to them.
       * @param arg Argument peeled from parameter pack that can now be added to
       * final tuple
       * @param curr_arg The current argument to be evaluated
       * @param args objects to be filtered and mapped
       */
      template <template <typename> class Filter, typename F,
                typename... AccArgs, typename PassArg, typename CheckArg,
                typename... OtherArgs>
      inline static decltype(auto) run(F&& f, std::tuple<AccArgs...>&& acc_tuple,
                                       PassArg&& arg, CheckArg&& next_arg,
                                       OtherArgs&&... args) {
          return filter_map_evaluator<Filter<CheckArg>::value>::template run<
              Filter>(
              std::forward<F>(f),
              std::tuple_cat(std::move(acc_tuple),
                             std::make_tuple(f(std::forward<PassArg>(arg)))),
              std::forward<CheckArg>(next_arg), std::forward<OtherArgs>(args)...);
      }
  };

  template <>
  struct filter_map_evaluator<false> {
      /**
       * Base case with no more arguments left to check so the accumulated tuple
       * is returned.
       */
      template <template <typename> class Filter, typename F,
                typename... AccArgs, typename FailArg>
      inline static auto run(F&& /* f */, std::tuple<AccArgs...>&& acc_tuple,
                             FailArg&& /* arg */) {
          return std::move(acc_tuple);
      }

      /*
       * `false` path where we do not add the argument to the end parameter pack
       *
       * @tparam Filter a struct that accepts one template parameter and has a
       * static constexpr bool member named value
       * @tparam F Type of functor
       * @tparam AccArgs Parameter pack of arguments that will be returned in the
       * final tuple
       * @tparam FailArg First template parameter peeled off from `OtherArgs` that
       * failed the check
       * @tparam CheckArg Second template parameter peeled off from `OtherArgs` to
       * check
       * @tparam OtherArgs Rest of the parameter pack of arguments
       * @param f functor callable
       * @param acc_tuple Tuple of arguments that passed filter and had `f`
       * applied to them.
       * @param curr_arg The current argument to be evaluated
       * @param args objects to be filtered and mapped
       */
      template <template <typename> class Filter, typename F,
                typename... AccArgs, typename FailArg, typename CheckArg,
                typename... OtherArgs>
      inline static decltype(auto) run(F&& f, std::tuple<AccArgs...>&& acc_tuple,
                                       FailArg&& /* arg */, CheckArg&& curr_arg,
                                       OtherArgs&&... args) {
          return filter_map_evaluator<Filter<CheckArg>::value>::template run<
              Filter>(std::forward<F>(f), std::move(acc_tuple),
                      std::forward<CheckArg>(curr_arg),
                      std::forward<OtherArgs>(args)...);
      }
  };

  /**
   * Base case for no arguments to filter which returns an empty tuple. This is
   * only called when `filter_map` receives no inputs to traverse so it can just
   * return an empty tuple.
   * @tparam Filter a struct that accepts one template parameter and has a static
   *   constexpr bool member named value
   * @tparam F Type of functor
   * @tparam AccArgs A parameter pack of arguments
   */
  template <template <typename> class Filter, typename F, typename... AccArgs>
  inline constexpr decltype(auto) filter_map_impl(
      F&& /* */, std::tuple<AccArgs...>&& /* acc_tuple */) {
      return std::tuple<>{};
  }

  /*
   * Subset a parameter pack by a compile time filter on the types and return a
   * tuple with an unary functor applied to each element where the filter was
   * true:
   *
   * @tparam Filter a struct that accepts one template parameter and has a static
   *   constexpr bool member named value
   * @tparam F Type of functor
   * @tparam AccArgs A parameter pack of arguments that will be returned in the
   * final tuple
   * @tparam CheckArg The first parameter in the parameter pack of `filter_map`
   * @tparam OtherArgs A parameter pack of arguments
   * @param f functor callable
   * @param acc_tuple Tuple of arguments that passed filter and had `f` applied to
   * them.
   * @param check_arg The first argument in the parameter pack given to
   * `filter_map`.
   * @param args objects to be filtered and mapped
   */
  template <template <typename> class Filter, typename F, typename... AccArgs,
            typename CheckArg, typename... OtherArgs>
  inline constexpr decltype(auto) filter_map_impl(
      F&& f, std::tuple<AccArgs...>&& acc_tuple, CheckArg&& first_arg,
      OtherArgs&&... args) {
      return filter_map_evaluator<Filter<CheckArg>::value>::template run<Filter>(
          std::forward<F>(f), std::move(acc_tuple),
          std::forward<CheckArg>(first_arg), std::forward<OtherArgs>(args)...);
  }

  }  // namespace internal

  /*
   * Subset a parameter pack by a compile time filter on the types and return a
   * tuple with an unary functor applied to each element where the filter was
   * true:
   *
   * @tparam Filter a struct that accepts one template parameter and has a static
   *   constexpr bool member named value
   * @tparam F Type of functor
   * @tparam Args A parameter pack of arguments
   * @param f functor callable
   * @param args parameter pack of args
   */
  template <template <typename> class Filter, typename F, typename... Args>
  inline constexpr decltype(auto) filter_map(F&& f, Args&&... args) {
      return internal::filter_map_impl<Filter>(std::forward<F>(f), std::tuple<>{},
                                               std::forward<Args>(args)...);
  }

  /*
   * Subset a tuple by a compile time filter on the types and return a
   * tuple with an unary functor applied to each element where the filter was
   * true:
   *
   * @tparam Filter a struct that accepts one template parameter and has a static
   *   constexpr bool member named value
   * @tparam F Type of functor
   * @tparam Args A parameter pack of arguments
   * @param f functor callable
   * @param args parameter pack of args
   */
  template <template <typename> class Filter, typename F, typename... Args>
  inline constexpr decltype(auto) filter_map(F&& f, std::tuple<Args...>&& args) {
      return stan::math::apply(
          [&f](auto&&... args) {
              return internal::filter_map_impl<Filter>(f, std::tuple<>{},
                                                       std::move(args)...);
          },
          std::move(args));
  }

  /*
   * Subset a tuple by a compile time filter on the types and return a
   * tuple with an unary functor applied to each element where the filter was
   * true:
   *
   * @tparam Filter a struct that accepts one template parameter and has a static
   *   constexpr bool member named value
   * @tparam F Type of functor
   * @tparam Args A parameter pack of arguments
   * @param f functor callable
   * @param args parameter pack of args
   */
  template <template <typename> class Filter, typename F, typename... Args>
  inline constexpr decltype(auto) filter_map(F&& f,
                                             const std::tuple<Args...>& args) {
      return stan::math::apply(
          [&f](auto&&... args) {
              return internal::filter_map_impl<Filter>(f, std::tuple<>{},
                                                       args...);
          },
          args);
  }

}  // namespace math
}  // namespace stan

#endif
