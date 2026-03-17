#ifndef STAN_MATH_PRIM_FUNCTOR_FILTER_MAP_HPP
#define STAN_MATH_PRIM_FUNCTOR_FILTER_MAP_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/make_holder_tuple.hpp>
#include <stan/math/prim/functor/tuple_concat.hpp>
#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

namespace internal {

template <template <typename...> class Filter, typename T>
struct inspect_tuple {
  static constexpr bool value = Filter<T>::value;
};

template <template <typename...> class Filter, typename... Types>
struct inspect_tuple<Filter, std::tuple<Types...>> {
  static constexpr bool value = Filter<std::tuple<Types...>>::value
                                || (inspect_tuple<Filter, Types>::value || ...);
};

template <template <typename...> class Filter, typename T, typename... VecArgs>
struct inspect_tuple<Filter, std::vector<T, VecArgs...>> {
  static constexpr bool value
      = Filter<std::vector<T, VecArgs...>>::value
        || inspect_tuple<Filter, std::decay_t<T>>::value;
};

/**
 * Check if a tuple or type contains a tuple that passes the filter.
 * @tparam Filter a struct that accepts one template parameter and has a static
 *  constexpr bool member named value that is true if the type should be
 *  included in the output tuple.
 * @tparam T type to check
 */
template <template <typename...> class Filter, typename T>
inline constexpr bool inspect_tuple_v
    = internal::inspect_tuple<Filter, std::decay_t<T>>::value;

/**
 * Filter a tuple and apply a functor to each element that passes the filter.
 * @note The `Filter` must have a static constexpr bool member named `value`
 * that is true if the type should be included in the output tuple.
 * Note that this function automatically inspects into tuples and
 * `std::vector<T>::value_type`'s. The `filter_map` will recursively apply
 * itself to inner containers as long as it sees a tuple in type type.
 *  So for instance if your type is a
 * `tuple<vector<tuple<vector<vector<double>>>>` your functor `f` must support
 * operationg on `vector<vector<double>>` types.
 * @tparam Filter a struct that accepts one template parameter and has a static
 *  constexpr bool member named value that is true if the type should be
 *  included in the output tuple.
 * @tparam InVector For internal use. If true then we assume we are inside of a
 *  `std::vector` and the return type should not be wrapped in a tuple.
 * @tparam InTuple For internal use. If true then we assume we are inside of a
 *  tuple and any subtuples should be double wrapped so that tuple_concat
 *  produces a tuple for this element.
 * @tparam F Type of functor
 * @tparam T Any type
 * @param f functor callable
 * @param x Any type
 * @return a tuple with the functor applied to each element which passed the
 * filter.
 */
template <template <typename...> class Filter, bool InVector = false,
          bool InTuple = false, typename F, typename T>
inline constexpr decltype(auto) filter_map(F&& f, T&& x) {
  if constexpr (inspect_tuple_v<Filter, T>) {
    if constexpr (is_tuple_v<T>) {
      auto ret = stan::math::apply(
          [&f](auto&&... args) {
            return stan::math::tuple_concat(filter_map<Filter, false, true>(
                f, std::forward<decltype(args)>(args))...);
          },
          std::forward<T>(x));
      /**
       * If we are in at this stage, we want tuple_concat to return a tuple here
       * So we return a tuple(tuple()) so that tuple_cat concats
       * the first layer of tuple.
       * For example, if our input is a tuple(double, tuple(double,
       * vec<double>)) with an identity filter we want tuple_concat to return a
       * tuple(double, tuple(double, vec<double>)).
       * Without the double tuple we would get back a tuple(double, double,
       * vec<double>).
       */
      if constexpr (InTuple) {
        return make_holder_tuple(std::move(ret));
      } else {
        return ret;
      }
    } else if constexpr (is_std_vector_v<T>) {
      /* 3 cases for vectors
       * 1. value_type is a tuple
       * 2. value_type is a scalar or Eigen matrix
       * 3. value_type is a std::vector which can hold either (1) or (2)
       */
      if constexpr (contains_tuple<T>::value) {
        std::vector<decltype(filter_map<Filter, true>(f, x[0]))> ret;
        for (size_t i = 0; i < x.size(); ++i) {
          ret.push_back(filter_map<Filter, true>(f, x[i]));
        }
        /*
         * If we are in a vector, return the raw type, otherwise we are in
         * a tuple and we want to return a tuple of the vector.
         */
        if constexpr (InVector) {
          return ret;
        } else {
          return std::make_tuple(std::move(ret));
        }
      } else {
        if constexpr (InVector) {
          return std::forward<F>(f)(std::forward<T>(x));
        } else {
          return make_holder_tuple(std::forward<F>(f)(std::forward<T>(x)));
        }
      }
    } else {
      if constexpr (InVector) {
        return std::forward<F>(f)(std::forward<T>(x));
      } else {
        return make_holder_tuple(std::forward<F>(f)(std::forward<T>(x)));
      }
    }
  } else {
    return std::make_tuple();
  }
}
}  // namespace internal
/**
 * Filter a tuple and apply a functor to each element that passes the filter.
 * @note The `Filter` will only check `T` and if `T` is a tuple, it will
 * recursively check each element of the tuple. But it will not inspect into
 * `std::vector` elements automatically. If you want to inspect the inner
 * element of an `std::vector` your type trait must do that itself.
 * @tparam Filter a struct that accepts one template parameter and has a static
 *  constexpr bool member named value that is true if the type should be
 *  included in the output tuple.
 * @tparam F Type of functor
 * @tparam T A tuple
 * @param f functor callable
 * @param x tuple of arguments
 * @return a tuple with the functor applied to each element which passed the
 * filter.
 */
template <template <typename...> class Filter, typename F, typename T,
          require_tuple_t<T>* = nullptr>
inline constexpr decltype(auto) filter_map(F&& f, T&& x) {
  return internal::filter_map<Filter>(std::forward<F>(f), std::forward<T>(x));
}
}  // namespace math
}  // namespace stan

#endif
