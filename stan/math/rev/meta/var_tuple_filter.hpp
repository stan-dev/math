#ifndef STAN_MATH_REV_META_CONDITIONAL_SEQUENCE_HPP
#define STAN_MATH_REV_META_CONDITIONAL_SEQUENCE_HPP

#include <stan/math/rev/meta/is_var.hpp>
namespace stan {
// NOTE: This should probably be in Stan namespace not math
namespace math {

struct tuple_cat_caller {
    template <typename... Ts>
   auto operator()(Ts&&... args) {
       return std::tuple_cat(std::forward<Ts>(args)...);
   }
};

template <template<typename...> class Pred, ,template <typename...> class Filter, typename T>
using var_filter_helper = std::conditional_t<Pred<T>::value,
  std::tuple<Filter<T>>, std::tuple<>>;

template <typename T>
using var_vari_value_t = ::stan::internal::get_var_vari_value_t<scalar_type_t<T>>;

template <typename T>
using container_var_vari_value_t = std::conditional_t<is_container<T>::value,
  ::stan::internal::get_var_vari_value_t<scalar_type_t<T>>**,
  ::stan::internal::get_var_vari_value_t<scalar_type_t<T>>*>;

template <typename... Ts>
using var_to_vari_filter_t = std::result_of_t<
 tuple_cat_caller(var_filter_helper<is_var_value, container_var_vari_value_t, Ts>...)>;


}
}
#endif
