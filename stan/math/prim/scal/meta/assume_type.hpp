#ifndef STAN_MATH_PRIM_SCAL_META_ASSUME_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_ASSUME_TYPE_HPP

#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

template <typename T_desired, typename T_actual, typename = std::enable_if_t<std::is_convertible<T_actual, T_desired>::value && !is_eigen<T_desired>::value>>
inline T_actual&& assume_type(T_actual&& a){
  return std::forward<T_actual>(a);
}

template <typename T_desired, typename T_actual, typename = std::enable_if_t<std::is_convertible<T_actual, T_desired>::value &&
    static_cast<int>(T_desired::RowsAtCompileTime) == static_cast<int>(T_actual::RowsAtCompileTime) &&
    static_cast<int>(T_desired::ColsAtCompileTime) == static_cast<int>(T_actual::ColsAtCompileTime)>, typename = void>
inline T_actual&& assume_type(T_actual&& a){
  return std::forward<T_actual>(a);
}

template <typename T_desired, typename T_actual, typename = std::enable_if_t<!std::is_convertible<T_actual, T_desired>::value>>
inline T_desired assume_type(const T_actual& a){
  throw std::runtime_error("Wrong type assumed! Please file a bug report.");
}

template <typename T_desired, typename T_actual, typename = std::enable_if_t<!std::is_convertible<T_actual, T_desired>::value ||
    static_cast<int>(T_desired::RowsAtCompileTime) != static_cast<int>(T_actual::RowsAtCompileTime) ||
    static_cast<int>(T_desired::ColsAtCompileTime) != static_cast<int>(T_actual::ColsAtCompileTime)>, typename = void>
inline T_desired assume_type(const T_actual& a){
  throw std::runtime_error("Wrong type assumed! Please file a bug report.");
}


}
}

#endif
