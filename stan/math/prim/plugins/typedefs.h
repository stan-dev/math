#ifndef STAN_MATH_PRIM_PLUGINS_TYPEDEFS_H
#define STAN_MATH_PRIM_PLUGINS_TYPEDEFS_H

/**
 * Reimplements is_fvar without requiring external math headers
 *
 * decltype((void)(T::d_)) is a pre C++17 replacement for
 * std::void_t<decltype(T::d_)>
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template <typename, typename = void>
struct is_fvar : std::false_type { };
template <typename T>
struct is_fvar<T, decltype((void)(T::d_))> : std::true_type { };

/**
 * Reimplements is_var without requiring external math headers
 */
template <typename, typename = void>
struct is_var : std::false_type { };
template <typename T>
struct is_var<T, decltype((void)(std::decay_t<T>::vi_))> : std::true_type { };

/**
 * Reimplements is_vari without requiring external math headers
 */
template <typename, typename = void>
struct is_vari : std::false_type { };
template <typename T>
struct is_vari<T, decltype((void)(std::decay_t<std::remove_pointer_t<T>>::adj_))> : std::true_type
{ };

/**
 * Struct for determining the appropriate return type for a given input type.
 * The type mappings are:
 * - arithmetic -> arithmetic
 * - var_value<T> -> T
 * - vari_value<T> -> T
 * - fvar<T> -> T
 */
template <typename T, typename Enable = void>
struct val_return { };

template <typename T>
struct val_return<T, std::enable_if_t<std::is_arithmetic<T>::value>> {
  using type = T;
};

template <typename T>
struct val_return<T, std::enable_if_t<is_var<T>::value>> {
  using type = typename T::value_type;
};

template <typename T>
struct val_return<T, std::enable_if_t<is_vari<T>::value>> {
  using type = typename std::remove_pointer_t<T>::value_type;
};

template <typename T>
struct val_return<T, std::enable_if_t<is_fvar<T>::value>> {
  using type = typename T::Scalar;
};

template <typename T>
using val_return_t = typename val_return<std::decay_t<T>>::type;

// Typedef for determining the type within a vari_value<T> and returning
// with a pointer type
template <typename T>
using vi_return_t = std::add_pointer_t<typename T::vari_type>;

#endif
