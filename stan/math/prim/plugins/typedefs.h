#ifndef STAN_MATH_PRIM_PLUGINS_TYPEDEFS_H
#define STAN_MATH_PRIM_PLUGINS_TYPEDEFS_H

template<typename, typename = void>
struct is_fvar : std::false_type
{ };
template<typename T>
struct is_fvar<T, decltype((void)(T::d_))> : std::true_type
{ };

template<typename, typename = void>
struct is_var : std::false_type
{ };
template<typename T>
struct is_var<T, decltype((void)(std::decay_t<std::remove_pointer_t<T>>::vi_))> : std::true_type
{ };

template<typename, typename = void>
struct is_vari : std::false_type
{ };
template<typename T>
struct is_vari<T, decltype((void)(std::decay_t<std::remove_pointer_t<T>>::adj_))> : std::true_type
{ };

template <typename T, typename Enable = void>
struct val_return { };

template <typename T>
struct val_return<T, std::enable_if_t<!is_fvar<T>::value>> {
  using type = double;
};
template <typename T>
struct val_return<T, std::enable_if_t<is_fvar<T>::value>> {
  using type = decltype(T::d_);
};

template <typename T>
using val_return_t = typename val_return<std::decay_t<T>>::type;

template <typename T>
using vi_return_t = decltype(T::vi_);

template <typename T>
using d_return_t = decltype(T::d_);

#endif
