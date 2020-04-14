/**
 * Reimplements is_fvar without requiring external math headers
 *
 * decltype((void)(T::d_)) is a pre C++17 replacement for
 * std::void_t<decltype(T::d_)>
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template<class, class = void>
struct is_fvar : std::false_type
{ };
template<class T>
struct is_fvar<T, decltype((void)(T::d_))> : std::true_type
{ };

template<class, class = void>
struct is_var : std::false_type
{ };
template<class T>
struct is_var<T, decltype((void)(T::vi_))> : std::true_type
{ };

//TODO(Andrew): Replace std::is_const<>::value with std::is_const_v<> after move to C++17
template<typename T>
using double_return_t = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                                         const double,
                                         double>;
template<typename T>
using reverse_return_t = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                                         const double&,
                                         double&>;

template<typename T>
using vari_return_t = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                                          const decltype(T::vi_)&,
                                          decltype(T::vi_)&>;

template<typename T>
using forward_return_t = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                                         const decltype(T::val_)&,
                                         decltype(T::val_)&>;

template <typename T>
using require_double_t = std::enable_if_t<std::is_floating_point<T>::value,double_return_t<T>>;

template <typename T>
using require_var_t = std::enable_if_t<is_var<T>::value,reverse_return_t<T>>;

template <typename T>
using require_fvar_t = std::enable_if_t<is_fvar<T>::value,forward_return_t<T>>;

template <typename T>
using require_vari_t = std::enable_if_t<std::is_pointer<T>::value,reverse_return_t<T>>;
