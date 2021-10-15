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

template <typename T>
using eigen_base_filter_t = typename internal::global_math_functions_filtering_base<T>::type;