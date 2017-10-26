#ifndef TEST_UNIT_MATH_PRIM_SCAL_META_APPLY_TEMPLATE_PERMUTATIONS_HPP
#define TEST_UNIT_MATH_PRIM_SCAL_META_APPLY_TEMPLATE_PERMUTATIONS_HPP

#include <tuple>

/*
 * apply_template_permutations_helper is the entry point for the template recursion
 * needed for apply_template_permutations. It calls func<type[I], type[J], type[K]>(rig),
 * and continues the recursion by decrementing K.
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_param Parameter to pass to functor
 * tparam I Loop index for first type
 * tparam J Loop index for second type
 * tparam K Loop index for third type
 */
template<typename T_typelist, typename T_functor, typename T_param, int I,
         int J, int K>
struct apply_template_permutations_helper {
  void operator()(const T_functor& func, const T_param& param) const {
    func.template operator()<typename std::tuple_element<I, T_typelist>::type,
                             typename std::tuple_element<J, T_typelist>::type,
                             typename std::tuple_element<K, T_typelist>::type>
      (param);

    apply_template_permutations_helper
      <T_typelist, T_functor, T_param, I, J, K - 1>{}(func, param);
  }
};

/*
 * This edge case catches when the paramht-most argument has completed one
 * iteration through all types. It computes func<type[I], type[J], type[0]>,
 * carries from the second argument, and continues the recursion.
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_param Parameter to pass to functor
 * tparam I Loop index for first type
 * tparam J Loop index for second type
 * tparam K Loop index for third type
 */
template<typename T_typelist, typename T_functor, typename T_param, int I,
         int J>
struct apply_template_permutations_helper<T_typelist, T_functor, T_param, I, J,
                                          0> {
  void operator()(const T_functor& func, const T_param& param) const {
    func.template operator()<typename std::tuple_element<I, T_typelist>::type,
                             typename std::tuple_element<J, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type>
      (param);

    apply_template_permutations_helper
      <T_typelist, T_functor, T_param, I, J - 1,
       std::tuple_size<T_typelist>::value - 1>{}(func, param);
  }
};

/*
 * This edge case catches when the second argument has completed one iteration
 * through all types. It computes func<type[I], type[0], type[0]>, carries from
 * the left-most argument and continues the recursion.
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_param Parameter to pass to functor
 * tparam I Loop index for first type
 */
template<typename T_typelist, typename T_functor, typename T_param, int I>
struct apply_template_permutations_helper<T_typelist, T_functor, T_param, I, 0,
                                          0> {
  void operator()(const T_functor& func, const T_param& param) const {
    func.template operator()<typename std::tuple_element<I, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type>
      (param);

    apply_template_permutations_helper<T_typelist, T_functor, T_param, I - 1,
                             std::tuple_size<T_typelist>::value - 1,
                             std::tuple_size<T_typelist>::value - 1>{}(func,
                                                                       param);
  }
};

/*
 * This edge case catches when all types have been iterated through for each
 * template argument. It computes func<type[0], type[0], type[0]> and ends the
 * recursion.
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_param Parameter to pass to functor
 */
template<typename T_typelist, typename T_functor, typename T_param>
struct apply_template_permutations_helper<T_typelist, T_functor, T_param, 0, 0,
                                          0> {
  void operator()(const T_functor& func, const T_param& param) const {
    func.template operator()<typename std::tuple_element<0, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type,
                             typename std::tuple_element<0, T_typelist>::type>
      (param);
  }
};


/*
 * apply_template_permutations uses template recursion to call the functor func
 * with all 3-element permutations of the types in T_typelist
 *
 * Roughly, this corresponds to:
 *  for (i in (I - 1):0) // indexes go backwards, like: I - 1, I - 2, ... 0
 *   for (j in (J - 1):0)
 *    for (k in (K - 1):0)
 *     func<type[i], type[j], type[k]>(param)
 *
 * tparam T_typelist Tuple templated with list of types to iterate through
 * tparam T_functor Templated functor
 * tparam T_param Parameter to pass to functor
 */
template<typename T_typelist, typename T_functor, typename T_param>
void apply_template_permutations(const T_functor& func, const T_param& param) {
  apply_template_permutations_helper
    <T_typelist, T_functor, T_param,
     std::tuple_size<T_typelist>::value - 1,
     std::tuple_size<T_typelist>::value - 1,
     std::tuple_size<T_typelist>::value - 1>{}(func, param);
}

#endif
