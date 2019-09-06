#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_TYPE_STR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_TYPE_STR_HPP
#ifdef STAN_OPENCL

#include <string>

/**
 * Determins a string name of a type. Unsupported types fail static assert.
 */
template<typename T>
struct type_str{
    static_assert(sizeof(T)==-1, "Unsupported type in type_str");
};

#define ADD_TYPE_TO_TYPE_STR(t) template<> struct type_str<t>{static const std::string name;}; const std::string type_str<t>::name(#t);
ADD_TYPE_TO_TYPE_STR(double)
ADD_TYPE_TO_TYPE_STR(int)
ADD_TYPE_TO_TYPE_STR(bool)
#undef ADD_TYPE_TO_TYPE_STR

#endif
#endif
