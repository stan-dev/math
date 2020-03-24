#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>

#include <stan/math/opencl/kernel_generator/load.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/kernel_generator/binary_operation.hpp>
#include <stan/math/opencl/kernel_generator/unary_function_cl.hpp>
#include <stan/math/opencl/kernel_generator/block.hpp>
#include <stan/math/opencl/kernel_generator/select.hpp>
#include <stan/math/opencl/kernel_generator/rowwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/colwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/transpose.hpp>

#include <stan/math/opencl/kernel_generator/multi_result_kernel.hpp>
#include <stan/math/opencl/kernel_generator/get_kernel_source_for_evaluating_into.hpp>
#include <stan/math/opencl/kernel_generator/evaluate_into.hpp>

#include <stan/math/opencl/kernel_generator/matrix_cl_conversion.hpp>

#endif
#endif
