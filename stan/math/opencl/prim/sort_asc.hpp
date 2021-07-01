#ifndef STAN_MATH_OPENCL_PRIM_SORT_ASC_HPP
#define STAN_MATH_OPENCL_PRIM_SORT_ASC_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/mergesort.hpp>
#include <algorithm>

#include <stan/math/opencl/copy.hpp>

namespace stan {
namespace math {

/**
 */
template <typename T,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto sort_asc(T&& input) {
  using T_val = value_type_t<T>;
  matrix_cl<T_val> a = std::forward<T>(input);
  matrix_cl<T_val> b(a.rows(), a.cols());
  matrix_cl<T_val>* in_ptr = &a;
  matrix_cl<T_val>* out_ptr = &b;
  for (int run_len = 1; run_len < input.size(); run_len *= 2) {
    //      int run_len =(input.size() + runs - 1) / runs;
    int runs = (input.size() + run_len - 1) / run_len;
    //    std::cout << "runs " << runs << std::endl;
    opencl_kernels::merge_step(cl::NDRange(std::min((runs + 1) / 2, 32*15*16)), *out_ptr, *in_ptr,
                               run_len, input.size(), (runs + 1) / 2);
    //    std::cout << stan::math::from_matrix_cl(*out_ptr) << std::endl <<
    //    std::endl;
    std::swap(in_ptr, out_ptr);
//    if (run_len == 32*32) {
//      break;
//    }
  }
  return *in_ptr;
}
}  // namespace math
}  // namespace stan
#endif
#endif
