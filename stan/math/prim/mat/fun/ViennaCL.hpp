#ifndef STAN_MATH_PRIM_MAT_FUN_VIENNACL_HPP
#define STAN_MATH_PRIM_MAT_FUN_VIENNACL_HPP

#ifdef STAN_GPU
  #ifndef VIENNACL_WITH_OPENCL
    #define VIENNACL_WITH_OPENCL 1
  #endif

  #ifndef VIENNACL_WITH_EIGEN
    #define VIENNACL_WITH_EIGEN 1
  #endif

  #include <viennacl/ocl/backend.hpp>
  #include <viennacl/vector.hpp>
  #include <viennacl/matrix.hpp>
  #include <viennacl/linalg/direct_solve.hpp>
  #include <viennacl/linalg/lu.hpp>
  #include <viennacl/linalg/ilu_operations.hpp>

  #include <iostream>
  #include <string>
  #include <fstream>
    std::ifstream in("stan/math/prim/mat/fun/custom_kernels.cl",
    std::ios::in | std::ios::binary);
  std::string custom_kernels =
    std::string((std::istreambuf_iterator<char>(in)),
    std::istreambuf_iterator<char>());

#endif

#endif

