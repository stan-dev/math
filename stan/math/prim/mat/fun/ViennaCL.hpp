#ifndef STAN_MATH_PRIM_MAT_FUN_VIENNACL_HPP
  #define STAN_MATH_PRIM_MAT_FUN_VIENNACL_HPP

  #ifdef STAN_GPU
    #ifndef VIENNACL_WITH_OPENCL
      #define VIENNACL_WITH_OPENCL 1
    #endif
  #endif

  #ifndef VIENNACL_WITH_EIGEN
    #define VIENNACL_WITH_EIGEN 1
  #endif

  #include <viennacl/vector.hpp>
  #include <viennacl/matrix.hpp>
  #include <viennacl/linalg/direct_solve.hpp>
  #include <viennacl/linalg/lu.hpp>
  #include <viennacl/linalg/ilu_operations.hpp>

  #ifdef STAN_GPU
    std::ifstream in("stan/math/prim/mat/fun/custom_kernels.cl", std::ios::in | std::ios::binary);
    std::string custom_kernels = std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(custom_kernels, "custom_kernels");
    viennacl::ocl::kernel & my_kernel_mul = my_prog.get_kernel("copy_lower_tri_upper_tri");
    viennacl::ocl::kernel & my_kernel_diag = my_prog.get_kernel("diagonal_mul");
    viennacl::ocl::kernel & my_kernel_inv1 = my_prog.get_kernel("inverse_step1");
    viennacl::ocl::kernel & my_kernel_inv2 = my_prog.get_kernel("inverse_step2");
    viennacl::ocl::kernel & my_kernel_inv3 = my_prog.get_kernel("inverse_step3");
  #endif

#endif

