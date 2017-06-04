#ifdef STAN_GPU
  #ifndef VIENNACL_WITH_OPENCL
    #define VIENNACL_WITH_OPENCL 1
  #endif
#else
  #define STAN_GPU 0
#endif

#ifndef STAN_MATH_PRIM_MAT_FUN_VIENNACL_HPP
    #define STAN_MATH_PRIM_MAT_FUN_VIENNACL_HPP
    

    #ifndef VIENNACL_WITH_EIGEN
      #define VIENNACL_WITH_EIGEN 1
    #endif

    #include "viennacl/vector.hpp"
    #include "viennacl/matrix.hpp"
    #include "viennacl/linalg/direct_solve.hpp"
    #include "viennacl/linalg/lu.hpp"
    #include "viennacl/linalg/ilu_operations.hpp"

    #if STAN_GPU
        // a custom openCL kernel for copying the lower triangle to the upper
        static const char * custom_kernel_lower_upper =
        "__kernel void copy_lower_tri_upper_tri(\n"
        "          __global double * vec1,\n"
        "          unsigned int M,\n"
        "          unsigned int N, \n"
        "          unsigned int Npad) \n"
        "{ \n"
        "    int i = get_global_id(0); \n"
        "    int j = get_global_id(1); \n"
        "    if(i < M &&  j< N && i<j ){\n"
        "         vec1[i*Npad+j] = vec1[j*Npad+i];\n"
        "    }\n"
        "};\n";

        viennacl::ocl::program & my_prog = viennacl::ocl::current_context()
          .add_program(custom_kernel_lower_upper, "custom_kernel_lower_upper");
        viennacl::ocl::kernel & my_kernel_mul =
          my_prog.get_kernel("copy_lower_tri_upper_tri");
    #endif
#endif

