# OpenCL (`opencl`)

- All code is inside `#ifdef STAN_OPENCL`. To build it, set
  `STAN_OPENCL=true` in `make/local`.
- Data lives in `matrix_cl<T>`. Write new operations with the kernel
  generator (`opencl/kernel_generator/`) before hand-writing a kernel in
  `opencl/kernels/`.
- `opencl/prim` and `opencl/rev` mirror `stan/math/prim` and
  `stan/math/rev`; reuse the matching CPU function's checks and math.
- Guide: `doxygen/contributor_help_pages/add_new_opencl_kernel.md`.
- Tests need an OpenCL device. If you could not run them, say so.
