#ifndef STAN_MATH_OPENCL_STRINGIFY_HPP
#define STAN_MATH_OPENCL_STRINGIFY_HPP

// Used for importing the OpenCL kernels at compile time.
// There has been much discussion about the best ways to do this:
// https://github.com/bstatcomp/math/pull/7
// and https://github.com/stan-dev/math/pull/966
#ifndef STRINGIFY
#define STRINGIFY(...) #__VA_ARGS__
#endif

#endif  // STAN_MATH_OPENCL_STRINGIFY_HPP
