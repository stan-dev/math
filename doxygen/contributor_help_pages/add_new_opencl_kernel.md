## OpenCL Developer Guide {#opencl_guide}

The OpenCL integration with Stan is meant to be seamless and user friendly - setting a flag moves the computation of supported routines to the GPU, with no need to change Stan code. The API provides experts with a simple way of implementing their GPU kernels or to use existing GPU kernels as building blocks. In the below we detail the main pieces of the architecture and how to add your own OpenCL routines.

Stan's OpenCL backend uses a single context to receive data and routines from individual nodes in the expression tree. Ensuring there is only one context and queue per device for the life of the program makes context management simpler. The implementation of the OpenCL context which manages the device queue follows the Meyers singleton pattern and sits in the class `opencl_context_base`.

Instead of calling `opencl_context_base::getInstance().method()`, developers can access the context through a friend adapter class `opencl_context` which provides an API for accessing the base context. If the OpenCL implementation supports asynchronous operations, then the context asynchronously executes kernels. Asynchronous operations are particularly useful in conjunction with threading as the individual threads will be able to enqueue operations, which will execute while threads do other calculations using CPU.

## Matrices On The GPU

The base matrix class `matrix_cl` holds the device memory buffer, meta-information on the matrix, and methods for reading and writing event stacks for asynchronous computation. When a kernel receives a `matrix_cl`, the kernel's event is attached to the appropriate read or write event stack. Reading and writing to OpenCL buffers uses the generic `enqueueWriteBuffer` and `enqueueWriteBuffer` methods. Because Stan Math heavily relies on Eigen matrices, constructors and methods are available for passing data back and forth.

Developers can pass in Eigen matrices directly to the `matrix_cl` constructor or use the `to_matrix_cl()` or `from_matrix_cl()` methods.

```cpp

Eigen::MatrixXd m(2, 2);
m << 1, 2, 3, 4;

matrix_cl<double> A(m);
matrix_cl<double> B(2, 2);

B = to_matrix_cl(m);
Eigen::MatrixXd C = from_matrix_cl(B);
```

Similar constructors for the `matrix_cl` class are available for standard vectors `std::vector<T>` and arrays of doubles.

We can reduce the data transfers of triangular matrices by only transferring the non-zero parts of the matrix in a packed form. The kernel `unpack` deals with unpacking the packed form shown on the right-hand side on to a flat matrix shown on the left-hand side. For lower (upper) triangular matrices, the upper (lower) triangular fill with zeros. The kernel `pack` packs the flat matrix to packed form for the transfer back to the host's global memory.

```cpp
  matrix_cl<double> L = packed_copy(L_val_cpu, M_, matrix_view::Lower);
```

## Adding New OpenCL Kernels

The OpenCL specification demands that strings are used to represent OpenCL kernels. However, having a large number of files comprised of strings is unwieldy and difficult to maintain. Stan wraps its kernels inside of a STRINGIFY macro, which gives developers access to the standard set of developer tools such as code highlighting, linting, Doxygen, and auto-formatting. This style makes the kernel code easier to maintain compared to having files full of strings.

```cpp
const char *example_kernel_code = STRINGIFY(

  __kernel void example(double *A, double *B, int *val) {
  	  // kernel code...
  }

);

const kernel_cl<out_buffer, in_buffer, int> example(
  "example", example_kernel_code, {"THREAD_BLOCK_SIZE", 32);
```


In the above, a developer uses `STRINGIFY` to create a `const char*` that holds the kernel code. That string passes into the `kernel_cl` struct templated by the kernel argument types and with arguments giving the name of the kernel, the kernel code, and optional kernel macros they would like to have defined in the kernel.

Internally, we keep track of OpenCL events via queues on each `matrix_cl` object that we use to conservatively prevent race conditions and provide ordering where necessary. `out_buffer` and `in_buffer` are empty structs that we pass as template arguments to configure the kernel during construction to indicate the directionality of each input buffer. At runtime, the kernel will check the correct event queues on its arguments for events it needs to wait for and then attach the event representing the kernel's completion to each buffer's queues correctly. That way we ensure that an operation that, for example, writes to a buffer, is completed before we allow the OpenCL runtime to read from that buffer.


When the `kernel_cl` struct is constructed, it compiles the kernel and developers call the kernel with

```cpp
matrix_cl<double> foo = //...
matrix_cl<double> goo;
example(cl::NDRange(...), goo, foo, 10);
```

Depending on the`in/out_buffer` passed when constructing the kernel, events will be added to the appropriate `matrix_cl` read and/or write event stack. For instance, in the above, `goo` in the output and will have the kernel's event attached to it's `write_stack`. While `foo` will have the kernel's event attached to it's `read_stack`. Later kernel calls that write to `foo` will know to wait for all the event's in `foo`'s `read_stack` and `write_stack` while kernels that use `goo` as input will know to wait for the event's in `goo`'s `write_stack`.
