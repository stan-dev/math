# TBB Threading {#tbb_threading}

By default the `stan-math` library is not thread safe which is due to the requirement that the autodiff stack uses a global autodiff tape which records all operations of functions being evaluated. Starting with version 2.18 of `stan-math` threading support can be switched on using the compile-time switch `STAN_THREADS`. Defining this variable at compile time can be achieved by adding to `make/local`
```
CXXFLAGS += -DSTAN_THREADS
```
Once this is set `stan-math` will use the C++11 `thread_local` facility such that the autodiff stack is maintained per thread and not anymore globally. This allows the use of autodiff in a threaded application as this enables the calculation of the derivatives of a function inside a thread (the function itself may not use threads). Only if the function to be evaluated is composed of independent tasks it maybe possible to evaluate derivatives of a function in a threaded approach. An example of a threaded gradient evaluation is the `map_rect` function in `stan-math`.

In addition to making `stan-math` thread safe this also turns on parallel execution support of the `map_rect` function. Currently, the maximal number of threads being used by the function is controlled by the environment variable `STAN_NUM_THREADS` at runtime. Setting this variable to a positive integer number defines the maximal number of threads being used. In case the variable is set to the special value of `-1` this requests that as many threads as physical cores are being used. If the variable is not set a single thread is used. Any illegal value (not an integer, zero, other negative) will cause an exception to be thrown.

# Intel Threading Building Blocks

The Intel TBB library is used in stan-math since version 2.21.0. The Intel TBB library uses a threadpool internally and distributes work through a task-based approach. The tasks are dispatched to the threadpool via a the Intel TBB work-stealing scheduler. For example, whenever threading is enabled via `STAN_THREADS` the `map_rect` function in stan-math will use the `tbb::parallel_for` of the TBB. This will execute the work chunks given to `map_rect` with scheduling and thus load-balance CPU core utilization.

By default stan-math builds only the main `tbb` library by defining the `makefile` variable

```
TBB_LIBRARIES=tbb
```

The Intel TBB provides in addition to the main library memory allocators which are specifically designed to speedup threaded programs. These speedups have so far only been observed on MacOS systems for Stan programs such that on MacOS the default is set to

```
TBB_LIBRARIES=tbb tbbmalloc tbbmalloc_proxy
```

Users may override the default choices by defining `TBB_LIBRARIES` in the `make/local` file manually. Please refer to the [pull request](https://github.com/stan-dev/math/pull/1376) which merged the Intel TBB for further details on the performance evaluations.

# Requirements

Threading support requires a fully C++11 compliant compiler which has a working `thread_local` implementation. Below you find for each operating system what is known to work. Known to work configurations refers to run successfully by developers.

The compiler support for the C++11 `thread_local` keyword for the major open-source compilers is available since these versions:

- GNU g++ 4.8.1, [see here](https://gcc.gnu.org/projects/cxx-status.html); please also add `-pthread` to the `CXXFLAGS` variable
- clang++ 3.3, [see here](https://clang.llvm.org/cxx_status.html)  
  Note: clang + linux long had issues with `thread_local` which should be fixed with clang >=4.0

## Mac OS X

Known to work:
- macOS R toolchain, clang 4, [see here](https://github.com/coatless/r-macos-rtools)
- Apple's clang 9.1.0 (Xcode 9.1) on macOS High Sierra
- Apple's clang 11.0.0 (Xcode 9.1) on macOS Mojave
- g++ 6.4.0 from macports

Should work:
- Apple's clang compilers support the `thread_local` keyword since Mac OS Sierra (Xcode 8.0)

## Linux

Known to work:
- GNU g++ 4.9 - this configuration is also tested

With `clang` on linux there are issues during the linking step of programs which happens on old distributions like ubuntu trusty 14.04 (the 16.04 LTS is fine). A solution can be [found here.](https://stackoverflow.com/questions/29322666/undefined-reference-to-cxa-thread-atexitcxxabi-when-compiling-with-libc#30437761) It is likely that clang 4.0 *if used with libc++ 4.0* will work just fine, but developers have not yet confirmed this.

## Windows

Known to work:
- RTools 3.5 for Windows which uses mingw g++ 4.9.1 (since stan-math 2.20.0)
- RTools 4.0 for Windows with a port of GNU g++ 8.2 see [here](https://github.com/stan-dev/rstan/wiki/Using-RStan-with-the-R-3.6.0-Prerelease-on-Windows) but that compiler can also be used with CmdStan
