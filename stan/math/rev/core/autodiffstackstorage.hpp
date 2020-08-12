#ifndef STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP
#define STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP

#include <stan/math/memory/stack_alloc.hpp>
#include <vector>

namespace stan {
namespace math {

// Internal macro used to modify global pointer definition to the
// global AD instance.
#ifdef STAN_THREADS
// Whenever STAN_THREADS is set a TLS keyword is used. For reasons
// explained below we use the GNU compiler extension __thread if
// supported by the compiler while the generic thread_local C++11
// keyword is used otherwise.
#ifdef __GNUC__
#define STAN_THREADS_DEF __thread
#else
#define STAN_THREADS_DEF thread_local
#endif
#else
// In case STAN_THREADS is not set, then no modifier is needed.
#define STAN_THREADS_DEF
#endif

/**
 * This struct always provides access to the autodiff stack using
 * the singleton pattern. Read warnings below!
 *
 * The singleton <code>instance_</code> is a global static pointer,
 * which is thread local (TLS) if the STAN_THREADS preprocess variable
 * is defined.
 *
 * The use of a pointer is motivated by performance reasons for the
 * threading case. When a TLS is used, initialization with a constant
 * expression at compile time is required for fast access to the
 * TLS. As the autodiff storage struct is non-POD, its initialization
 * is a dynamic expression at compile time. These dynamic expressions
 * are wrapped, in the TLS case, by a TLS wrapper function which slows
 * down its access. Using a pointer instead allows to initialize at
 * compile time to <code>nullptr</code>, which is a compile time
 * constant. In this case, the compiler avoids the use of a TLS
 * wrapper function.
 *
 * For performance reasons we use the __thread keyword on compilers
 * which support it. The __thread keyword is a GNU compiler-specific
 * (gcc, clang, Intel) extension which requires initialization with a
 * compile time constant expression. The C++11 keyword thread_local
 * does allow for constant and dynamic initialization of the
 * TLS. Thus, only the __thread keyword guarantees that constant
 * initialization and its implied speedup, is used.
 *
 * The initialization of the AD instance at run-time is handled by the
 * lifetime of a AutodiffStackSingleton object. More specifically, the
 * first instance of the AutodiffStackSingleton object will initialize
 * the AD instance and take ownership (it is the only one instance
 * with the private member own_instance_ being true). Thus, whenever
 * the first instance of the AutodiffStackSingleton object gets
 * destructed, the AD tape will be destructed as well. Within
 * stan-math the initialization of the AD instance for the main thread
 * of the program is handled by instantiating the singleton once in
 * the init_chainablestack.hpp file. Whenever STAN_THREADS is defined
 * then all created child threads must instantiate a
 * AutodiffStackSingleton object within the child thread before
 * accessing the AD system in order to initialize the TLS AD tape
 * within the child thread.
 *
 * The design of a globally held (optionally TLS) pointer, which is
 * globally initialized, allows the compiler to apply necessary
 * inlining to get maximal performance. However, the design suffers
 * from "the static init order fiasco"[0]. Whenever the static init
 * order fiasco occurs, the C++ client of the library may instantiate
 * an AutodiffStackSingleton object at the adequate code position prior
 * to any AD tape access to ensure proper initialization order. In
 * exchange, we get a more performant singleton pattern with automatic
 * initialization of the AD stack for the main thread. There has been
 * some discussion on earlier designs using the Meyer singleton
 * approach; see [1] and [2] and the discussions those PRs link to as
 * well.
 *
 * [0] https://isocpp.org/wiki/faq/ctors#static-init-order
 * [1] https://github.com/stan-dev/math/pull/840
 * [2] https://github.com/stan-dev/math/pull/826
 * [3]
 * http://discourse.mc-stan.org/t/potentially-dropping-support-for-older-versions-of-apples-version-of-clang/3780/
 */
template <typename ChainableT, typename ChainableAllocT>
struct AutodiffStackSingleton {
  using AutodiffStackSingleton_t
      = AutodiffStackSingleton<ChainableT, ChainableAllocT>;

  AutodiffStackSingleton() : own_instance_(init()) {}
  ~AutodiffStackSingleton() {
    if (own_instance_) {
      delete instance_;
      instance_ = nullptr;
    }
  }

  struct AutodiffStackStorage {
    AutodiffStackStorage &operator=(const AutodiffStackStorage &) = delete;

    std::vector<ChainableT *> var_stack_;
    std::vector<ChainableT *> var_nochain_stack_;
    std::vector<ChainableAllocT *> var_alloc_stack_;
    stack_alloc memalloc_;

    // nested positions
    std::vector<size_t> nested_var_stack_sizes_;
    std::vector<size_t> nested_var_nochain_stack_sizes_;
    std::vector<size_t> nested_var_alloc_stack_starts_;
  };

  explicit AutodiffStackSingleton(AutodiffStackSingleton_t const &) = delete;
  AutodiffStackSingleton &operator=(const AutodiffStackSingleton_t &) = delete;

  static STAN_THREADS_DEF AutodiffStackStorage *instance_;

 private:
  static bool init() {
    static STAN_THREADS_DEF bool is_initialized = false;
    if (!is_initialized) {
      is_initialized = true;
      instance_ = new AutodiffStackStorage();
      return true;
    }
    if (!instance_) {
      is_initialized = true;
      instance_ = new AutodiffStackStorage();
      return true;
    }
    return false;
  }

  bool own_instance_;
};

template <typename ChainableT, typename ChainableAllocT>
STAN_THREADS_DEF
    typename AutodiffStackSingleton<ChainableT,
                                    ChainableAllocT>::AutodiffStackStorage
        *AutodiffStackSingleton<ChainableT, ChainableAllocT>::instance_;

}  // namespace math
}  // namespace stan
#endif
