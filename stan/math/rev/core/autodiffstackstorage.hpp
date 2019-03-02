#ifndef STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP
#define STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP

#include <stan/math/memory/stack_alloc.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Provides a thread_local singleton if needed. Read warnings below!
 * For performance reasons the singleton is a global static pointer
 * which is stored as thread local (TLS) if STAN_THREADS is
 * defined. The use of a pointer is motivated by performance reasons
 * for the threading case. When a TLS is used then initialization with
 * a constant expression is required for fast access to the TLS. As
 * the AD storage struct is non-POD it must be initialized as a
 * dynamic expression such that compilers will wrap any access to the
 * TLS by a TLS wrapper function which causes a significant
 * slow-down. A pointer to the AD storage instance can be initialized
 * to a compile-time constant expression of nullptr. In this case the
 * compiler avoids the use of a TLS wrapper function. Furthermore we
 * use the __thread keyword on compilers which support it. The
 * __thread keyword is a compiler-specific (gcc, clang, Intel)
 * extension which requires initialization with a compile time
 * constant expression. The C++11 keyword thread_local does allow for
 * constant and dynamic initialization of the TLS. Thus, only the
 * __thread keyword gurantees that constant initialization and it's
 * implied speedup, is used.
 *
 * The initialzation of the AD instance is handled by the init()
 * method which gets called by the constructor of the
 * singleton. Within stan-math the initialization of the AD instance
 * for the main thread of the program is handled by instantiating once
 * the singleton once in the init_chainablestack.hpp file. Whenever
 * STAN_THREADS is defined then all created child threads must call the
 * init() method of the AD singleton in order to initialize the TLS if
 * child threads want to perform AD operations.
 *
 * The design of a globally held (optionally TLS) pointer, which is
 * globally initialized, allows the compiler to apply necessary
 * inlining to get maximal performance. However, the design suffers
 * from "the static init order fiasco"[0]. Whenever the static init
 * order fiasco occurs, the C++ client of the library may call the
 * init method as needed to ensure proper initialization order. In
 * exchange, we get a more performant singleton pattern with automatic
 * initialization of the AD stack for the main thread. There has been
 * some discussion on earlier designs using the Mayer singleton
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
  typedef AutodiffStackSingleton<ChainableT, ChainableAllocT>
      AutodiffStackSingleton_t;

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

  static constexpr inline AutodiffStackStorage &instance() {
    return *instance_;
  }

 private:
  static bool init() {
    if (!instance_) {
      instance_ = new AutodiffStackStorage();
      return true;
    }
    return false;
  }

  static
#ifdef STAN_THREADS
#ifdef __GNUC__
      __thread
#else
      thread_local
#endif
#endif
      AutodiffStackStorage *instance_;

  bool own_instance_;
};

template <typename ChainableT, typename ChainableAllocT>
#ifdef STAN_THREADS
#ifdef __GNUC__
__thread
#else
thread_local
#endif
#endif
    typename AutodiffStackSingleton<ChainableT,
                                    ChainableAllocT>::AutodiffStackStorage
        *AutodiffStackSingleton<ChainableT, ChainableAllocT>::instance_
    = nullptr;

}  // namespace math
}  // namespace stan
#endif
