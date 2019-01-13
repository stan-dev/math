#ifndef STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP
#define STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP

#include <stan/math/memory/stack_alloc.hpp>
#include <vector>
#include <mutex>

#include "tbb/enumerable_thread_specific.h"
#include "tbb/task_arena.h"

namespace stan {
namespace math {

/**
 * Provides a thread_local singleton if needed. Read warnings below!
 * For performance reasons the singleton is a global static for the
 * case of no threading which is returned by a function. This design
 * should allow the compiler to apply necessary inlining to get
 * maximal performance. However, this design suffers from "the static
 * init order fiasco"[0].  Anywhere this is used, we must be
 * absolutely positive that it doesn't matter when the singleton will
 * get initialized relative to other static variables.  In exchange,
 * we get a more performant singleton pattern for the non-threading
 * case. In the threading case we use the defacto standard C++11
 * singleton pattern relying on a function wrapping a static local
 * variable. This standard pattern is expected to be well supported
 * by the major compilers (as its standard), but it does incur some
 * performance penalty.  There has been some discussion on this; see
 * [1] and [2] and the discussions those PRs link to as well.
 *
 * These are thread_local only if the user asks for it with
 * -DSTAN_THREADS. This is primarily because Apple clang compilers
 * before 2016 don't support thread_local and the additional
 * performance cost. We have proposed removing support for those[3],
 * and at that time we should evaluate the performance of a switch to
 * thread_local.  If there is no loss in performance, we can remove
 * this ifdef.
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

  static inline std::size_t get_new_id() {
    static std::mutex id_mutex;
    static std::size_t id = 0;
    std::lock_guard<std::mutex> id_lock(id_mutex);
    std::size_t new_id = ++id;
    return new_id;
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

    const std::size_t id_ = get_new_id();
  };

  AutodiffStackSingleton() = delete;
  explicit AutodiffStackSingleton(AutodiffStackSingleton_t const &) = delete;
  AutodiffStackSingleton &operator=(const AutodiffStackSingleton_t &) = delete;

  static std::vector<AutodiffStackStorage> thread_tapes_;
  //static int num_tapes_;

  static inline AutodiffStackStorage &instance() {
    // TBB TLS
    //return instance_.local();
    return thread_tapes_[tbb::this_task_arena::current_thread_index()];
  }

  // TBB TLS
  /*
  typedef tbb::enumerable_thread_specific<
      AutodiffStackStorage, tbb::cache_aligned_allocator<AutodiffStackStorage>,
      tbb::ets_key_per_instance>
      AutodiffStackStorage_tls_t;
  static AutodiffStackStorage_tls_t instance_;
  */
};

/*
template <typename ChainableT, typename ChainableAllocT>
typename AutodiffStackSingleton<ChainableT,
                                ChainableAllocT>::AutodiffStackStorage_tls_t
    AutodiffStackSingleton<ChainableT, ChainableAllocT>::instance_;
*/

template <typename ChainableT, typename ChainableAllocT>
std::vector<typename AutodiffStackSingleton<ChainableT,
                                   ChainableAllocT>::AutodiffStackStorage>
AutodiffStackSingleton<ChainableT, ChainableAllocT>::thread_tapes_(tbb::this_task_arena::max_concurrency());

/*
template <typename ChainableT, typename ChainableAllocT>
int AutodiffStackSingleton<ChainableT,
                                ChainableAllocT>::num_tapes_ = tbb::this_task_arena::max_concurrency();
*/

}  // namespace math
}  // namespace stan
#endif
