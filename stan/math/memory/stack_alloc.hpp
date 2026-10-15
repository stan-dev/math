#ifndef STAN_MATH_MEMORY_STACK_ALLOC_HPP
#define STAN_MATH_MEMORY_STACK_ALLOC_HPP

// TODO(Bob): <cstddef> replaces this ifdef in C++11, until then this
//            is best we can do to get safe pointer casts to uints.
#include <stdint.h>
#include <stan/math/prim/meta.hpp>
#include <cstdlib>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <vector>
#ifdef _WIN32
#include <malloc.h>
#endif

/**
 * `STAN_ALIGNED_ALLOC(alignment, size)` allocates `size` bytes aligned to
 * `alignment` bytes and `STAN_ALIGNED_FREE(ptr)` releases that memory.
 *
 * MSVC and MinGW do not provide `std::aligned_alloc`, so on Windows these
 * use `_aligned_malloc` and `_aligned_free`. Memory from `_aligned_malloc`
 * must not be released with `free()`. Defining `STAN_ALIGNED_ALLOC` before
 * including Stan replaces both, so `STAN_ALIGNED_FREE` must be defined too.
 */
#ifndef STAN_ALIGNED_ALLOC
#ifdef _WIN32
#define STAN_ALIGNED_ALLOC(alignment, size) _aligned_malloc((size), (alignment))
#define STAN_ALIGNED_FREE(ptr) _aligned_free(ptr)
#else
#define STAN_ALIGNED_ALLOC(alignment, size) \
  std::aligned_alloc((alignment), (size))
#define STAN_ALIGNED_FREE(ptr) std::free(ptr)
#endif
#endif

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the specified pointer is aligned
 * on the number of bytes.
 *
 * This doesn't really make sense other than for powers of 2.
 *
 * @param ptr Pointer to test.
 * @param bytes_aligned Number of bytes of alignment required.
 * @return <code>true</code> if pointer is aligned.
 * @tparam Type of object to which pointer points.
 */
template <typename T>
inline bool is_aligned(T* ptr, unsigned int bytes_aligned) {
  return (reinterpret_cast<uintptr_t>(ptr) % bytes_aligned) == 0U;
}

namespace internal {
const size_t DEFAULT_INITIAL_NBYTES = 1 << 16;  // 64KB

/**
 * Allocate a block of memory aligned to `alignment` bytes with
 * `STAN_ALIGNED_ALLOC`. The block must be released with
 * `STAN_ALIGNED_FREE`.
 *
 * @param size Number of bytes to allocate; rounded up to a multiple of
 * `alignment` as `std::aligned_alloc` requires.
 * @param alignment Alignment in bytes, a power of two.
 * @return Pointer to the block, or `nullptr` if the allocation failed.
 */
inline char* aligned_malloc(size_t size, size_t alignment) {
  const size_t rounded_size = (size + alignment - 1) / alignment * alignment;
  return static_cast<char*>(STAN_ALIGNED_ALLOC(alignment, rounded_size));
}

/**
 * Return the first address at or after `ptr` that is aligned to
 * `alignment` bytes.
 *
 * @param ptr Pointer to align.
 * @param alignment Alignment in bytes, a power of two.
 */
inline char* align_up(char* ptr, size_t alignment) {
  const uintptr_t address = reinterpret_cast<uintptr_t>(ptr);
  const uintptr_t aligned
      = (address + alignment - 1) & ~static_cast<uintptr_t>(alignment - 1);
  return ptr + (aligned - address);
}
}  // namespace internal

/**
 * An instance of this class provides a memory pool through
 * which blocks of raw memory may be allocated and then collected
 * simultaneously.
 *
 * This class is useful in settings where large numbers of small
 * objects are allocated and then collected all at once.  This may
 * include objects whose destructors have no effect.
 *
 * Memory is allocated on a stack of blocks.  Each block allocated
 * is twice as large as the previous one.  The memory may be
 * recovered, with the blocks being reused, or all blocks may be
 * freed, resetting the stack of blocks to its original state.
 *
 * Every allocation is at least 8 byte aligned. Larger alignments can be
 * requested through the `alignment` argument of `alloc()`.
 */
class stack_alloc {
 private:
  std::vector<char*> blocks_;  // storage for blocks,
                               // may be bigger than cur_block_
  std::vector<size_t> sizes_;  // could store initial & shift for others
  size_t cur_block_;           // index into blocks_ for next alloc
  char* cur_block_end_;        // ptr to cur_block_ptr_ + sizes_[cur_block_]
  char* next_loc_;             // ptr to next available spot in cur
                               // block
  // next three for keeping track of nested allocations on top of stack:
  std::vector<size_t> nested_cur_blocks_;
  std::vector<char*> nested_next_locs_;
  std::vector<char*> nested_cur_block_ends_;

  /**
   * Return `true` if block `i` has room for `len` bytes starting at an
   * address aligned to `alignment` bytes.
   */
  inline bool block_fits(size_t i, size_t len, size_t alignment) const {
    const size_t offset
        = internal::align_up(blocks_[i], alignment) - blocks_[i];
    return offset + len <= sizes_[i];
  }

  /**
   * Moves us to the next block of memory, allocating that block
   * if necessary, and allocates len bytes of memory within that
   * block, aligned to `alignment` bytes.
   *
   * @param len Number of bytes to allocate.
   * @param alignment Alignment in bytes, a power of two. New blocks come
   * from `std::aligned_alloc` with this alignment.
   * @return A pointer to the allocated memory.
   */
  STAN_STRONG_INLINE inline char* move_to_next_block(size_t len,
                                                     size_t alignment = 8) {
    ++cur_block_;
    // Find the next block (if any) with room for len aligned bytes.
    while ((cur_block_ < blocks_.size())
           && !block_fits(cur_block_, len, alignment)) {
      ++cur_block_;
    }
    // Allocate a new block if necessary.
    if (unlikely(cur_block_ >= blocks_.size())) {
      // New block should be max(2*size of last block, len) bytes.
      size_t newsize = sizes_.back() * 2;
      if (newsize < len) {
        newsize = len;
      }
      blocks_.push_back(internal::aligned_malloc(newsize, alignment));
      if (!blocks_.back()) {
        throw std::bad_alloc();
      }
      sizes_.push_back(newsize);
    }
    char* result = internal::align_up(blocks_[cur_block_], alignment);
    // Get the object's state back in order.
    next_loc_ = result + len;
    cur_block_end_ = blocks_[cur_block_] + sizes_[cur_block_];
    return result;
  }

 public:
  /**
   * Construct a resizable stack allocator initially holding the
   * specified number of bytes.
   *
   * @param initial_nbytes Initial number of bytes for the
   * allocator.  Defaults to <code>(1 << 16) = 64KB</code> initial bytes.
   * @throws std::bad_alloc if the initial block cannot be allocated.
   */
  explicit stack_alloc(size_t initial_nbytes = internal::DEFAULT_INITIAL_NBYTES)
      : blocks_(1, internal::aligned_malloc(initial_nbytes, 8)),
        sizes_(1, initial_nbytes),
        cur_block_(0),
        cur_block_end_(blocks_[0] + initial_nbytes),
        next_loc_(blocks_[0]) {
    if (!blocks_[0]) {
      throw std::bad_alloc();  // no msg allowed in bad_alloc ctor
    }
  }

  /**
   * Destroy this memory allocator.
   *
   * This is implemented as a no-op as there is no destruction
   * required.
   */
  ~stack_alloc() {
    // free ALL blocks
    for (auto& block : blocks_) {
      if (block) {
        STAN_ALIGNED_FREE(block);
      }
    }
  }

  /**
   * Return a newly allocated block of memory of the appropriate size
   * managed by the stack allocator, aligned to `alignment` bytes.
   *
   * Every allocation is at least 8 byte aligned; an `alignment` below 8
   * gives 8 byte alignment. The reserved space is padded up to a multiple
   * of 8 bytes, so allocations after this one stay 8 byte aligned.
   *
   * This function may call C++'s <code>malloc()</code> function,
   * with any exceptions percolated through this function.
   *
   * @param[in] len Number of bytes to allocate.
   * @param[in] alignment Alignment in bytes, a power of two.
   * @return A pointer to the allocated memory.
   */
  STAN_STRONG_INLINE inline void* alloc(size_t len, size_t alignment = 8) {
    const size_t padded_len = (len + 7) & ~static_cast<size_t>(7);
    // Every block and every allocation is 8 byte aligned, so only
    // larger alignments need to round up the next location.
    char* result
        = alignment <= 8 ? next_loc_ : internal::align_up(next_loc_, alignment);
    // Typically, just return and increment the next location.
    next_loc_ = result + padded_len;
    // Occasionally, we have to switch blocks.
    if (unlikely(next_loc_ >= cur_block_end_)) {
      result = move_to_next_block(padded_len, alignment);
    }
    return result;
  }

  /**
   * Allocate an array on the arena of the specified size to hold
   * values of the specified template parameter type.
   *
   * @tparam T type of entries in allocated array.
   * @param[in] n size of array to allocate.
   * @param[in] alignment Alignment in bytes, a power of two.
   * @return new array allocated on the arena.
   */
  template <typename T>
  STAN_STRONG_INLINE inline T* alloc_array(size_t n, size_t alignment = 8) {
    return static_cast<T*>(alloc(n * sizeof(T), alignment));
  }

  /**
   * Recover all the memory used by the stack allocator.  The stack
   * of memory blocks allocated so far will be available for further
   * allocations.  To free memory back to the system, use the
   * function free_all().
   */
  inline void recover_all() {
    cur_block_ = 0;
    next_loc_ = blocks_[0];
    cur_block_end_ = next_loc_ + sizes_[0];
  }

  /**
   * Store current positions before doing nested operation so can
   * recover back to start.
   */
  inline void start_nested() {
    nested_cur_blocks_.push_back(cur_block_);
    nested_next_locs_.push_back(next_loc_);
    nested_cur_block_ends_.push_back(cur_block_end_);
  }

  /**
   * recover memory back to the last start_nested call.
   */
  inline void recover_nested() {
    if (unlikely(nested_cur_blocks_.empty())) {
      recover_all();
    }

    cur_block_ = nested_cur_blocks_.back();
    nested_cur_blocks_.pop_back();

    next_loc_ = nested_next_locs_.back();
    nested_next_locs_.pop_back();

    cur_block_end_ = nested_cur_block_ends_.back();
    nested_cur_block_ends_.pop_back();
  }

  /**
   * Free all memory used by the stack allocator other than the
   * initial block allocation back to the system.  Note:  the
   * destructor will free all memory.
   */
  inline void free_all() {
    // frees all BUT the first (index 0) block
    for (size_t i = 1; i < blocks_.size(); ++i) {
      if (blocks_[i]) {
        STAN_ALIGNED_FREE(blocks_[i]);
      }
    }
    sizes_.resize(1);
    blocks_.resize(1);
    recover_all();
  }

  /**
   * Return number of bytes allocated to this instance by the heap.
   * This is not the same as the number of bytes allocated through
   * calls to memalloc_.  The latter number is not calculable
   * because space is wasted at the end of blocks if the next
   * alloc request doesn't fit.  (Perhaps we could trim down to
   * what is actually used?)
   *
   * @return number of bytes allocated to this instance
   */
  inline size_t bytes_allocated() const {
    size_t sum = 0;
    for (size_t i = 0; i <= cur_block_; ++i) {
      sum += sizes_[i];
    }
    return sum;
  }

  /**
   * Return the approximate number of bytes handed out through calls to
   * alloc(). Blocks before the current one are counted in full, so any
   * space wasted at the end of a block when moving to the next one is
   * included.
   *
   * @return approximate number of bytes used on this instance
   */
  inline size_t approx_bytes_used() const {
    size_t sum = 0;
    for (size_t i = 0; i < cur_block_; ++i) {
      sum += sizes_[i];
    }
    return sum + static_cast<size_t>(next_loc_ - blocks_[cur_block_]);
  }

  /**
   * Indicates whether the memory in the pointer
   * is in the stack.
   *
   * @param[in] ptr memory location
   * @return true if the pointer is in the stack,
   *    false otherwise.
   */
  inline bool in_stack(const void* ptr) const {
    for (size_t i = 0; i < cur_block_; ++i) {
      if (ptr >= blocks_[i] && ptr < blocks_[i] + sizes_[i]) {
        return true;
      }
    }
    if (ptr >= blocks_[cur_block_] && ptr < next_loc_) {
      return true;
    }
    return false;
  }
};

}  // namespace math
}  // namespace stan
#endif
