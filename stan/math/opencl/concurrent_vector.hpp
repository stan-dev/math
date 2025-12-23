#ifndef STAN_MATH_OPENCL_CONCURRENT_VECTOR_HPP
#define STAN_MATH_OPENCL_CONCURRENT_VECTOR_HPP

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>
#include <stdexcept>
#include <cassert>

namespace stan {
namespace math {
namespace internal {

/**
 * Minimal segmented concurrent_vector.
 *
 * Key properties:
 *  - concurrent emplace_back/push_back using an atomic size counter.
 *  - segmented storage => no moving elements during growth, stable addresses.
 *  - segments allocated lazily; allocation uses CAS to avoid locks.
 *
 * Important constraints / notes:
 *  - operator[] is safe if you only read indices < size() that are known to be
 * constructed.
 *  - For "publish-then-read" correctness: size() is updated before construction
 * finishes. So consumers must not read index i just because i < size(); they
 * must have a stronger protocol (e.g., producer hands out index, or you add a
 * "constructed" bitmap). This matches common usage where only the pushing
 * thread uses the returned index.
 *  - clear()/destruction are NOT concurrent with pushes.
 *
 * If you need "readers can iterate up to size() safely while writers push",
 * add a constructed flag per element (see comment near emplace_back()).
 */
template <typename T, std::size_t BaseSegmentSize = 1024,
          std::size_t MaxSegments = 32>
class concurrent_vector {
  static_assert(BaseSegmentSize > 0, "BaseSegmentSize must be > 0");
  static_assert((BaseSegmentSize & (BaseSegmentSize - 1)) == 0,
                "BaseSegmentSize must be a power of two (helps mapping).");

 public:
  concurrent_vector() : size_(0) {
    segments_.resize(MaxSegments);
    for (auto& p : segments_)
      p.store(nullptr, std::memory_order_relaxed);
  }

  concurrent_vector(const concurrent_vector&) = delete;
  concurrent_vector& operator=(const concurrent_vector&) = delete;

  ~concurrent_vector() { destroy_all_(); }

  std::size_t size() const noexcept {
    return size_.load(std::memory_order_acquire);
  }

  bool empty() const noexcept { return size() == 0; }

  // Non-concurrent: safe only when no other threads are pushing/reading.
  void clear() {
    destroy_all_();
    size_.store(0, std::memory_order_release);
  }

  // Concurrent append (construct in place). Returns the index.
  template <typename... Args>
  std::size_t emplace_back(Args&&... args) {
    // Claim an index
    const std::size_t idx = size_.fetch_add(1, std::memory_order_acq_rel);

    // Ensure the segment exists
    T* seg = ensure_segment_for_index_(idx);

    // Placement-new into the correct slot
    const std::size_t off = offset_in_segment_(idx);
    T* slot = seg + off;

    // Construct element
    ::new (static_cast<void*>(slot)) T(std::forward<Args>(args)...);

    // If you need "safe iteration by other threads that use size()",
    // you must publish construction completion separately, e.g.:
    // constructed_[idx].store(true, release);
    // and readers check constructed_[i].load(acquire).
    return idx;
  }

  std::size_t push_back(const T& v) { return emplace_back(v); }
  std::size_t push_back(T&& v) { return emplace_back(std::move(v)); }

  // Returns pointer to element at i (no bounds check).
  // Safe if element i is fully constructed and lifetime is valid.
  T* data_at(std::size_t i) noexcept {
    T* seg = segment_ptr_(segment_index_(i));
    return seg + offset_in_segment_(i);
  }
  const T* data_at(std::size_t i) const noexcept {
    const T* seg = segment_ptr_(segment_index_(i));
    return seg + offset_in_segment_(i);
  }

  // Bounds-checked access (still not concurrent-safe unless you have a
  // protocol).
  T& at(std::size_t i) {
    if (i >= size())
      throw std::out_of_range("concurrent_vector::at");
    return *data_at(i);
  }
  const T& at(std::size_t i) const {
    if (i >= size())
      throw std::out_of_range("concurrent_vector::at");
    return *data_at(i);
  }

  // Unchecked access
  T& operator[](std::size_t i) noexcept { return *data_at(i); }
  const T& operator[](std::size_t i) const noexcept { return *data_at(i); }

  // Capacity is segmented and unbounded until MaxSegments is exceeded.
  // This is the max number of elements representable by the segment scheme.
  static constexpr std::size_t max_size() noexcept {
    // Total capacity = Base * (2^MaxSegments - 1)
    // but beware overflow for large MaxSegments.
    return BaseSegmentSize * ((std::size_t{1} << MaxSegments) - 1);
  }

 private:
  // Segment k has size BaseSegmentSize * 2^k
  static constexpr std::size_t segment_size_(std::size_t k) noexcept {
    return BaseSegmentSize << k;
  }

  // Prefix count before segment k:
  // Base * (2^k - 1)
  static constexpr std::size_t segment_prefix_(std::size_t k) noexcept {
    return BaseSegmentSize * ((std::size_t{1} << k) - 1);
  }

  // Map global index -> segment index.
  // Let q = idx / Base. Then segment = floor(log2(q + 1)).
  static std::size_t segment_index_(std::size_t idx) noexcept {
    const std::size_t q = idx / BaseSegmentSize;
    const std::size_t x = q + 1;

#if defined(__GNUG__) || defined(__clang__)
    // floor(log2(x)) via clz
    return (sizeof(std::size_t) * 8 - 1)
           - static_cast<std::size_t>(__builtin_clzl(x));
#else
    // portable fallback
    std::size_t s = 0;
    std::size_t t = x;
    while (t >>= 1)
      ++s;
    return s;
#endif
  }

  static std::size_t offset_in_segment_(std::size_t idx) noexcept {
    const std::size_t s = segment_index_(idx);
    return idx - segment_prefix_(s);
  }

  T* segment_ptr_(std::size_t s) noexcept {
    return static_cast<T*>(segments_[s].load(std::memory_order_acquire));
  }
  const T* segment_ptr_(std::size_t s) const noexcept {
    return static_cast<const T*>(segments_[s].load(std::memory_order_acquire));
  }

  T* ensure_segment_for_index_(std::size_t idx) {
    const std::size_t s = segment_index_(idx);
    if (s >= MaxSegments) {
      throw std::length_error("concurrent_vector: exceeded MaxSegments");
    }

    T* seg = segment_ptr_(s);
    if (seg)
      return seg;

    // Allocate segment lazily (raw storage for T objects)
    const std::size_t n = segment_size_(s);
    void* raw = ::operator new(sizeof(T) * n);
    T* fresh = static_cast<T*>(raw);

    // CAS install; if another thread won, free ours.
    void* expected = nullptr;
    if (!segments_[s].compare_exchange_strong(expected, fresh,
                                              std::memory_order_release,
                                              std::memory_order_acquire)) {
      ::operator delete(raw);
      seg = static_cast<T*>(segments_[s].load(std::memory_order_acquire));
      assert(seg != nullptr);
      return seg;
    }

    return fresh;
  }

  // Destroy constructed elements and free segments.
  // Not concurrent with pushes or reads.
  void destroy_all_() noexcept {
    const std::size_t n = size_.load(std::memory_order_acquire);

    // Destroy elements that were constructed.
    // NOTE: This assumes indices [0, n) are all constructed.
    // If you allow exceptions or partial construction, track constructed flags.
    for (std::size_t i = 0; i < n; ++i) {
      data_at(i)->~T();
    }

    // Free segments
    for (std::size_t s = 0; s < segments_.size(); ++s) {
      void* p = segments_[s].load(std::memory_order_acquire);
      if (p) {
        ::operator delete(p);
        segments_[s].store(nullptr, std::memory_order_relaxed);
      }
    }
  }

  std::atomic<std::size_t> size_;
  std::vector<std::atomic<void*>> segments_;
};
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
