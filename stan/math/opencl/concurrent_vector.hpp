#ifndef STAN_MATH_OPENCL_CONCURRENT_VECTOR_HPP
#define STAN_MATH_OPENCL_CONCURRENT_VECTOR_HPP

#include <atomic>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * Segmented concurrent_vector.
 *
 * Properties:
 *  - concurrent emplace_back/push_back via atomic size counter
 *  - segmented storage => no moving elements during growth, stable addresses
 *  - segments allocated lazily; allocation uses CAS to avoid locks
 *  - movable (required for return-by-value usage such as vec_concat)
 *
 * Notes:
 *  - Intended for append-then-read patterns.
 *  - size_ increments before construction finishes. If readers iterate up to size()
 *    concurrently with writers, you need a constructed/published protocol.
 *  - clear()/destruction are NOT concurrent with pushes/reads.
 */
template <typename T, std::size_t BaseSegmentSize = 1024,
          std::size_t MaxSegments = 32>
class concurrent_vector {
  static_assert(BaseSegmentSize > 0, "BaseSegmentSize must be > 0");
  static_assert((BaseSegmentSize & (BaseSegmentSize - 1)) == 0,
                "BaseSegmentSize must be a power of two");
  static_assert(MaxSegments > 0, "MaxSegments must be > 0");

 public:
  concurrent_vector() noexcept : size_(0) {
    for (auto& p : segments_) p.store(nullptr, std::memory_order_relaxed);
  }

  concurrent_vector(const concurrent_vector& other) : concurrent_vector() {
    copy_from_(other);
  }

  concurrent_vector& operator=(const concurrent_vector& other) {
    if (this != &other) {
      clear();
      copy_from_(other);
    }
    return *this;
  }
  
  // Movable (needed so Stan can return-by-value)
  concurrent_vector(concurrent_vector&& other) noexcept : size_(0) {
    for (auto& p : segments_) p.store(nullptr, std::memory_order_relaxed);
    move_from_(other);
  }

  concurrent_vector& operator=(concurrent_vector&& other) noexcept {
    if (this != &other) {
      destroy_all_();
      for (auto& p : segments_) p.store(nullptr, std::memory_order_relaxed);
      size_.store(0, std::memory_order_relaxed);
      move_from_(other);
    }
    return *this;
  }

  ~concurrent_vector() noexcept { destroy_all_(); }

  std::size_t size() const noexcept {
    return size_.load(std::memory_order_acquire);
  }

  bool empty() const noexcept { return size() == 0; }

  // Not concurrent with pushes/reads.
  void clear() {
    destroy_all_();
    size_.store(0, std::memory_order_release);
  }

  // Pre-allocate enough segments to back indices [0, capacity-1].
  // Safe to call concurrently with emplace_back (may race allocating segments; losers free).
  void reserve(std::size_t capacity) {
    if (capacity == 0) return;
    const std::size_t last = capacity - 1;
    const std::size_t last_seg = segment_index_(last);
    if (last_seg >= MaxSegments) {
      throw std::length_error("concurrent_vector::reserve: exceeds MaxSegments");
    }
    for (std::size_t s = 0; s <= last_seg; ++s) {
      ensure_segment_(s);
    }
  }

  template <typename... Args>
  std::size_t emplace_back(Args&&... args) {
    const std::size_t idx = size_.fetch_add(1, std::memory_order_acq_rel);
    T* seg = ensure_segment_for_index_(idx);
    T* slot = seg + offset_in_segment_(idx);
    ::new (static_cast<void*>(slot)) T(std::forward<Args>(args)...);
    return idx;
  }

  std::size_t push_back(const T& v) { return emplace_back(v); }
  std::size_t push_back(T&& v) { return emplace_back(std::move(v)); }

  // Pointer helper (no bounds check).
  T* data_at(std::size_t i) noexcept {
    T* seg = segment_ptr_(segment_index_(i));
    return seg + offset_in_segment_(i);
  }
  const T* data_at(std::size_t i) const noexcept {
    const T* seg = segment_ptr_(segment_index_(i));
    return seg + offset_in_segment_(i);
  }

  // Bounds-checked access.
  T& at(std::size_t i) {
    if (i >= size()) throw std::out_of_range("concurrent_vector::at");
    return *data_at(i);
  }
  const T& at(std::size_t i) const {
    if (i >= size()) throw std::out_of_range("concurrent_vector::at");
    return *data_at(i);
  }

  // Unchecked access.
  T& operator[](std::size_t i) noexcept { return *data_at(i); }
  const T& operator[](std::size_t i) const noexcept { return *data_at(i); }

  // -------------------------
  // Iterators (InputIterator is enough for std::vector(first,last))
  // -------------------------

  class iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

    iterator() : v_(nullptr), i_(0) {}
    iterator(concurrent_vector* v, std::size_t i) : v_(v), i_(i) {}

    reference operator*() const { return (*v_)[i_]; }
    pointer operator->() const { return &(*v_)[i_]; }

    iterator& operator++() { ++i_; return *this; }
    iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }

    friend bool operator==(const iterator& a, const iterator& b) {
      return a.v_ == b.v_ && a.i_ == b.i_;
    }
    friend bool operator!=(const iterator& a, const iterator& b) { return !(a == b); }

   private:
    concurrent_vector* v_;
    std::size_t i_;
  };

  class const_iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = const T&;

    const_iterator() : v_(nullptr), i_(0) {}
    const_iterator(const concurrent_vector* v, std::size_t i) : v_(v), i_(i) {}

    reference operator*() const { return (*v_)[i_]; }
    pointer operator->() const { return &(*v_)[i_]; }

    const_iterator& operator++() { ++i_; return *this; }
    const_iterator operator++(int) { const_iterator tmp = *this; ++(*this); return tmp; }

    friend bool operator==(const const_iterator& a, const const_iterator& b) {
      return a.v_ == b.v_ && a.i_ == b.i_;
    }
    friend bool operator!=(const const_iterator& a, const const_iterator& b) { return !(a == b); }

   private:
    const concurrent_vector* v_;
    std::size_t i_;
  };

  iterator begin() noexcept { return iterator(this, 0); }
  iterator end() noexcept { return iterator(this, size()); }  // snapshot at call time

  const_iterator begin() const noexcept { return const_iterator(this, 0); }
  const_iterator end() const noexcept { return const_iterator(this, size()); }

  const_iterator cbegin() const noexcept { return const_iterator(this, 0); }
  const_iterator cend() const noexcept { return const_iterator(this, size()); }

  T& back() {
    const std::size_t n = size();
    if (n == 0) throw std::out_of_range("concurrent_vector::back on empty");
    return (*this)[n - 1];
  }

  const T& back() const {
    const std::size_t n = size();
    if (n == 0) throw std::out_of_range("concurrent_vector::back on empty");
    return (*this)[n - 1];
  }
  // -------------------------

 private:
  // Segment k has size BaseSegmentSize * 2^k
  static constexpr std::size_t segment_size_(std::size_t k) noexcept {
    return BaseSegmentSize << k;
  }

  // Prefix elements before segment k: Base * (2^k - 1)
  static constexpr std::size_t segment_prefix_(std::size_t k) noexcept {
    return BaseSegmentSize * ((std::size_t{1} << k) - 1);
  }

  // Map global index -> segment index
  // Let q = idx / Base. Then segment = floor(log2(q + 1)).
  static std::size_t segment_index_(std::size_t idx) noexcept {
    const std::size_t q = idx / BaseSegmentSize;
    const std::size_t x = q + 1;

#if defined(__GNUG__) || defined(__clang__)
    if constexpr (sizeof(std::size_t) == 8) {
      return 63u - static_cast<std::size_t>(
                       __builtin_clzll(static_cast<unsigned long long>(x)));
    } else {
      return 31u - static_cast<std::size_t>(
                       __builtin_clzl(static_cast<unsigned long>(x)));
    }
#else
    std::size_t s = 0;
    std::size_t t = x;
    while (t >>= 1) ++s;
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

  T* ensure_segment_(std::size_t s) {
    T* seg = segment_ptr_(s);
    if (seg) return seg;

    const std::size_t n = segment_size_(s);
    void* raw = ::operator new(sizeof(T) * n);
    T* fresh = static_cast<T*>(raw);

    void* expected = nullptr;
    if (!segments_[s].compare_exchange_strong(expected, fresh,
                                              std::memory_order_release,
                                              std::memory_order_acquire)) {
      ::operator delete(raw);
      seg = segment_ptr_(s);
      assert(seg != nullptr);
      return seg;
    }
    return fresh;
  }

  T* ensure_segment_for_index_(std::size_t idx) {
    const std::size_t s = segment_index_(idx);
    if (s >= MaxSegments) {
      throw std::length_error("concurrent_vector: exceeded MaxSegments");
    }
    return ensure_segment_(s);
  }

  void destroy_all_() noexcept {
    const std::size_t n = size_.load(std::memory_order_acquire);

    // Assumes [0, n) constructed.
    for (std::size_t i = 0; i < n; ++i) {
      data_at(i)->~T();
    }

    for (auto& a : segments_) {
      void* p = a.exchange(nullptr, std::memory_order_acq_rel);
      if (p) ::operator delete(p);
    }
  }

  void move_from_(concurrent_vector& other) noexcept {
    // Steal size
    const std::size_t n = other.size_.exchange(0, std::memory_order_acq_rel);
    size_.store(n, std::memory_order_release);

    // Steal segments
    for (std::size_t s = 0; s < MaxSegments; ++s) {
      void* p = other.segments_[s].exchange(nullptr, std::memory_order_acq_rel);
      segments_[s].store(p, std::memory_order_release);
    }
  }

  void copy_from_(const concurrent_vector& other) {
    const std::size_t n = other.size();
    if (n == 0) return;

    reserve(n);
    // Important: we want size_ to match, but we must construct elements.
    // Use emplace_back so construction happens in this container.
    for (std::size_t i = 0; i < n; ++i) {
      emplace_back(other[i]);
    }
  }
  
  std::atomic<std::size_t> size_;
  std::array<std::atomic<void*>, MaxSegments> segments_;
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
