#ifndef STAN_MATH_PRIM_META_INDEX_HPP
#define STAN_MATH_PRIM_META_INDEX_HPP

#include <stan/math/prim/meta/is_vector.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Structure for an empty (size zero) index list.
 */
struct nil_index_list {};

/**
 * Template structure for an index list consisting of a head and
 * tail index.
 *
 * @tparam H type of index stored as the head of the list.
 * @tparam T type of index list stored as the tail of the list.
 */
template <typename H, typename T>
struct cons_index_list {
  std::decay_t<H> head_;
  std::decay_t<T> tail_;

  /**
   * Construct a non-empty index list with the specified index for
   * a head and specified index list for a tail.
   *
   * @param head Index for head.
   * @param tail Index list for tail.
   */
  template <typename Head, typename Tail>
  explicit constexpr cons_index_list(Head&& head, Tail&& tail) noexcept
      : head_(std::forward<Head>(head)), tail_(std::forward<Tail>(tail)) {}
};


/**
 * Structure for an indexing consisting of a single index.
 * Applying this index reduces the dimensionality of the container
 * to which it is applied by one.
 */
struct index_uni {
  int n_;

  /**
   * Construct a single indexing from the specified index.
   *
   * @param n single index.
   */
  explicit constexpr index_uni(int n) noexcept : n_(n) {}
};

// MULTIPLE INDEXING (does not reduce dimensionality)

/**
 * Structure for an indexing consisting of multiple indexes.  The
 * indexes do not need to be unique or in order.
 */
struct index_multi {
  std::vector<int> ns_;

  /**
   * Construct a multiple indexing from the specified indexes.
   *
   * @param ns multiple indexes.
   */
  template <typename T, require_std_vector_vt<std::is_integral, T>* = nullptr>
  explicit constexpr index_multi(T&& ns) noexcept : ns_(std::forward<T>(ns)) {}
};

/**
 * Structure for an indexing that consists of all indexes for a
 * container.  Applying this index is a no-op.
 */
struct index_omni {};

/**
 * Structure for an indexing from a minimum index (inclusive) to
 * the end of a container.
 */
struct index_min {
  int min_;

  /**
   * Construct an indexing from the specified minimum index (inclusive).
   *
   * @param min minimum index (inclusive).
   */
  explicit constexpr index_min(int min) noexcept : min_(min) {}
};

/**
 * Structure for an indexing from the start of a container to a
 * specified maximum index (inclusive).
 */
struct index_max {
  int max_;

  /**
   * Construct an indexing from the start of the container up to
   * the specified maximum index (inclusive).
   *
   * @param max maximum index (inclusive).
   */
  explicit constexpr index_max(int max) noexcept : max_(max) {}
};

/**
 * Structure for an indexing from a minimum index (inclusive) to a
 * maximum index (inclusive).
 */
struct index_min_max {
  int min_;
  int max_;
  /**
   * Return whether the index is positive or negative
   */
  bool is_positive_idx() const {
    return min_ <= max_;
  }
  /**
   * Construct an indexing from the specified minimum index
   * (inclusive) and maximum index (inclusive).
   *
   * @param min minimum index (inclusive).
   * @param max maximum index (inclusive).
   */
  constexpr index_min_max(int min, int max) noexcept
      : min_(min), max_(max) {}
};

}  // namespace model
}  // namespace stan
#endif
