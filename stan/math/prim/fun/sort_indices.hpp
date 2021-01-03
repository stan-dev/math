#ifndef STAN_MATH_PRIM_FUN_SORT_INDICES_HPP
#define STAN_MATH_PRIM_FUN_SORT_INDICES_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

/**
 * A comparator that works for any container type that has the
 * brackets operator.
 *
 * @tparam ascending true if sorting in ascending order
 * @tparam C container type
 */
namespace internal {
template <bool ascending, typename C>
class index_comparator {
  const C& xs_;

 public:
  /**
   * Construct an index comparator holding a reference
   * to the specified container.
   *
   * @param xs Container
   */
  explicit index_comparator(const C& xs) : xs_(xs) {}

  /**
   * Return true if the value at the first index is sorted in
   * front of the value at the second index;  this will depend
   * on the template parameter <code>ascending</code>.
   *
   * @param i Index of first value for comparison
   * @param j Index of second value for comparison
   */
  bool operator()(int i, int j) const {
    if (ascending) {
      return xs_[i - 1] < xs_[j - 1];
    } else {
      return xs_[i - 1] > xs_[j - 1];
    }
  }
};

}  // namespace internal

/**
 * Return an integer array of indices of the specified container
 * sorting the values in ascending or descending order based on
 * the value of the first template parameter.
 *
 * @tparam ascending true if sort is in ascending order
 * @tparam C type of container
 * @param xs Container to sort
 * @return sorted version of container
 */
template <bool ascending, typename C>
std::vector<int> sort_indices(const C& xs) {
  using idx_t = index_type_t<C>;
  idx_t xs_size = xs.size();
  std::vector<int> idxs;
  idxs.resize(xs_size);
  for (idx_t i = 0; i < xs_size; ++i) {
    idxs[i] = i + 1;
  }
  internal::index_comparator<ascending, ref_type_t<C>> comparator(to_ref(xs));
  std::sort(idxs.begin(), idxs.end(), comparator);
  return idxs;
}

}  // namespace math
}  // namespace stan

#endif
