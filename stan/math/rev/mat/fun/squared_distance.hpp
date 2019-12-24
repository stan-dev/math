#ifndef STAN_MATH_REV_MAT_FUN_SQUARED_DISTANCE_HPP
#define STAN_MATH_REV_MAT_FUN_SQUARED_DISTANCE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/err/check_vector.hpp>
#include <stan/math/prim/arr/err/check_matching_sizes.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

namespace internal {

class squared_distance_vv_vari : public vari {
 protected:
  vari** v1_;
  vari** v2_;
  size_t length_;

  template <int R1, int C1, int R2, int C2>
  inline static double var_squared_distance(
      const Eigen::Matrix<var, R1, C1>& v1,
      const Eigen::Matrix<var, R2, C2>& v2) {
    using idx_t = typename index_type<matrix_v>::type;

    return (Eigen::Ref<const vector_v>(v1).val()
            - Eigen::Ref<const vector_v>(v2).val())
        .squaredNorm();
  }

 public:
  template <int R1, int C1, int R2, int C2>
  squared_distance_vv_vari(const Eigen::Matrix<var, R1, C1>& v1,
                           const Eigen::Matrix<var, R2, C2>& v2)
      : vari(var_squared_distance(v1, v2)), length_(v1.size()) {
    v1_ = reinterpret_cast<vari**>(
        ChainableStack::instance_->memalloc_.alloc(length_ * sizeof(vari*)));
    v2_ = reinterpret_cast<vari**>(
        ChainableStack::instance_->memalloc_.alloc(length_ * sizeof(vari*)));
    Eigen::Map<vector_vi>(v1_, length_) = v1.vi();
    Eigen::Map<vector_vi>(v2_, length_) = v2.vi();
  }
  virtual void chain() {
    Eigen::Map<vector_vi> v1_map(v1_, length_);
    Eigen::Map<vector_vi> v2_map(v2_, length_);
    vector_d di = 2 * adj_ * (v1_map.val() - v2_map.val());
    v1_map.adj() += di;
    v2_map.adj() -= di;
  }
};
class squared_distance_vd_vari : public vari {
 protected:
  vari** v1_;
  double* v2_;
  size_t length_;

  template <int R1, int C1, int R2, int C2>
  inline static double var_squared_distance(
      const Eigen::Matrix<var, R1, C1>& v1,
      const Eigen::Matrix<double, R2, C2>& v2) {
    using idx_t = typename index_type<matrix_d>::type;

    return (Eigen::Ref<const vector_v>(v1).val()
            - Eigen::Ref<const vector_d>(v2))
        .squaredNorm();
  }

 public:
  template <int R1, int C1, int R2, int C2>
  squared_distance_vd_vari(const Eigen::Matrix<var, R1, C1>& v1,
                           const Eigen::Matrix<double, R2, C2>& v2)
      : vari(var_squared_distance(v1, v2)), length_(v1.size()) {
    v1_ = reinterpret_cast<vari**>(
        ChainableStack::instance_->memalloc_.alloc(length_ * sizeof(vari*)));
    v2_ = reinterpret_cast<double*>(
        ChainableStack::instance_->memalloc_.alloc(length_ * sizeof(double)));
    Eigen::Map<vector_vi>(v1_, length_) = v1.vi();
    Eigen::Map<vector_d>(v2_, length_) = v2;
  }
  virtual void chain() {
    Eigen::Map<vector_vi> v1_map(v1_, length_);
    v1_map.adj()
        += 2 * adj_ * (v1_map.val() - Eigen::Map<vector_d>(v2_, length_));
  }
};
}  // namespace internal

template <int R1, int C1, int R2, int C2>
inline var squared_distance(const Eigen::Matrix<var, R1, C1>& v1,
                            const Eigen::Matrix<var, R2, C2>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return var(new internal::squared_distance_vv_vari(v1, v2));
}
template <int R1, int C1, int R2, int C2>
inline var squared_distance(const Eigen::Matrix<var, R1, C1>& v1,
                            const Eigen::Matrix<double, R2, C2>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return var(new internal::squared_distance_vd_vari(v1, v2));
}
template <int R1, int C1, int R2, int C2>
inline var squared_distance(const Eigen::Matrix<double, R1, C1>& v1,
                            const Eigen::Matrix<var, R2, C2>& v2) {
  check_vector("squared_distance", "v1", v1);
  check_vector("squared_distance", "v2", v2);
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return var(new internal::squared_distance_vd_vari(v2, v1));
}

}  // namespace math
}  // namespace stan
#endif
