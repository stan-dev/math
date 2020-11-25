#ifndef STAN_MATH_REV_FUN_SQUARED_DISTANCE_HPP
#define STAN_MATH_REV_FUN_SQUARED_DISTANCE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/squared_distance.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

class scal_squared_distance_vv_vari : public op_vv_vari {
 public:
  scal_squared_distance_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(squared_distance(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    double diff = avi_->val_ - bvi_->val_;
    avi_->adj_ += adj_ * 2.0 * diff;
    bvi_->adj_ -= adj_ * 2.0 * diff;
  }
};
class scal_squared_distance_vd_vari : public op_vd_vari {
 public:
  scal_squared_distance_vd_vari(vari* avi, double b)
      : op_vd_vari(squared_distance(avi->val_, b), avi, b) {}
  void chain() { avi_->adj_ += adj_ * 2 * (avi_->val_ - bd_); }
};

/**
 * Returns the squared distance.
 */
inline var squared_distance(const var& a, const var& b) {
  return {new scal_squared_distance_vv_vari(a.vi_, b.vi_)};
}

/**
 * Returns the squared distance.
 */
inline var squared_distance(const var& a, double b) {
  return {new scal_squared_distance_vd_vari(a.vi_, b)};
}

/**
 * Returns the squared distance.
 */
inline var squared_distance(double a, const var& b) {
  return {new scal_squared_distance_vd_vari(b.vi_, a)};
}

namespace internal {

class squared_distance_vv_vari : public vari {
 protected:
  vari** v1_;
  vari** v2_;
  size_t length_;

 public:
  template <
      typename EigVecVar1, typename EigVecVar2,
      require_all_eigen_vector_vt<is_var, EigVecVar1, EigVecVar2>* = nullptr>
  squared_distance_vv_vari(const EigVecVar1& v1, const EigVecVar2& v2)
      : vari((as_column_vector_or_scalar(v1).val()
              - as_column_vector_or_scalar(v2).val())
                 .squaredNorm()),
        length_(v1.size()) {
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

 public:
  template <typename EigVecVar, typename EigVecArith,
            require_eigen_vector_vt<is_var, EigVecVar>* = nullptr,
            require_eigen_vector_vt<std::is_arithmetic, EigVecArith>* = nullptr>
  squared_distance_vd_vari(const EigVecVar& v1, const EigVecArith& v2)
      : vari((as_column_vector_or_scalar(v1).val()
              - as_column_vector_or_scalar(v2))
                 .squaredNorm()),
        length_(v1.size()) {
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

template <
    typename EigVecVar1, typename EigVecVar2,
    require_all_eigen_vector_vt<is_var, EigVecVar1, EigVecVar2>* = nullptr>
inline var squared_distance(const EigVecVar1& v1, const EigVecVar2& v2) {
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return {new internal::squared_distance_vv_vari(to_ref(v1), to_ref(v2))};
}

template <typename EigVecVar, typename EigVecArith,
          require_eigen_vector_vt<is_var, EigVecVar>* = nullptr,
          require_eigen_vector_vt<std::is_arithmetic, EigVecArith>* = nullptr>
inline var squared_distance(const EigVecVar& v1, const EigVecArith& v2) {
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return {new internal::squared_distance_vd_vari(to_ref(v1), to_ref(v2))};
}

template <typename EigVecArith, typename EigVecVar,
          require_eigen_vector_vt<std::is_arithmetic, EigVecArith>* = nullptr,
          require_eigen_vector_vt<is_var, EigVecVar>* = nullptr>
inline var squared_distance(const EigVecArith& v1, const EigVecVar& v2) {
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return {new internal::squared_distance_vd_vari(to_ref(v2), to_ref(v1))};
}

/**
 * Compute the squared distance between the elements in
 * two inputs.
 *
 * This overload handles arguments where one of T1 or T2 are
 * `var_value<T>` where `T` is an Eigen type. The other type can
 * also be a `var_value` or it can be a matrix type that inherits
 * from EigenBase
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param A first argument
 * @param B second argument
 * @return sum of squared difference of A and B
 */
template <typename T1, typename T2,
	  require_all_matrix_t<T1, T2>* = nullptr,
          require_any_var_matrix_t<T1, T2>* = nullptr>
inline var squared_distance(const T1& A, const T2& B) {
  check_matching_sizes("squared_distance", "A", A.val(), "B", B.val());

  if (A.size() == 0)
    return 0.0;

  var res;

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<var, T2>> arena_B = B;

    double res_val = 0.0;
    for(size_t i = 0; i < arena_A.size(); ++i) {
      double diff = arena_A.val().coeff(i) - arena_B.val().coeff(i);
      res_val += diff * diff;
    }
    res = res_val;

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      double res_adj = res.adj();
      
      for(size_t i = 0; i < arena_A.size(); ++i) {
	double diff = arena_A.val().coeff(i) - arena_B.val().coeff(i);
	arena_A.adj().coeffRef(i) += 2.0 * res_adj * diff;
	arena_B.adj().coeffRef(i) -= 2.0 * res_adj * diff;
      }
    });
  } else if (!is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<double, T2>> arena_B = value_of(B);

    double res_val = 0.0;
    for(size_t i = 0; i < arena_A.size(); ++i) {
      double diff = arena_A.val().coeff(i) - arena_B.coeff(i);
      res_val += diff * diff;
    }
    res = res_val;

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      double res_adj = res.adj();
      
      for(size_t i = 0; i < arena_A.size(); ++i) {
	double diff = arena_A.val().coeff(i) - arena_B.coeff(i);
	arena_A.adj().coeffRef(i) += 2.0 * res_adj * diff;
      }
    });
  } else {
    arena_t<promote_scalar_t<double, T1>> arena_A = value_of(A);
    arena_t<promote_scalar_t<var, T2>> arena_B = B;

    double res_val = 0.0;
    for(size_t i = 0; i < arena_A.size(); ++i) {
      double diff = arena_A.coeff(i) - arena_B.val().coeff(i);
      res_val += diff * diff;
    }
    res = res_val;

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      double res_adj = res.adj();
      
      for(size_t i = 0; i < arena_A.size(); ++i) {
	double diff = arena_A.coeff(i) - arena_B.val().coeff(i);
	arena_B.adj().coeffRef(i) -= 2.0 * res_adj * diff;
      }
    });
  }

  return res;
}

}  // namespace math
}  // namespace stan
#endif
