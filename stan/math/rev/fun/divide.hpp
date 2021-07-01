#ifndef STAN_MATH_REV_FUN_DIVIDE_HPP
#define STAN_MATH_REV_FUN_DIVIDE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_var.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {
namespace internal {

template <int R, int C>
class matrix_scalar_divide_dv_vari : public vari {
 public:
  int rows_;
  int cols_;
  vari* adjCRef_;
  vari** adjResultRef_;
  double invc_;

  explicit matrix_scalar_divide_dv_vari(const Eigen::Matrix<double, R, C>& m,
                                        const var& c)
      : vari(0),
        rows_(m.rows()),
        cols_(m.cols()),
        adjCRef_(c.vi_),
        adjResultRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())),
        invc_(1.0 / c.val()) {
    Eigen::Map<matrix_vi>(adjResultRef_, rows_, cols_)
        = (invc_ * m).unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    Eigen::Map<matrix_vi> adjResult(adjResultRef_, rows_, cols_);
    adjCRef_->adj_
        -= invc_ * (adjResult.adj().array() * adjResult.val().array()).sum();
  }
};

template <int R, int C>
class matrix_scalar_divide_vd_vari : public vari {
 public:
  int rows_;
  int cols_;
  vari** adjMRef_;
  vari** adjResultRef_;
  double invc_;

  explicit matrix_scalar_divide_vd_vari(const Eigen::Matrix<var, R, C>& m,
                                        const double& c)
      : vari(0),
        rows_(m.rows()),
        cols_(m.cols()),
        adjMRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())),
        adjResultRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())),
        invc_(1.0 / c) {
    Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_) = m.vi();
    Eigen::Map<matrix_vi>(adjResultRef_, rows_, cols_)
        = (invc_ * m.val()).unaryExpr([](double x) {
            return new vari(x, false);
          });
  }

  virtual void chain() {
    Eigen::Map<matrix_vi> adjM(adjMRef_, rows_, cols_);
    Eigen::Map<matrix_vi> adjResult(adjResultRef_, rows_, cols_);
    adjM.adj() += invc_ * adjResult.adj();
  }
};

template <int R, int C>
class matrix_scalar_divide_vv_vari : public vari {
 public:
  int rows_;
  int cols_;
  vari** adjMRef_;
  vari* adjC_;
  vari** adjResultRef_;
  double invc_;

  explicit matrix_scalar_divide_vv_vari(const Eigen::Matrix<var, R, C>& m,
                                        const var& c)
      : vari(0),
        rows_(m.rows()),
        cols_(m.cols()),
        adjMRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())),
        adjC_(c.vi_),
        adjResultRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())),
        invc_(1.0 / c.val()) {
    Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_) = m.vi();
    Eigen::Map<matrix_vi>(adjResultRef_, rows_, cols_)
        = (invc_ * m.val()).unaryExpr([](double x) {
            return new vari(x, false);
          });
  }

  virtual void chain() {
    Eigen::Map<matrix_vi> adjM(adjMRef_, rows_, cols_);
    Eigen::Map<matrix_vi> adjResult(adjResultRef_, rows_, cols_);
    adjC_->adj_
        -= invc_ * (adjResult.adj().array() * adjResult.val().array()).sum();
    adjM.adj() += invc_ * adjResult.adj();
  }
};

}  // namespace internal

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat type of the matrix or expression
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, require_eigen_vt<std::is_arithmetic, Mat>* = nullptr>
inline auto divide(const Mat& m, const var& c) {
  auto inv_c = (1.0 / c.val());
  arena_t<promote_scalar_t<var, Mat>> res = inv_c * m.array();
  reverse_pass_callback([c, inv_c, res]() mutable {
    c.adj() -= inv_c * (res.adj().array() * res.val().array()).sum();
  });
  return promote_scalar_t<var, Mat>(res);
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat type of the matrix or expression
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, require_matrix_st<is_var, Mat>* = nullptr>
inline auto divide(const Mat& m, const double& c) {
  arena_t<promote_scalar_t<var, Mat>> arena_m = m;
  auto inv_c = (1.0 / c);
  arena_t<promote_scalar_t<var, Mat>> res = inv_c * arena_m.val();
  reverse_pass_callback([c, inv_c, arena_m, res]() mutable {
    arena_m.adj() += inv_c * res.adj_op();
  });
  return promote_scalar_t<var, Mat>(res);
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat type of the matrix or expression
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, require_matrix_st<is_var, Mat>* = nullptr>
inline auto divide(const Mat& m, const var& c) {
  arena_t<plain_type_t<Mat>> arena_m = m;
  auto inv_c = (1.0 / c.val());
  arena_t<plain_type_t<Mat>> res = inv_c * arena_m.val();
  reverse_pass_callback([c, inv_c, arena_m, res]() mutable {
    c.adj() -= inv_c * (res.adj().array() * res.val().array()).sum();
    arena_m.adj() += inv_c * res.adj();
  });
  return plain_type_t<Mat>(res);
}

template <typename Mat1, typename Mat2, require_any_matrix_st<is_var, Mat1, Mat2>* = nullptr>
inline auto divide(const Mat1& m1, const Mat2& m2) {
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    auto inv_m2 = to_arena(arena_m2.val().array().inverse());
    using val_ret = decltype((inv_m2 * arena_m1.val().array()).matrix().eval());
    using ret_type = return_var_matrix_t<val_ret, Mat1, Mat2>;
    arena_t<ret_type> res = (inv_m2.array() * arena_m1.val().array()).matrix();
    reverse_pass_callback([inv_m2, arena_m1, arena_m2, res]() mutable {
      auto inv_times_res = (inv_m2 * res.adj().array()).eval();
      arena_m1.adj().array() += inv_times_res;
      arena_m2.adj().array() -= inv_times_res * res.val().array();
    });
    return ret_type(res);
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_m1 = value_of(m1);
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    auto inv_m2 = to_arena(arena_m2.val().array().inverse());
    using val_ret = decltype((inv_m2 * arena_m1.array()).matrix().eval());
    using ret_type = return_var_matrix_t<val_ret, Mat1, Mat2>;
    arena_t<ret_type> res = (inv_m2.array() * arena_m1.array()).matrix();
    reverse_pass_callback([inv_m2, arena_m1, arena_m2, res]() mutable {
      arena_m2.adj().array() -= inv_m2 * res.adj().array() * res.val().array();
    });
    return ret_type(res);
  } else {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<double, Mat2>> arena_m2 = value_of(m2);
    auto inv_m2 = to_arena(arena_m2.array().inverse());
    using val_ret = decltype((inv_m2 * arena_m1.val().array()).matrix().eval());
    using ret_type = return_var_matrix_t<val_ret, Mat1, Mat2>;
    arena_t<ret_type> res = (inv_m2.array() * arena_m1.val().array()).matrix();
    reverse_pass_callback([inv_m2, arena_m1, arena_m2, res]() mutable {
      arena_m1.adj().array() += inv_m2 * res.adj().array();
    });
    return ret_type(res);
  }
}




}  // namespace math
}  // namespace stan
#endif
