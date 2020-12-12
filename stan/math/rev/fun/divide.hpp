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
template <typename Mat, typename = require_eigen_vt<std::is_arithmetic, Mat>>
inline auto divide(const Mat& m, const var& c) {
  auto* baseVari
      = new internal::matrix_scalar_divide_dv_vari<Mat::RowsAtCompileTime,
                                                   Mat::ColsAtCompileTime>(m,
                                                                           c);
  Eigen::Matrix<var, Mat::RowsAtCompileTime, Mat::ColsAtCompileTime> result(
      m.rows(), m.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, m.rows(), m.cols());
  return result;
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat type of the matrix or expression
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, typename = require_eigen_vt<is_var, Mat>>
inline auto divide(const Mat& m, const double& c) {
  auto* baseVari
      = new internal::matrix_scalar_divide_vd_vari<Mat::RowsAtCompileTime,
                                                   Mat::ColsAtCompileTime>(m,
                                                                           c);
  Eigen::Matrix<var, Mat::RowsAtCompileTime, Mat::ColsAtCompileTime> result(
      m.rows(), m.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, m.rows(), m.cols());
  return result;
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat type of the matrix or expression
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, typename = require_eigen_vt<is_var, Mat>,
          typename = void>
inline auto divide(const Mat& m, const var& c) {
  auto* baseVari
      = new internal::matrix_scalar_divide_vv_vari<Mat::RowsAtCompileTime,
                                                   Mat::ColsAtCompileTime>(m,
                                                                           c);
  Eigen::Matrix<var, Mat::RowsAtCompileTime, Mat::ColsAtCompileTime> result(
      m.rows(), m.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, m.rows(), m.cols());
  return result;
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat type of the matrix
 * @tparam Scal type of the scalar
 * @param[in] m input matrix
 * @param[in] c input scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, typename Scal, require_var_matrix_t<Mat>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
inline auto divide(const Mat& m, const Scal& c) {
  double invc = 1.0 / value_of(c);

  plain_type_t<Mat> res = invc * m.val();

  reverse_pass_callback([m, c, res, invc]() mutable {
    m.adj() += invc * res.adj();
    if (!is_constant<Scal>::value)
      forward_as<var>(c).adj()
          -= invc * (res.adj().array() * res.val().array()).sum();
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
