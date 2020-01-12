#ifndef STAN_MATH_REV_FUN_DIVIDE_HPP
#define STAN_MATH_REV_FUN_DIVIDE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_var.hpp>
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
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param[in] m specified matrix
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <int R, int C>
inline Eigen::Matrix<var, R, C> divide(const Eigen::Matrix<double, R, C>& m,
                                       const var& c) {
  internal::matrix_scalar_divide_dv_vari<R, C>* baseVari
      = new internal::matrix_scalar_divide_dv_vari<R, C>(m, c);
  Eigen::Matrix<var, R, C> result(m.rows(), m.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, m.rows(), m.cols());
  return result;
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param[in] m specified matrix
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <int R, int C>
inline Eigen::Matrix<var, R, C> divide(const Eigen::Matrix<var, R, C>& m,
                                       const double& c) {
  internal::matrix_scalar_divide_vd_vari<R, C>* baseVari
      = new internal::matrix_scalar_divide_vd_vari<R, C>(m, c);
  Eigen::Matrix<var, R, C> result(m.rows(), m.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, m.rows(), m.cols());
  return result;
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param[in] m specified matrix
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <int R, int C>
inline Eigen::Matrix<var, R, C> divide(const Eigen::Matrix<var, R, C>& m,
                                       const var& c) {
  internal::matrix_scalar_divide_vv_vari<R, C>* baseVari
      = new internal::matrix_scalar_divide_vv_vari<R, C>(m, c);
  Eigen::Matrix<var, R, C> result(m.rows(), m.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, m.rows(), m.cols());
  return result;
}

}  // namespace math
}  // namespace stan
#endif
