#ifndef STAN_MATH_REV_MAT_FUN_DIVIDE_HPP
#define STAN_MATH_REV_MAT_FUN_DIVIDE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

namespace internal {
template <int R, int C>
class matrix_scalar_divide_vari_dv : public vari {
 public:
  int rows_;
  int cols_;
  vari* adjCRef_;
  vari** adjResultRef_;
  double invc;

  explicit matrix_scalar_divide_vari_dv(const Eigen::Matrix<double, R, C>& m,
                                        const var& c)
      : vari(0.0),
        rows_(m.rows()),
        cols_(m.cols()),
        adjCRef_(c.vi_),
        adjResultRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())) {
    double c_d = c.val();
    invc = 1.0 / c_d;
    Eigen::Matrix<double, R, C> result = invc * m;
    Eigen::Map<matrix_vi>(adjResultRef_, rows_, cols_)
        = result.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    Eigen::Map<matrix_vi> adjResult(adjResultRef_, rows_, cols_);
    adjCRef_->adj_
        += -1.0 * invc
           * (adjResult.adj().array() * adjResult.val().array()).sum();
  }
};

template <int R, int C>
class matrix_scalar_divide_vari_vd : public vari {
 public:
  int rows_;
  int cols_;
  vari** adjMRef_;
  vari** adjResultRef_;
  double invc;

  explicit matrix_scalar_divide_vari_vd(const Eigen::Matrix<var, R, C>& m,
                                        const double& c)
      : vari(0.0),
        rows_(m.rows()),
        cols_(m.cols()),
        adjMRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())),
        adjResultRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())) {
    Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_) = m.vi();
    invc = 1.0 / c;
    Eigen::Matrix<double, R, C> result = invc * m.val();
    Eigen::Map<matrix_vi>(adjResultRef_, rows_, cols_)
        = result.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    Eigen::Map<matrix_vi> adjM(adjMRef_, rows_, cols_);
    Eigen::Map<matrix_vi> adjResult(adjResultRef_, rows_, cols_);
    adjM.adj() += invc * adjResult.adj();
  }
};

template <int R, int C>
class matrix_scalar_divide_vari_vv : public vari {
 public:
  int rows_;
  int cols_;
  vari** adjMRef_;
  vari* adjC_;
  vari** adjResultRef_;
  double invc;

  explicit matrix_scalar_divide_vari_vv(const Eigen::Matrix<var, R, C>& m,
                                        const var& c)
      : vari(0.0),
        rows_(m.rows()),
        cols_(m.cols()),
        adjMRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())),
        adjC_(c.vi_),
        adjResultRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            m.rows() * m.cols())) {
    Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_) = m.vi();
    double c_d = c.val();
    invc = 1.0 / c_d;
    Eigen::Matrix<double, R, C> result = invc * m.val();
    Eigen::Map<matrix_vi>(adjResultRef_, rows_, cols_)
        = result.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    Eigen::Map<matrix_vi> adjM(adjMRef_, rows_, cols_);
    Eigen::Map<matrix_vi> adjResult(adjResultRef_, rows_, cols_);
    adjC_->adj_ += -1.0 * invc
                   * (adjResult.adj().array() * adjResult.val().array()).sum();
    adjM.adj() += invc * adjResult.adj();
  }
};

}  // namespace internal

/**
 * Return matrix divided by scalar.
 * @tparam R number of rows or Eigen::Dynamic
 * @tparam C number of columns or Eigen::Dynamic
 * @param[in] m specified matrix
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <int R, int C>
inline Eigen::Matrix<var, R, C> divide(const Eigen::Matrix<double, R, C>& m,
                                       const var& c) {
  internal::matrix_scalar_divide_vari_dv<R, C>* baseVari
      = new internal::matrix_scalar_divide_vari_dv<R, C>(m, c);
  Eigen::Matrix<var, R, C> result(m.rows(), m.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, m.rows(), m.cols());
  return result;
}

/**
 * Return matrix divided by scalar.
 * @tparam R number of rows or Eigen::Dynamic
 * @tparam C number of columns or Eigen::Dynamic
 * @param[in] m specified matrix
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <int R, int C>
inline Eigen::Matrix<var, R, C> divide(const Eigen::Matrix<var, R, C>& m,
                                       const double& c) {
  internal::matrix_scalar_divide_vari_vd<R, C>* baseVari
      = new internal::matrix_scalar_divide_vari_vd<R, C>(m, c);
  Eigen::Matrix<var, R, C> result(m.rows(), m.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, m.rows(), m.cols());
  return result;
}

/**
 * Return matrix divided by scalar.
 * @tparam R number of rows or Eigen::Dynamic
 * @tparam C number of columns or Eigen::Dynamic
 * @param[in] m specified matrix
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <int R, int C>
inline Eigen::Matrix<var, R, C> divide(const Eigen::Matrix<var, R, C>& m,
                                       const var& c) {
  internal::matrix_scalar_divide_vari_vv<R, C>* baseVari
      = new internal::matrix_scalar_divide_vari_vv<R, C>(m, c);
  Eigen::Matrix<var, R, C> result(m.rows(), m.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, m.rows(), m.cols());
  return result;
}

}  // namespace math
}  // namespace stan
#endif
