#ifndef STAN_MATH_EIGEN_PLUGINS_H
#define STAN_MATH_EIGEN_PLUGINS_H

#include "plugins/typedefs.h"
#include "plugins/adj_view.h"
#include "plugins/val_view.h"
#include "plugins/d_view.h"
#include "plugins/vi_view.h"

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const Scalar& coeffRef(Index row, Index col) const {
    eigen_internal_assert(row >= 0 && row < rows()
                          && col >= 0 && col < cols());
    return internal::evaluator<Derived>(derived()).coeffRef(row, col);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
const Scalar& coeffRef(Index index) const {
    eigen_internal_assert(index >= 0 && index < size());
    return internal::evaluator<Derived>(derived()).coeffRef(index);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
Scalar& coeffRef(Index row, Index col) {
    eigen_internal_assert(row >= 0 && row < rows()
                          && col >= 0 && col < cols());
    return internal::evaluator<Derived>(derived()).coeffRef(row, col);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
Scalar& coeffRef(Index index) {
    eigen_internal_assert(index >= 0 && index < size());
    return internal::evaluator<Derived>(derived()).coeffRef(index);
}

#define EIGEN_STAN_DENSEBASE_PLUGIN
#endif
