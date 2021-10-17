
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

#define EIGEN_STAN_MATRIXBASE_PLUGIN
#define EIGEN_STAN_ARRAYBASE_PLUGIN
