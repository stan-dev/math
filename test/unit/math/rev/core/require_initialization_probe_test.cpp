// Canary for the interaction between Eigen::NumTraits<var> and Eigen 5's
// memset fill fast path. Stock Eigen 5 turns Matrix<var>::setZero() into
// memset(0) because var is trivially copyable, nulling every vari pointer;
// the vendored Eigen carries a Stan-local Fill.h patch (see
// lib/eigen_5.0.1/STAN_CHANGES.md) that disables that path. With
// NumTraits<var>::RequireInitialization = 0 this patch is the only
// protection, so this test must be re-checked whenever Eigen is
// re-vendored: if it fails, the Fill.h patch was lost or the upstream
// gate changed.
#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>

TEST_F(AgradRev, require_initialization_probe) {
  using stan::math::var;
  // preconditions for Eigen's memset path
  static_assert(std::is_trivially_copyable<var>::value,
                "var is trivially copyable");
  constexpr bool memset_path = Eigen::internal::eigen_memset_helper<
      Eigen::Matrix<var, Eigen::Dynamic, 1>>::value;
  ::testing::Test::RecordProperty("memset_path", memset_path);

  Eigen::Matrix<var, Eigen::Dynamic, 1> m(3);
  m.setZero();
  for (int i = 0; i < 3; ++i) {
    EXPECT_NE(m(i).vi_, nullptr) << "setZero() nulled vari at " << i
                                 << " (memset path = " << memset_path << ")";
  }

  Eigen::Matrix<var, Eigen::Dynamic, 1> c(3);
  c.setConstant(var(2.5));
  for (int i = 0; i < 3; ++i) {
    ASSERT_NE(c(i).vi_, nullptr);
    EXPECT_DOUBLE_EQ(c(i).val(), 2.5);
  }
}
