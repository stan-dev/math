#ifndef TEST_UNIT_MATH_PRIM_FUN_TERNARY_SCALAR_TESTER_HPP
#define TEST_UNIT_MATH_PRIM_FUN_TERNARY_SCALAR_TESTER_HPP
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

namespace stan {
namespace test {

/**
 * Implementation function which checks that the ternary vectorisation
 * framework returns the same value as the function with scalar inputs,
 * for all valid combinations of scalar/vector/nested vector.
 *
 * Specialization for use with matrix inputs
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first matrix.
 * @tparam T2 Type of second matrix.
 * @tparam T3 Type of third matrix.
 * @param x First matrix input to which operation is applied.
 * @param y Second matrix input to which operation is applied.
 * @param z Third matrix input to which operation is applied.
 * @param f functor to apply to inputs.
 */
template <typename F, typename T1, typename T2, typename T3,
          require_all_not_vector_t<T1, T2, T3>* = nullptr>
void ternary_scalar_tester_impl(const F& f, const T1& x, const T2& y,
                                const T3& z) {
  auto mat_mat_mat = math::eval(f(x, y, z));
  auto mat_mat_scal = math::eval(f(x, y, z(0)));
  auto mat_scal_scal = math::eval(f(x, y(0), z(0)));
  auto mat_scal_mat = math::eval(f(x, y(0), z));
  auto scal_scal_mat = math::eval(f(x(0), y(0), z));
  auto scal_mat_scal = math::eval(f(x(0), y, z(0)));
  auto scal_mat_mat = math::eval(f(x(0), y, z));

  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(f(x(i), y(i), z(i)), mat_mat_mat(i));
    EXPECT_FLOAT_EQ(f(x(i), y(i), z(0)), mat_mat_scal(i));
    EXPECT_FLOAT_EQ(f(x(i), y(0), z(0)), mat_scal_scal(i));
    EXPECT_FLOAT_EQ(f(x(i), y(0), z(i)), mat_scal_mat(i));
    EXPECT_FLOAT_EQ(f(x(0), y(0), z(i)), scal_scal_mat(i));
    EXPECT_FLOAT_EQ(f(x(0), y(i), z(0)), scal_mat_scal(i));
    EXPECT_FLOAT_EQ(f(x(0), y(i), z(i)), scal_mat_mat(i));
  }
  plain_type_t<T1> x_zero;
  plain_type_t<T2> y_zero;
  plain_type_t<T3> z_zero;
  EXPECT_THROW(f(x_zero, y, z), std::invalid_argument);
  EXPECT_THROW(f(x, y_zero, z), std::invalid_argument);
  EXPECT_THROW(f(x, y, z_zero), std::invalid_argument);

  std::vector<T1> nest_x{x, x, x};
  std::vector<T2> nest_y{y, y, y};
  std::vector<T3> nest_z{z, z, z};
  auto nestmat_nestmat_nestmat = f(nest_x, nest_y, nest_z);
  auto nestmat_nestmat_scal = f(nest_x, nest_y, z(0));
  auto nestmat_scal_scal = f(nest_x, y(0), z(0));
  auto nestmat_scal_nestmat = f(nest_x, y(0), nest_z);
  auto scal_scal_nestmat = f(x(0), y(0), nest_z);
  auto scal_nestmat_scal = f(x(0), nest_y, z(0));
  auto scal_nestmat_nestmat = f(x(0), nest_y, nest_z);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < x.size(); ++j) {
      EXPECT_FLOAT_EQ(f(nest_x[i](j), nest_y[i](j), nest_z[i](j)),
                      nestmat_nestmat_nestmat[i](j));
      EXPECT_FLOAT_EQ(f(nest_x[i](j), nest_y[i](j), nest_z[i](0)),
                      nestmat_nestmat_scal[i](j));
      EXPECT_FLOAT_EQ(f(nest_x[i](j), nest_y[i](0), nest_z[i](0)),
                      nestmat_scal_scal[i](j));
      EXPECT_FLOAT_EQ(f(nest_x[i](j), nest_y[i](0), nest_z[i](j)),
                      nestmat_scal_nestmat[i](j));
      EXPECT_FLOAT_EQ(f(nest_x[i](0), nest_y[i](0), nest_z[i](j)),
                      scal_scal_nestmat[i](j));
      EXPECT_FLOAT_EQ(f(nest_x[i](0), nest_y[i](j), nest_z[i](0)),
                      scal_nestmat_scal[i](j));
      EXPECT_FLOAT_EQ(f(nest_x[i](0), nest_y[i](j), nest_z[i](j)),
                      scal_nestmat_nestmat[i](j));
    }
  }
  std::vector<T1> nest_x_small{x, x};
  std::vector<T2> nest_y_small{y, y};
  std::vector<T3> nest_z_small{z, z};
  EXPECT_THROW(f(nest_x_small, nest_y, nest_z), std::invalid_argument);
  EXPECT_THROW(f(nest_x, nest_y_small, nest_z), std::invalid_argument);
  EXPECT_THROW(f(nest_x, nest_y, nest_z_small), std::invalid_argument);

  std::vector<std::vector<T1>> nest_nest_x{nest_x, nest_x, nest_x};
  std::vector<std::vector<T2>> nest_nest_y{nest_y, nest_y, nest_y};
  std::vector<std::vector<T3>> nest_nest_z{nest_z, nest_z, nest_z};
  auto nestnestmat_nestnestmat_nestnestmat
      = f(nest_nest_x, nest_nest_y, nest_nest_z);
  auto nestnestmat_nestnestmat_scal = f(nest_nest_x, nest_nest_y, z(0));
  auto nestnestmat_scal_scal = f(nest_nest_x, y(0), z(0));
  auto nestnestmat_scal_nestnestmat = f(nest_nest_x, y(0), nest_nest_z);
  auto scal_scal_nestnestmat = f(x(0), y(0), nest_nest_z);
  auto scal_nestnestmat_scal = f(x(0), nest_nest_y, z(0));
  auto scal_nestnestmat_nestnestmat = f(x(0), nest_nest_y, nest_nest_z);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < x.size(); ++k) {
        EXPECT_FLOAT_EQ(
            f(nest_nest_x[i][j](k), nest_nest_y[i][j](k), nest_nest_z[i][j](k)),
            nestnestmat_nestnestmat_nestnestmat[i][j](k));
        EXPECT_FLOAT_EQ(f(nest_nest_x[i][j](k), nest_nest_y[i][j](k), z(0)),
                        nestnestmat_nestnestmat_scal[i][j](k));
        EXPECT_FLOAT_EQ(f(nest_nest_x[i][j](k), y(0), z(0)),
                        nestnestmat_scal_scal[i][j](k));
        EXPECT_FLOAT_EQ(f(nest_nest_x[i][j](k), y(0), nest_nest_z[i][j](k)),
                        nestnestmat_scal_nestnestmat[i][j](k));
        EXPECT_FLOAT_EQ(f(x(0), y(0), nest_nest_z[i][j](k)),
                        scal_scal_nestnestmat[i][j](k));
        EXPECT_FLOAT_EQ(f(x(0), nest_nest_y[i][j](k), z(0)),
                        scal_nestnestmat_scal[i][j](k));
        EXPECT_FLOAT_EQ(f(x(0), nest_nest_y[i][j](k), nest_nest_z[i][j](k)),
                        scal_nestnestmat_nestnestmat[i][j](k));
      }
    }
  }
  std::vector<std::vector<T1>> nest_nest_x_small{nest_x, nest_x};
  std::vector<std::vector<T2>> nest_nest_y_small{nest_y, nest_y};
  std::vector<std::vector<T3>> nest_nest_z_small{nest_z, nest_z};
  EXPECT_THROW(f(nest_nest_x_small, nest_nest_y, nest_nest_z),
               std::invalid_argument);
  EXPECT_THROW(f(nest_nest_x, nest_nest_y_small, nest_nest_z),
               std::invalid_argument);
  EXPECT_THROW(f(nest_nest_x, nest_nest_y, nest_nest_z_small),
               std::invalid_argument);
}

/**
 * Implementation function which checks that the ternary vectorisation
 * framework returns the same value as the function with scalar inputs,
 * for all valid combinations of scalar/vector/nested vector.
 *
 * This is a specialisation for vector inputs.
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first vector.
 * @tparam T2 Type of second vector.
 * @tparam T3 Type of third vector.
 * @param x First vector input to which operation is applied.
 * @param y Second vector input to which operation is applied.
 * @param z Third vector input to which operation is applied.
 * @param f functor to apply to inputs.
 */
template <typename F, typename T1, typename T2, typename T3,
          require_all_vector_t<T1, T2, T3>* = nullptr>
void ternary_scalar_tester_impl(const F& f, const T1& x, const T2& y,
                                const T3& z) {
  auto mat_mat_mat = math::eval(f(x, y, z));
  auto mat_mat_scal = math::eval(f(x, y, z[0]));
  auto mat_scal_scal = math::eval(f(x, y[0], z[0]));
  auto mat_scal_mat = math::eval(f(x, y[0], z));
  auto scal_scal_mat = math::eval(f(x[0], y[0], z));
  auto scal_mat_scal = math::eval(f(x[0], y, z[0]));
  auto scal_mat_mat = math::eval(f(x[0], y, z));

  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(f(x[i], y[i], z[i]), mat_mat_mat[i]);
    EXPECT_FLOAT_EQ(f(x[i], y[i], z[0]), mat_mat_scal[i]);
    EXPECT_FLOAT_EQ(f(x[i], y[0], z[0]), mat_scal_scal[i]);
    EXPECT_FLOAT_EQ(f(x[i], y[0], z[i]), mat_scal_mat[i]);
    EXPECT_FLOAT_EQ(f(x[0], y[0], z[i]), scal_scal_mat[i]);
    EXPECT_FLOAT_EQ(f(x[0], y[i], z[0]), scal_mat_scal[i]);
    EXPECT_FLOAT_EQ(f(x[0], y[i], z[i]), scal_mat_mat[i]);
  }
  plain_type_t<T1> x_zero;
  plain_type_t<T2> y_zero;
  plain_type_t<T3> z_zero;
  EXPECT_THROW(f(x_zero, y, z), std::invalid_argument);
  EXPECT_THROW(f(x, y_zero, z), std::invalid_argument);
  EXPECT_THROW(f(x, y, z_zero), std::invalid_argument);

  std::vector<T1> nest_x{x, x, x};
  std::vector<T2> nest_y{y, y, y};
  std::vector<T3> nest_z{z, z, z};
  auto nestmat_nestmat_nestmat = f(nest_x, nest_y, nest_z);
  auto nestmat_nestmat_scal = f(nest_x, nest_y, z[0]);
  auto nestmat_scal_scal = f(nest_x, y[0], z[0]);
  auto nestmat_scal_nestmat = f(nest_x, y[0], nest_z);
  auto scal_scal_nestmat = f(x[0], y[0], nest_z);
  auto scal_nestmat_scal = f(x[0], nest_y, z[0]);
  auto scal_nestmat_nestmat = f(x[0], nest_y, nest_z);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < x.size(); ++j) {
      EXPECT_FLOAT_EQ(f(nest_x[i][j], nest_y[i][j], nest_z[i][j]),
                      nestmat_nestmat_nestmat[i][j]);
      EXPECT_FLOAT_EQ(f(nest_x[i][j], nest_y[i][j], z[0]),
                      nestmat_nestmat_scal[i][j]);
      EXPECT_FLOAT_EQ(f(nest_x[i][j], y[0], z[0]), nestmat_scal_scal[i][j]);
      EXPECT_FLOAT_EQ(f(nest_x[i][j], y[0], nest_z[i][j]),
                      nestmat_scal_nestmat[i][j]);
      EXPECT_FLOAT_EQ(f(x[0], y[0], nest_z[i][j]), scal_scal_nestmat[i][j]);
      EXPECT_FLOAT_EQ(f(x[0], nest_y[i][j], z[0]), scal_nestmat_scal[i][j]);
      EXPECT_FLOAT_EQ(f(x[0], nest_y[i][j], nest_z[i][j]),
                      scal_nestmat_nestmat[i][j]);
    }
  }
  std::vector<T1> nest_x_small{x, x};
  std::vector<T2> nest_y_small{y, y};
  std::vector<T3> nest_z_small{z, z};
  EXPECT_THROW(f(nest_x_small, nest_y, nest_z), std::invalid_argument);
  EXPECT_THROW(f(nest_x, nest_y_small, nest_z), std::invalid_argument);
  EXPECT_THROW(f(nest_x, nest_y, nest_z_small), std::invalid_argument);

  std::vector<std::vector<T1>> nest_nest_x{nest_x, nest_x, nest_x};
  std::vector<std::vector<T2>> nest_nest_y{nest_y, nest_y, nest_y};
  std::vector<std::vector<T3>> nest_nest_z{nest_z, nest_z, nest_z};
  auto nestnestmat_nestnestmat_nestnestmat
      = f(nest_nest_x, nest_nest_y, nest_nest_z);
  auto nestnestmat_nestnestmat_scal = f(nest_nest_x, nest_nest_y, z[0]);
  auto nestnestmat_scal_scal = f(nest_nest_x, y[0], z[0]);
  auto nestnestmat_scal_nestnestmat = f(nest_nest_x, y[0], nest_nest_z);
  auto scal_scal_nestnestmat = f(x[0], y[0], nest_nest_z);
  auto scal_nestnestmat_scal = f(x[0], nest_nest_y, z[0]);
  auto scal_nestnestmat_nestnestmat = f(x[0], nest_nest_y, nest_nest_z);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < x.size(); ++k) {
        EXPECT_FLOAT_EQ(
            f(nest_nest_x[i][j][k], nest_nest_y[i][j][k], nest_nest_z[i][j][k]),
            nestnestmat_nestnestmat_nestnestmat[i][j][k]);
        EXPECT_FLOAT_EQ(f(nest_nest_x[i][j][k], nest_nest_y[i][j][k], z[0]),
                        nestnestmat_nestnestmat_scal[i][j][k]);
        EXPECT_FLOAT_EQ(f(nest_nest_x[i][j][k], y[0], z[0]),
                        nestnestmat_scal_scal[i][j][k]);
        EXPECT_FLOAT_EQ(f(nest_nest_x[i][j][k], y[0], nest_nest_z[i][j][k]),
                        nestnestmat_scal_nestnestmat[i][j][k]);
        EXPECT_FLOAT_EQ(f(x[0], y[0], nest_nest_z[i][j][k]),
                        scal_scal_nestnestmat[i][j][k]);
        EXPECT_FLOAT_EQ(f(x[0], nest_nest_y[i][j][k], z[0]),
                        scal_nestnestmat_scal[i][j][k]);
        EXPECT_FLOAT_EQ(f(x[0], nest_nest_y[i][j][k], nest_nest_z[i][j][k]),
                        scal_nestnestmat_nestnestmat[i][j][k]);
      }
    }
  }
  std::vector<std::vector<T1>> nest_nest_x_small{nest_x, nest_x};
  std::vector<std::vector<T2>> nest_nest_y_small{nest_y, nest_y};
  std::vector<std::vector<T3>> nest_nest_z_small{nest_z, nest_z};
  EXPECT_THROW(f(nest_nest_x_small, nest_nest_y, nest_nest_z),
               std::invalid_argument);
  EXPECT_THROW(f(nest_nest_x, nest_nest_y_small, nest_nest_z),
               std::invalid_argument);
  EXPECT_THROW(f(nest_nest_x, nest_nest_y, nest_nest_z_small),
               std::invalid_argument);
}

/**
 * Testing framework for checking that the vectorisation of ternary
 * functions returns the same results as the ternary function with
 * scalar inputs. This framework takes three Eigen column vectors of
 * inputs which are tested, and then transformed to Eigen row vectors,
 * Eigen matrices, Eigen arrays, Eigen expressions, and std::vectors to
 * also be tested.
 *
 * @tparam F Type of functor to apply.
 * @tparam T1 Type of first Eigen column-vector.
 * @tparam T2 Type of second Eigen column-vector.
 * @tparam T3 Type of third Eigen column-vector.
 * @param x First Eigen column-vector input to which operation is applied.
 * @param y Second Eigen column-vector input to which operation is applied.
 * @param z Third Eigen column-vector input to which operation is applied.
 * @param f functor to apply to inputs.
 */
template <typename F, typename T1, typename T2, typename T3,
          require_all_eigen_col_vector_t<T1, T2, T3>* = nullptr>
void ternary_scalar_tester(const F& f, const T1& x, const T2& y, const T3& z) {
  ternary_scalar_tester_impl(f, x, y, z);
  ternary_scalar_tester_impl(f, x.transpose().eval(), y.transpose().eval(),
                             z.transpose().eval());
  ternary_scalar_tester_impl(f, x.replicate(1, x.size()).eval(),
                             y.replicate(1, y.size()).eval(),
                             z.replicate(1, z.size()).eval());
  ternary_scalar_tester_impl(f, x.replicate(1, x.size()).array().eval(),
                             y.replicate(1, y.size()).array().eval(),
                             z.replicate(1, z.size()).array().eval());
  ternary_scalar_tester_impl(f, x.transpose(), y.transpose(), z.transpose());
  ternary_scalar_tester_impl(
      f, std::vector<typename T1::Scalar>(x.data(), x.data() + x.size()),
      std::vector<typename T2::Scalar>(y.data(), y.data() + y.size()),
      std::vector<typename T3::Scalar>(z.data(), z.data() + z.size()));
}

}  // namespace test
}  // namespace stan
#endif
