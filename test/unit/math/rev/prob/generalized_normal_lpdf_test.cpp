#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

namespace generalized_normal_lpdf_rev_test {

template <typename T_operand>
double adjoint_or_zero(const T_operand& operand) {
  if constexpr (stan::is_autodiff_v<T_operand>) {
    return operand.adj();
  } else {
    return 0.0;
  }
}

template <size_t autodiff_mask>
void expect_propto_combination(double observation_value, double shape_value) {
  using stan::math::generalized_normal_lpdf;
  using stan::math::var;
  constexpr bool observation_is_autodiff = (autodiff_mask & 1) != 0;
  constexpr bool location_is_autodiff = (autodiff_mask & 2) != 0;
  constexpr bool scale_is_autodiff = (autodiff_mask & 4) != 0;
  constexpr bool shape_is_autodiff = (autodiff_mask & 8) != 0;
  const std::conditional_t<observation_is_autodiff, var, double> observation
      = observation_value;
  const std::conditional_t<location_is_autodiff, var, double> location = 0.4;
  const std::conditional_t<scale_is_autodiff, var, double> scale = 1.7;
  const std::conditional_t<shape_is_autodiff, var, double> shape = shape_value;

  SCOPED_TRACE(autodiff_mask);
  auto propto_log_density
      = generalized_normal_lpdf<true>(observation, location, scale, shape);
  if constexpr (autodiff_mask == 0) {
    EXPECT_DOUBLE_EQ(0.0, propto_log_density);
  } else {
    propto_log_density.grad();
    const std::array<double, 4> propto_gradient{
        adjoint_or_zero(observation), adjoint_or_zero(location),
        adjoint_or_zero(scale), adjoint_or_zero(shape)};
    stan::math::set_zero_all_adjoints();

    auto full_log_density
        = generalized_normal_lpdf(observation, location, scale, shape);
    const double omitted_normalizers
        = stan::math::LOG_TWO + (scale_is_autodiff ? 0.0 : std::log(1.7))
          + (shape_is_autodiff ? 0.0
                               : stan::math::lgamma(1.0 + 1.0 / shape_value));
    EXPECT_NEAR(full_log_density.val() + omitted_normalizers,
                propto_log_density.val(), 1e-12);

    full_log_density.grad();
    EXPECT_NEAR(propto_gradient[0], adjoint_or_zero(observation), 1e-12);
    EXPECT_NEAR(propto_gradient[1], adjoint_or_zero(location), 1e-12);
    EXPECT_NEAR(propto_gradient[2], adjoint_or_zero(scale), 1e-12);
    EXPECT_NEAR(propto_gradient[3], adjoint_or_zero(shape), 1e-12);
  }
  stan::math::recover_memory();
}

template <size_t... autodiff_masks>
void expect_all_propto_combinations(double observation_value,
                                    double shape_value,
                                    std::index_sequence<autodiff_masks...>) {
  (expect_propto_combination<autodiff_masks>(observation_value, shape_value),
   ...);
}

}  // namespace generalized_normal_lpdf_rev_test

TEST_F(AgradRev, generalized_normal_lpdf_propto_combinations) {
  for (const double observation_value : {1.3, 0.4}) {
    for (const double shape_value : {0.7, 2.0, 4.0}) {
      generalized_normal_lpdf_rev_test::expect_all_propto_combinations(
          observation_value, shape_value, std::make_index_sequence<16>{});
    }
  }
}

TEST_F(AgradRev, generalized_normal_lpdf_tiny_location_gradients) {
  for (const double observation_value : {-1e-200, 1e-200}) {
    const stan::math::var observation = observation_value;
    const stan::math::var location = 0.0;
    auto log_density
        = stan::math::generalized_normal_lpdf(observation, location, 1.0, 2.0);
    log_density.grad();

    // The kernel underflows to zero, but its location gradient does not.
    EXPECT_DOUBLE_EQ(-2.0, observation.adj() / observation_value);
    EXPECT_DOUBLE_EQ(2.0, location.adj() / observation_value);
    stan::math::recover_memory();
  }
}

TEST_F(AgradRev, generalized_normal_lpdf_expression_inputs) {
  const Eigen::Vector3d observation_values(-1.0, 0.5, 2.0);
  const Eigen::Vector3d location_values(0.0, 0.2, -0.4);
  const Eigen::Vector3d scale_values(1.0, 2.5, 0.7);
  const Eigen::Vector3d shape_values(3.5, 4.0, 5.5);
  const Eigen::Matrix<stan::math::var, 3, 1> observations
      = observation_values.cast<stan::math::var>();
  const Eigen::Matrix<stan::math::var, 3, 1> locations
      = location_values.cast<stan::math::var>();
  const Eigen::Matrix<stan::math::var, 3, 1> scales
      = scale_values.cast<stan::math::var>();
  const Eigen::Matrix<stan::math::var, 3, 1> shapes
      = shape_values.cast<stan::math::var>();

  auto expression_log_density = stan::math::generalized_normal_lpdf(
      observations.head(2).array() + 0.2, locations.head(2).array() - 0.1,
      scales.head(2).array() + 0.7, shapes.head(2).array() + 0.6);
  expression_log_density.grad();
  const Eigen::Vector3d observation_gradient = observations.adj();
  const Eigen::Vector3d location_gradient = locations.adj();
  const Eigen::Vector3d scale_gradient = scales.adj();
  const Eigen::Vector3d shape_gradient = shapes.adj();
  stan::math::set_zero_all_adjoints();

  const Eigen::Array<stan::math::var, -1, 1> evaluated_observations
      = observations.head(2).array() + 0.2;
  const Eigen::Array<stan::math::var, -1, 1> evaluated_locations
      = locations.head(2).array() - 0.1;
  const Eigen::Array<stan::math::var, -1, 1> evaluated_scales
      = scales.head(2).array() + 0.7;
  const Eigen::Array<stan::math::var, -1, 1> evaluated_shapes
      = shapes.head(2).array() + 0.6;
  auto evaluated_log_density = stan::math::generalized_normal_lpdf(
      evaluated_observations, evaluated_locations, evaluated_scales,
      evaluated_shapes);
  EXPECT_NEAR(expression_log_density.val(), evaluated_log_density.val(), 1e-12);
  evaluated_log_density.grad();
  EXPECT_TRUE(observation_gradient.isApprox(observations.adj(), 1e-12));
  EXPECT_TRUE(location_gradient.isApprox(locations.adj(), 1e-12));
  EXPECT_TRUE(scale_gradient.isApprox(scales.adj(), 1e-12));
  EXPECT_TRUE(shape_gradient.isApprox(shapes.adj(), 1e-12));
}

TEST_F(AgradRev, generalized_normal_lpdf_propto_checks_constants) {
  using stan::math::generalized_normal_lpdf;
  const double infinity = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(generalized_normal_lpdf<true>(nan, 0.0, 1.0, 2.0),
               std::domain_error);
  EXPECT_THROW(generalized_normal_lpdf<true>(1.0, infinity, 1.0, 2.0),
               std::domain_error);
  EXPECT_THROW(generalized_normal_lpdf<true>(1.0, 0.0, 0.0, 2.0),
               std::domain_error);
  EXPECT_THROW(generalized_normal_lpdf<true>(1.0, 0.0, 1.0, infinity),
               std::domain_error);
}
