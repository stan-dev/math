# R script to run a stan trest of the interpolate_1d_cubic_spline 
# interpolation functionality when supplied as a local clone 
# of stan math. Probably to not include in the final pull 
# request. Martin Lindén 2026-10-13

require(tidyverse)
library(cmdstanr)
require(posterior)

yfun <- \(x)(sin(2 * pi * x) + sin(4 * pi * x) + 0.7*sin(6 * pi * x))
dydxfun <- \(x)(2*pi*cos(2 * pi * x) + 4*pi*cos(4 * pi * x) + 6*pi*0.7*cos(6 * pi * x))
xy_data <- tibble(x = seq(0, 1, 0.1)) |> mutate(y = yfun(x), dydx=dydxfun(x))

tibble(curve = "true", x = seq(-0.2, 1.2, 0.01)) |>
  mutate( y = yfun(x) ) |> 
  ggplot(aes(x = x, y = y, color = curve)) +
  geom_line() +
  geom_point(data = xy_data |> mutate(curve = "lookup points"))


tibble(curve = "true", x = seq(-0.2, 1.2, 0.01)) |>
  mutate( dydx = dydxfun(x) ) |> 
  ggplot(aes(x = x, y = dydx, color = curve)) +
  geom_line() +
  geom_point(data = xy_data |> mutate(curve = "lookup points"))


mod <- cmdstan_model("hermite_math_test.stan", 
  force_recompile = FALSE,
  user_header = "stan/math/prim/fun/interpolate_1d_cubic_spline.hpp",
  cpp_options = list( include_paths = "math" )
)

# Test quality of interpolation and extrapolation against R's natural cubic spline.
###############################################################################
xy_data <- tibble(x = sort(runif(20,0,1))) |> mutate(y = yfun(x), dydx=dydxfun(x)) # points do not need to be quidistant

x_eval <- seq(min(xy_data$x)+0.05, max(xy_data$x)-0.05, 0.001)

stan_data <- list(
  K = nrow(xy_data),
  xk = xy_data$x,
  yk = xy_data$y,
  dydxk= xy_data$dydx,
  N_eval = length(x_eval),
  x_eval = x_eval
)

fit <- mod$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 0,
  iter_sampling = 1,
  fixed_param = TRUE,
  seed = 123
)

h_func <- splinefunH(stan_data$xk, stan_data$yk, stan_data$dydxk)

test <- tibble(
  x = stan_data$x_eval,
  y.stan = fit |>
    as_draws("y_eval") |>
    as_draws_matrix() |>
    as.numeric(),
    y.r = h_func(stan_data$x_eval)
    #y.r = spline(x = xy_data$x, y = xy_data$y, method = "natural", xout=stan_data$x_eval)$y,
    #y.r = yfun(x_eval)
) |>
  mutate(
    y.diff_r = y.r - y.stan,
    #y.orig = yfun(x_eval)
  )

test |>
  pivot_longer(starts_with("y.")) |>
  mutate(gr = case_when(str_detect(name, "diff") ~ "diff", TRUE ~ "curves")) |>
  ggplot(aes(
    x = x,
    y = value,
    color = name,
    shape = name
  )) +
  facet_wrap(~gr, ncol = 1, scale = "free_y") +
  geom_vline(xintercept = xy_data$x ) +
  geom_line() +
  geom_point(data = xy_data |>
    mutate(gr = "curves", name = "tabulated", value = y))+
    labs(title="comparing Stan and R splines with non-equidistant points")

# support points do not need to be equidistant
xy_data <- tibble(x = sort(runif(20,0,1))) |> mutate(y = yfun(x), dydx=dydxfun(x))

stan_data <- list(
  K = nrow(xy_data),
  xk = xy_data$x,
  yk = xy_data$y,
  dydxk= xy_data$dydx,
  N_eval = length(x_eval),
  x_eval = x_eval
)

fit <- mod$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 0,
  iter_sampling = 1,
  fixed_param = TRUE,
  seed = 123
)
h_func <- splinefunH(stan_data$xk, stan_data$yk, stan_data$dydxk)

test <- tibble(
  x = stan_data$x_eval,
  y.stan = fit |>
    as_draws("y_eval") |>
    as_draws_matrix() |>
    as.numeric(),
    y.r = h_func(x),
) |>
  mutate(
    y.diff_r = y.r - y.stan,
    #y.orig = yfun(x_eval)
  )

test |>
  pivot_longer(starts_with("y.")) |>
  mutate(gr = case_when(str_detect(name, "diff") ~ "diff", TRUE ~ "curves")) |>
  ggplot(aes(
    x = x,
    y = value,
    color = name,
    shape = name
  )) +
  facet_wrap(~gr, ncol = 1, scale = "free_y") +
  geom_vline(xintercept = xy_data$x ) +
  geom_line() +
  geom_point(data = xy_data |>
               mutate(gr = "curves", name = "tabulated", value = y))+
  labs(title="comparing Stan and R splines with non-equidistant points")


# Edge cases #
##############
xy_data <- tibble(x = (runif(20,0,1))) |> mutate(y = yfun(x), dydx=dydxfun(x))
x_eval <- seq(0.05, 0.95, 0.001)

stan_data1 <- list(
  K = nrow(xy_data),
  xk = sample(xy_data$x), # violate strictly increasing order of x
  yk = xy_data$y,
  dydxk = xy_data$y,
  N_eval = length(x_eval),
  x_eval = x_eval
)

fit <- mod$sample(
  data = stan_data1,
  chains = 1,
  iter_warmup = 0,
  iter_sampling = 1,
  fixed_param = TRUE,
  seed = 123
)
