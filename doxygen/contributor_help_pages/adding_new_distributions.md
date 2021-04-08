## Adding A New Distribution {#new_distribution}

Before reading this section it's a good idea to at least
skim over the [getting started guide](@ref getting_started). Stan's univariate
distribution functions must work with mixes of scalars and vectors. This
requirement can make the code for the distributions look a bit daunting, but
in the below we'll cover each of the general steps needed to add a new distribution.

1. Get the distribution working in the Stan language.
2. Write out the partial derivatives for each input
3. Writing the function
4. Testing the function

For examples we will use the [loglogistic](https://github.com/stan-dev/math/issues/2320) distribution

### Get the Distribution Working In Stan

In the example for [`loglogistic`] we see the distribution has been written in
the Stan language using custom `lpdf`, `lcdf`, `cdf`, `lccdf` and `rng` functions
in the Stan language.

```stan
real loglogistic_lpdf (real y, real alpha, real beta){
  real lalpha = log(alpha);
  real numerator = log(beta) - lalpha + (beta - 1) * (log(y) - lalpha);
  return numerator - 2 * (log1p( (y/alpha)^beta ));
}

real loglogistic_lcdf (real y, real alpha, real beta) {
  return -log1p((y / alpha) ^-beta);
}

real loglogistic_cdf (real y, real alpha, real beta) {
  return 1/(1 + (y / alpha) ^-beta);
}

real loglogistic_lccdf (real y, real alpha, real beta){
  return -log1p( (y / alpha)^beta);
}

real loglogistic_rng (real alpha, real beta, real lb, real ub) {
  real p_ub = loglogistic_cdf(ub, alpha, beta);
  real p_lb = loglogistic_cdf(lb, alpha, beta);
  real r = uniform_rng(p_lb, p_ub);
  return alpha * (r / (1 - r) )^(1/beta);
}
```

This is nice because now we can test those functions against another implimentation in order
to verify their correctness. Also note in the example above there are custom
distribution functions for each of the `lpdf`, `lcdf`, `cdf`, `lccdf` and `rng`.
This is a necessary condition to accept a new distribution into Stan math and
them in the Stan language makes them easy to test and verify their correctness.
It also makes the next step of writing out the partials a bit easier.

At this point we are going to just focus on the `cdf` function, but the below should apply
to all of the other distribution functions (besides rng which we will handle seperatly)

### Writing out the Partials

For an efficient implimentation of the distribution we want to calculate
each of it's partials with respect to the distributions inputs. This sounds
rough, but for simpler distributions like the above
[matrixcalculus.org](http://www.matrixcalculus.org/) is your friend. There
we can plug in the cdf

$$\begin{aligned}
f = \text{loglogistic_cdf}(y, \alpha, \beta) &= \frac{1}{(1 + (\frac{y}{\alpha}) ^{-\beta})} \cr
\frac{\partial f}{\partial y} &= \frac{(\beta\cdot y^{(-(1+\beta))}\cdot \alpha^{\beta})}{(1+y^{(-\beta)}\cdot \alpha^{\beta})^{2}} \cr
\frac{\partial f}{\partial \alpha} &= -\frac{(\beta\cdot \alpha^{(\beta-1)}\cdot y^{(-\beta)})}{(1+y^{(-\beta)}\cdot \alpha^{\beta})^{2}} \cr
\frac{\partial f}{\partial \beta} &= -\frac{((y^{(-\beta)}\cdot \alpha^{\beta}\cdot \log(\alpha))}{(1+y^{(-\beta)}\cdot \alpha^{\beta})^{2}}-\frac{(\alpha^{\beta}\cdot y^{(-\beta)}\cdot \log(y))}{(1+y^{(-\beta)}\cdot \alpha^{\beta})^{2})}
\end{aligned}$$

Notice there are a good number of common subexpressiosn such as $$y^{-beta} \cdot \alpha^{-beta} $$ and $$ 1+y^{(-\beta)}\cdot \alpha^{\beta})^{2} $$. At this point it's good to sit down and simplify the equations and pull out common subexpressions.

### Writing the function


So now let's add the cdf function in `stan/math/prim/dist/loglogistic_cdf.hpp`. First we'll go over what we have to do before we start doing any math.

#### Working With Scalars and Vectors

```cpp
/** \ingroup prob_dists
 * Docs describing the templates and arguments
 */
template <typename T_y, typename T_alpha, typename T_beta> // (1)
return_type_t<T_y, T_alpha, T_beta> loglogistic_cdf(const T_y& y,
    const T_alpha& alpha,
    const T_beta& beta) {
  // (2) Getting the least upper bound partial
  using T_partials_return = partials_return_t<T_y, T_alpha, T_beta>;
  // (3) Setup types using type traits.
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  using T_beta_ref = ref_type_t<T_beta>;
  // (4) Check consistent sizes
  constexpr const char* function = "loglogistic_cdf";
  check_consistent_sizes(function, "y", y, "Shape", alpha, "Scale", beta);
  if (size_zero(n, N, alpha, beta)) {
    return 1.0;
  }

  // (5) Possibly evaluate Eigen expressions
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  // (6) Check distribution conditions.
  check_nonnegative(function, "y", y_ref);
  check_nonnegative(function, "Shape", alpha_ref);
  check_nonnegative(function, "Scale", beta_ref);
  // (7) logp and setting up operands and partials.
  T_partials_return logp(1.0);
  operands_and_partials<T_y_ref, T_alpha_ref, T_beta_ref> ops_partials(y_ref, alpha_ref, beta_ref);

  // (8) Make scalars and vectors iterable by sequence views.
  scalar_seq_view<T_y_ref> N_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t max_size_seq_view = max_size(y, alpha, beta);
// first half
}
```


Add details of each step here.


#### Assigning the Partials

```cpp
// Second Half
{
// Explicit return for extreme values
// The gradients are technically ill-defined, but treated as zero
for (size_t i = 0; i < stan::math::size(y); i++) {
  if (y_vec.val(i) == NEGATIVE_INFTY) {
    return ops_partials.build(0.0);
  }
}

for (size_t i = 0; i < max_size_seq_view; ++i) {
  auto y_val = y_vec.val(i);
  auto alpha_val = alpha_vec.val(i);
  auto beta_val = beta_vec.val(i);
  if (!is_constant_all<T_y>::value) {
    ops_partials.edge1_.partials_[i] += ;// Partial for y;
  }
  if (!is_constant_all<T_alpha>::value) {
    ops_partials.edge2_.partials_[i] += ;// Partial for alpha;
  }
  if (!is_constant_all<T_beta>::value) {
    ops_partials.edge3_.partials_[i] += ;// Partial for beta;
  }
}

return ops_partials.build(P);
}
```
