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


### Get the Distribution Working In Stan

We will use the normal distribution as an example and adding the lpdf
function, though note for acceptance into Stan math a function must have
its respective `lpdf`, `lcdf`, `cdf`, `lccdf` and `rng` implemented.
Though we will only be doing the `lpdf` in the below all of the notes here will apply
to the other functions.

So for the normal distribution probability density function

\f[
\text{Normal}(y|\mu,\sigma)=\frac{1}{\sigma \sqrt{2\pi}} e^{-\frac{1}{2}\left(\frac{y-\mu}{\sigma}\right)^2}
\f]

to get the log probability density function we log the above to get

\f[
\ln{\left(\text{Normal}(y|\mu,\sigma)\right)}=-\frac{1}{2} \left(\frac{y-\mu}{\sigma}\right)^2 - \ln{\left(\sigma\right)} - \frac{1}{2}\ln{\left(2\pi\right)}
\f]

Now we can directly plug this into Stan as a custom lpdf function

```stan
real new_normal_lpdf (real y, real mu, real sigma){
  return -0.5 * pow((y - mu) / sigma, 2) - log(sigma) - 0.5 * log(2*pi());
}
```

This is nice because now we can test this function against another implementations
to verify its correctness. At this point it is a good idea to post something on
[discourse](https://discourse.mc-stan.org/) or file an
[issue](https://github.com/stan-dev/math/issues) to let folks know you would like
to add this distribution.


### Writing out the Partials

For an efficient implimentation of the distribution we want to calculate
each of its partials with respect to the distributions inputs. This is easy for normal but can be rough for other distributions.
In that case then [matrixcalculus.org](http://www.matrixcalculus.org/) or [wolframalpha](https://www.wolframalpha.com/) is your friend.
There we can plug in the lpdf and get back each of the partials.

Note that for univariate distributions wolfram handles things a bit simpler,
though for multivariate distributions you'll have to use matrixcalculus. Though
for both you'll have a much nicer time if you plug in the log'd version of the
function. One other nice thing about the matrixcalculus site is that it can generate latex
which is nice for documentation.


\f{aligned}{
f = \text{ln}\left(\text{Normal}(y|\mu,\sigma)\right) &= -\frac{1}{2} \left(\frac{y-\mu}{\sigma}\right)^2 - \ln\left(\sigma\right) - \frac{1}{2}\ln{\left(2\pi\right)} \cr
\frac{\partial f}{\partial y} &= -\frac{y-\mu}{\sigma^{2}} \cr
\frac{\partial f}{\partial \mu} &= \frac{y-\mu}{\sigma^{2}} \cr
\frac{\partial f}{\partial \sigma} &= -\frac{1}{\sigma} + \frac{(y-\mu)^{2}}{\sigma^{3}}
\f}

It's a little early, but once we get the `lpdf` function working with the above we will want to get out a pen and paper to simplify and find common subexpressions we only need to calculate once.
For instance in the normal we can compute `y - mu` and `1/sigma`

\f{aligned}{
f(y|\mu,\sigma) = \text{ln}\left(\text{Normal}(y|\mu,\sigma)\right) &= -\frac{1}{2} \left(\frac{y-\mu}{\sigma}\right)^2 - \ln{\left(\sigma\right)} - \frac{1}{2}\ln{\left(2\pi\right)} \cr
\frac{\partial f}{\partial y} &= -t_3 \cr
\frac{\partial f}{\partial \mu} &= t_3 \cr
\frac{\partial f}{\partial \sigma} &= \frac{t_{2}^2}{t_1} \cdot t_0 - t_0 \cr
\text{Where} \cr
t_0 &= \frac{1}{\sigma} \cr
t_1 &= t_{0}^2 \cr
t_2 &= y - \mu \cr
t_3 &= \frac{t_2}{t_1}
\f}


### Writing the function


So now let's add the lpdf function in `stan/math/prim/prob/new_normal_lpdf.hpp`.
First we'll go over what we have to do before we start doing any math. We'll be
breaking down Stan's current `normal_lpdf` function which you can find [here](https://github.com/stan-dev/math/blob/develop/stan/math/prim/prob/normal_lpdf.hpp).


#### Distribution Signature

```cpp
/** \ingroup prob_dists
 * Docs describing the templates and arguments
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lpdf(T_y&& y, T_loc&& mu,
                                                      T_scale&& sigma) {}
```

Each of the input arguments represent the inputs to the `Normal` function we wrote out above.
The arguments are taken as forwarding references (`T_y&&`) so that temporaries
can be moved into the function and passed along with `std::forward` instead of copied.
The `require_*` template parameter keeps this overload from being chosen for
OpenCL kernel generator expressions, which have their own implementation.
The template parameters for univariate distributions are very general and they
must work for all of Stan's scalar types
`double`, `var`, and `fvar<T>` while accepting mixtures of scalars and vectors.
The return of the function is the joint log probability accumulated over all
of the inputs which is a scalar of the least upper bound of all the
parameters scalar types. That's a lot of big words, but in essence means that
if one of the inputs is a `var` and the others are double the return type needs to be
a `var`. If the input signature contained `fvar<var>`, `var`, `double` then the
return type would be `fvar<var>`. See the [Common pitfalls](@ref common_pitfalls) for an
explanation of `return_type_t`.

Notice the `bool propto` template parameter, this is used by the function to
decide whether or not the function needs to propagate constants to the joint
log probability we'll be calculating.

#### Preparing the Parameters

At the start of the function we need to take each argument, deduce whether
it is an unevaluated Eigen expression, and extract the values from them and then
convert them into `Eigen::Array` types. `ref_type_t` is the return type of `to_ref()`
which is explained in the [getting started guide](@ref getting_started).
`ref_type_if_not_constant_t<T>` will evaluate an Eigen expression only if `T`
is an Eigen expression _and_ it holds autodiff types (it is not constant), since
only then will the argument be read more than once. We also declare the name of
the function, which is used in the error messages of the `check_*` functions.

```cpp
// Making aliases for partials and unevaluated expressions
using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
using T_y_ref = ref_type_if_not_constant_t<T_y>;
using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
static constexpr const char* function = "normal_lpdf";
// Evaluating unevaluated eigen expressions
T_y_ref y_ref = std::forward<T_y>(y);
T_mu_ref mu_ref = std::forward<T_loc>(mu);
T_sigma_ref sigma_ref = std::forward<T_scale>(sigma);
// Extracting values from arguments
decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
```


#### Checking Correctness of the Inputs and Early Return

Then we need to check that all the vector inputs sizes match, and then check
that each of the inputs satisfies the conditions of the distribution. For the
normal distribution we need to check that `y` does not contain nan values,
`mu` is finite, and `sigma` is positive. In `normal_lpdf` the
`check_consistent_sizes()` call comes before the arguments are forwarded into
`y_ref`, `mu_ref` and `sigma_ref`, since it only needs their sizes.

```cpp
check_consistent_sizes(function, "Random variable", y, "Location parameter",
                       mu, "Scale parameter", sigma);
check_not_nan(function, "Random variable", y_val);
check_finite(function, "Location parameter", mu_val);
check_positive(function, "Scale parameter", sigma_val);
if (size_zero(y_ref, mu_ref, sigma_ref)) {
  return 0.0;
}
if constexpr (!include_summand<propto, T_y, T_loc, T_scale>::value) {
  return 0.0;
}

```

The `if` statements here are checking if

1. Any of the inputs are length zero
2. Either the function drops constants `propto=true` and all of the inputs are constant (aka if they are all of type `double`).

If either of the two conditions are met then there's no need to calculate the
rest of the lpdf function and we can return back zero.

#### Actually Doing The Math

Woof! That was a good bit of stuff just to get to the math, but here we are!
Our goal is to calculate the partial adjoints using our stuff above, but we only
want to bother ourselves to calculate adjoints of parameters
which are not constant (`double`). There's some more technical bits to building
the log joint probability, but those are all hidden away in the
`partials_propagator` class so we won't cover those here. For now you can take
the evaluated inputs and pass them to `make_partials_propagator()`

```cpp
  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);
```

This sets up each of the input operand's partials so that we only store and calculate
the ones we need.

#### The partials_propagator API

`make_partials_propagator()` lives in `stan/math/prim/functor/partials_propagator.hpp`,
with the reverse and forward mode specializations in
`stan/math/rev/functor/partials_propagator.hpp` and
`stan/math/fwd/functor/partials_propagator.hpp`. The object it returns has one
edge per operand, in the order they were passed. The edges are accessed with
free functions whose index is 0-based, so the first operand is edge 0.

- `partials<I>(ops_partials)` returns the partials of the `I`th operand. For a
  scalar operand this is a scalar; for a container it is a container of the
  same shape.
- `partials_vec<I>(ops_partials)` returns the partials of the `I`th operand
  when that operand is a `std::vector` of vectors or matrices, as in the
  multivariate distributions (see `dirichlet_lpdf`). Index it with `[i]` to get
  the partials of the `i`th element.
- `edge<I>(ops_partials)` returns the whole edge for the rare cases where you
  need more than its partials.
- `ops_partials.build(logp)` returns the final value with the partials
  attached, as a `double`, `var` or `fvar` depending on the operand types.

For operands which are constant (`double`) the edges are empty and writes to
them compile away. Older code uses `operands_and_partials` with `edge1_`,
`edge2_`, ... members. That class is legacy, so use `partials_propagator` in
new code.


There's two ways of doing the math, one using a simple loop and another utilizing
Eigen expressions. Below I will cover both starting with the loop version.

#### Doing The Math With A Loop

The loop version requires one other piece of overhead to make sure that vectors
and scalars can both be iterated over in the loop.

```cpp
  // Make scalars and vectors iterable
  scalar_seq_view<decltype(y_val)> y_vec(y_val);
  scalar_seq_view<decltype(mu_val)> mu_vec(mu_val);
  scalar_seq_view<decltype(sigma_val)> sigma_vec(sigma_val);
```

For vectors, `scalar_seq_view` simply holds a reference to the vector it's passed
and calling `scalar_seq_view`'s method `.val(i)` will return element i in the vector
after calling `value_of()` on the element. The actual element can be accessed
with `operator[]`. For scalars, `scalar_seq_view`'s `.val(i)` and `operator[]`
will just return the scalar no matter what index is passed.

But with that now we can get the maximum size of the input arguments and run a
loop calculating the partials for each input argument's values.

```cpp
size_t N = max_size(y_ref, mu_ref, sigma_ref);
// Stores the accumulated value from the lpdf from operands
T_partials_return logp(0.0);
constexpr double NEGATIVE_HALF = -0.5;
// Include constant if user asked for them.
if constexpr (include_summand<propto>::value) {
  logp += NEG_LOG_SQRT_TWO_PI * N;
}
for (size_t n = 0; n < N; n++) {
  // Do the intermediate calculations from above
  const T_partials_return y_dbl = y_vec.val(n);
  const T_partials_return mu_dbl = mu_vec.val(n);
  const T_partials_return inv_sigma = 1.0 / sigma_vec.val(n);
  const T_partials_return log_sigma = log(sigma_vec.val(n));

  const T_partials_return y_minus_mu_over_sigma
      = (y_dbl - mu_dbl) * inv_sigma;
  const T_partials_return y_minus_mu_over_sigma_squared
      = y_minus_mu_over_sigma * y_minus_mu_over_sigma;

  // Include constants if user asked for them.
  if constexpr (include_summand<propto, T_scale>::value) {
    logp -= log_sigma;
  }
  logp += NEGATIVE_HALF * y_minus_mu_over_sigma_squared;
  // Add partial calculations to each edge
  T_partials_return scaled_diff = inv_sigma * y_minus_mu_over_sigma;
  if constexpr (is_autodiff_v<T_y>) {
    partials<0>(ops_partials)[n] -= scaled_diff;
  }
  if constexpr (is_autodiff_v<T_loc>) {
    partials<1>(ops_partials)[n] += scaled_diff;
  }
  if constexpr (is_autodiff_v<T_scale>) {
    partials<2>(ops_partials)[n]
        += -inv_sigma + inv_sigma * y_minus_mu_over_sigma_squared;
  }
}
// return `logp` and handle rules for propagating partials for each autodiff type.
return ops_partials.build(logp);
```

The `logp` is used to accumulate the log probability density function's value,
where `propto` is used to decide whether or not that value should have constants
added or dropped.

The odd bits here are mostly the `if constexpr`s that include
`include_summand<propto>` and `is_autodiff_v<T_loc>`. We want to
only compute the partials and accumulate the constants if those values are not
constant (`double`). Since the conditionals are type traits whose values are known
at compile time, `if constexpr` discards the false branches before they are
compiled, so we don't pay for any of them when the inputs are constant.

We collect the partials for each of our inputs with `partials<I>(ops_partials)`,
where `I` is the 0-based position of the operand passed to
`make_partials_propagator()`. The first argument is `partials<0>`, the
second `partials<1>` and so on. One important question to ask here is, what if the edge is a
scalar? It seems odd that we are able to call `[n]` on the partials when the operand can be
either a vector or scalar. Under the hood, `partials_propagator` wraps the partials for `Scalar` types
in what's called a [`broadcast_array`](https://github.com/stan-dev/math/blob/develop/stan/math/prim/functor/broadcast_array.hpp)
which has an overloaded `operator[]` for scalars such that it just simply returns back the partials scalar.
Similarly, `broadcast_array` has an overloaded `operator=` which when assigning a vector to the partial the overloaded `operator=`
will sum the vector before assigning it to the partial.

```cpp
if constexpr (is_autodiff_v<T_loc>) {
  // pretend the partials are a scalar and scaled_diff is a vector
  partials<1>(ops_partials) = scaled_diff;
}
```

Finally once the loop is finished we call `ops_partials.build()` passing it
the joint log probability value. For reverse mode this will place a callback
on the callback stack that takes the partials of each edge and accumulates
them into the operands' adjoints.

The for loop version is nice and simple, but there's a few things for performance
that we can do better. For instance, in the for loop version we are constantly
reading and writing to memory from a bunch of different places. We can fix that by
rewriting the above to use multiple loops, but unless we have separate loops
that turn on and off for when combinations of partials need to be calculated
then we lose places where we can share calculations between partials.

For a more efficient version we can do the math for the partials with no loop
and _possibly_ sharing computation.


#### Doing The Math With Eigen

The below code replaces the loop above with Eigen. It uses most of the same tricks
we've used previously.

```cpp
  // Only evaluate inv_sigma here if it's going to be used more than once
  const auto& inv_sigma
      = to_ref_if<is_any_autodiff_v<T_y, T_scale, T_loc>>(inv(sigma_val));
  const auto& y_scaled = to_ref((y_val - mu_val) * inv_sigma);
  // Only evaluate y_scaled_sq here if T_scale is not constant
  const auto& y_scaled_sq
      = to_ref_if<is_autodiff_v<T_scale>>(y_scaled * y_scaled);

  size_t N = max_size(y_ref, mu_ref, sigma_ref);
  T_partials_return logp = -0.5 * sum(y_scaled_sq);
  if constexpr (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
  }
  // division by size of sigma is a trick to work with vectors and scalars
  if constexpr (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma);
  }

  if constexpr (is_any_autodiff_v<T_y, T_scale, T_loc>) {
    // Evaluate this if it's used more than once
    auto scaled_diff = to_ref_if<
        is_autodiff_v<
            T_y> + is_autodiff_v<T_scale> + is_autodiff_v<T_loc> >= 2>(
        inv_sigma * y_scaled);
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = -scaled_diff;
    }
    if constexpr (is_autodiff_v<T_scale>) {
      partials<2>(ops_partials) = inv_sigma * y_scaled_sq - inv_sigma;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = std::move(scaled_diff);
    }
  }
  return ops_partials.build(logp);
}
```

Some of the same tricks from the above sections are used here in clever ways.
For instance, when calculating `inv_sigma`, if that expression is used for
calculating multiple partials then we want to evaluate it once on that line
and reuse the precomputed operation multiple times in the preceding code. However
if it's not used multiple times then we just want an expression that will then
later be evaluated at its final destination. The same happens for `y_scaled_sq`
and `scaled_diff`.

One odd piece of code here is

```cpp
// division by size of sigma is a trick to work with vectors and scalars
if constexpr (include_summand<propto, T_scale>::value) {
  logp -= sum(log(sigma_val)) * N / math::size(sigma);
}
```

If `sigma_val` is a scalar then we want to decrement the joint log probability
by `log(sigma_val) * N`, but if it's a vector we want to decrement it by
`sum(log(sigma_val))`. In Stan, passing a scalar to `sum()` is a no-op
that just returns back the scalar, so we can use that to have the left hand side
of the multiply work with both vectors and scalars. We always multiply the
returned scalar from `sum()` by `N`, which we don't want if `sigma_val` is a vector.
So then we just divide by the `size()` of `sigma`, which for a scalar will be 1
and for a vector will be `N`.

The same logic can also be written with an explicit `if constexpr` branch

```cpp
if constexpr (include_summand<propto, T_scale>::value) {
  if constexpr (is_vector<T_scale>::value) {
    logp -= sum(log(sigma_val));
  } else {
    logp -= log(sigma_val) * N;
  }
}
```

But that's it, you can see the full `normal_lpdf` function [here](https://github.com/stan-dev/math/blob/develop/stan/math/prim/prob/normal_lpdf.hpp) in Stan that uses the Eigen version.

One other little piece you'll want to do is add a `normal_lpdf` function with the exact same signature
but without the `propto` parameter. Unless told otherwise we don't assume users want the proportional constants
added so we have a default signature that does not require setting the `propto` parameter.

```cpp
template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> normal_lpdf(T_y&& y, T_loc&& mu,
                                                      T_scale&& sigma) {
  return normal_lpdf<false>(std::forward<T_y>(y), std::forward<T_loc>(mu),
                            std::forward<T_scale>(sigma));
}
```


### Testing The Distribution

For testing the new distribution you'll be adding to the [Distribution Testing Framework](https://github.com/stan-dev/math/tree/develop/test/prob). In `test/prob` you will add a new folder with the name
of your distribution with the `.hpp` files for generating the tests inside of it.
Looking at the
[normal distribution](https://github.com/stan-dev/math/tree/develop/test/prob/normal)
testing files you'll see you need

- `mydist_test.hpp` containing a class inheriting from `AgradDistributionTest` containing methods
  - `valid_values()` which fills the `[in/out]` params argument with valued testing values for your distribution and the `[in/out] log_prob` argument with the log probability given those parameters.
  - `invalid_values()` which fills the `[in/out]` value argument with testing values that should fail.
  - Two `log_prob()` methods which calls your distribution, one with `propto` and the other without.
  - `log_prob_function()` which is the log probability density function's simple implementation much like the one you wrote in Stan.


The other test files are quite similar, though note the `normal` cdf test uses some extra numeric tricks
where for example the tests for [gumbel_lcdf](https://github.com/stan-dev/math/blob/develop/test/prob/gumbel/gumbel_cdf_test.hpp) are very similar to the `normal_lpdf` test.

Once you've added those files you can run your tests for your distribution with `./runTests.py` pointing
it to the folder containing your tests. This will generate all of the tests for each of Stan's scalar types
and mixes of scalars and vectors.

```
./runTests.py ./test/prob/mydist/
```

For more details see the [distribution testing docs](@ref dist_tests).

### That's it!

The above should cover all of the `lpdf`, `lcdf`, `cdf`, `lccdf` and `rng`
functions needed to add a new distribution to Stan. If you have further questions
please reach out on our [discourse](https://discourse.mc-stan.org/) or file an
[issue](https://github.com/stan-dev/math/issues)!
