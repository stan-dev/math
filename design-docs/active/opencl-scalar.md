# GPU-resident `stan::math::opencl::ScalarCl<T>`

Status: implemented on `feature/opencl-scalar`; see Implementation notes

## Goal

Keep scalar results of OpenCL computations (reductions, log densities,
their adjoints) on the device. Values and adjoints live in 1x1 device
buffers; autodiff nodes and callbacks stay on the CPU. Every host/device
crossing is explicit and visible in the source.

## Decisions

| Topic | Decision |
|-|-|
| Name / namespace | `stan::math::opencl::ScalarCl<double>` and `ScalarCl<var>`. No alias or `using` in `stan::math`. |
| Host conversion | None implicit. Only `opencl::to_host(s)` (blocking) and `explicit ScalarCl(double)` / `explicit ScalarCl<var>(var)`. Accidental crossings must fail to compile. |
| Scalar-only functions | `exp`, `log`, `+`, ... on `ScalarCl` live only in `stan::math::opencl` (reachable by ADL or `opencl::` qualification). No `using` declarations into `stan::math`. |
| Return types | Every OpenCL function that returned a scalar `var` / `double` returns `ScalarCl<var>` / `ScalarCl<double>`: `sum`, `dot_product`, `dot_self`, `log_sum_exp`, `prod`, `trace`, `variance`, `sd`, `squared_distance`, and all OpenCL distributions (`*_lpdf`, `*_lpmf`, `*_cdf`, `*_lcdf`, `*_lccdf`, GLMs, `log_mix`, ...). Integer sums keep a host `int`. |
| Stanc compatibility | Breakage of generated code (e.g. `lp_accum__.add(normal_lpdf(y_cl, ...))`) is accepted until compiler work lands. No accumulator shims. |
| Scalar arguments | Everywhere a host scalar is accepted next to OpenCL arguments, `ScalarCl<double>` / `ScalarCl<var>` is accepted too: kernel generator operands, function arguments (`rep_matrix`, `multiply`, `divide`, `add_diag`, `pow`, `fmax`, ..., constraint bounds), distribution parameters, and handwritten kernels. |
| Host literals | `m * 2.0` keeps the existing by-value `scalar_` path; no upload. |
| Existing CPU-var x GPU ops | Unchanged when given a plain CPU `var` (e.g. `multiply(var, var_value<matrix_cl>)`). |
| Constraint `lp` | Add `ScalarCl<var>& lp` overloads; keep the `var& lp` overloads. |
| GP covariance | `ScalarCl<double>` hyperparameters only (pointer-based kernel variants). `ScalarCl<var>` hyperparameters need new derivative kernels: separate PR. |
| `-cl-std=CL1.2` | Test-only compile check of the new kernel sources. Not added globally. |
| Transfer tracing | Not built. The no-implicit-conversion rule makes crossings compile-time visible instead. |

## Design

### `ScalarCl<double>` (`stan/math/opencl/scalar_cl.hpp`)

- Composition over a fixed 1x1 `matrix_cl<double>`. Must not satisfy
  `is_matrix_cl` or `is_kernel_expression_and_not_scalar` (so it never
  picks up the matrix `var_value` specialization, `arena_type`,
  `accumulator::add`, `value_type` forwarding, ...).
- Default ctor: device zero. `explicit ScalarCl(double)`: must upload through
  the rvalue path of `matrix_cl`'s scalar ctor; the lvalue path enqueues a
  non-blocking write from host memory that may be dead by the time it runs.
- Copy duplicates on device; move transfers ownership.
- Duck-typed buffer/event API (`buffer()`, `write_events()`,
  `read_write_events()`, `add_*_event`) delegating to the backing matrix.
- Trait `stan::is_scalar_cl<T>` in `stan/math/prim/meta/is_scalar_cl.hpp`
  (unguarded, forward declaration only), added to the scalar-inclusive
  `is_kernel_expression` only.

### Kernel generator

- New leaf op `scalar_buf_` (`kernel_generator/scalar_buf.hpp`): reads
  `buf[0]`, `rows()/cols()` = dynamic, so it broadcasts like `scalar_`.
  One kernel arg (the buffer). As an LHS it throws on non-1x1 instead of
  reallocating (unlike `load_::check_assign_dimensions`).
- `as_operation_cl` overloads map `ScalarCl` to `scalar_buf_`.
- Mixed scalar/matrix expressions stay lazy with no changes to the
  binary-operation macros.
- Scalar-only expressions are materialized eagerly by one evaluator that
  reuses `multi_result_kernel` code generation and launches a 1x1 grid.
  Assigning a matrix-sized expression to a `ScalarCl` throws.
- Scalar-only element-wise functions: one macro line per function in
  `stan::math::opencl`, wrapping the existing `fun_<>` op classes.

### Handwritten kernels

- `kernel_cl.hpp`: `ScalarCl` overloads of `get_kernel_args`,
  `select_events`, `assign_event` so it can be passed with
  `in_buffer` / `out_buffer` / `in_out_buffer` as `__global double*`.
- GP covariance kernels (`gp_exp_quad_cov`, `gp_exponential_cov`,
  `gp_matern32_cov`, `gp_matern52_cov`) gain pointer-parameter variants
  for `ScalarCl<double>` hyperparameters. `fill_strict_tri` and
  `check_symmetric` take internal constants and are unchanged.

### Reduction

- `sum(expr)`: `sum_2d` partial sums, then a new single-workgroup kernel
  (strided reads, local-memory tree, one writer) into element 0.
  Variants: overwrite and accumulate, plus an optional host-constant
  offset so distributions fold terms like `NEG_LOG_SQRT_TWO_PI * N` into
  the reduction. Empty input gives device zero. The small-input CPU
  shortcut is removed for floating-point sums. No atomics.

### `ScalarCl<var>` (`stan/math/opencl/rev/scalar_cl.hpp`)

- Wraps a 1x1 `var_value<matrix_cl<double>>`. Copies share identity.
  Assignment and compound assignment rebind to new nodes; never mutate a
  value a callback may have captured.
- `.val()` read-only / `.adj()` writable device views sharing events.
  `value_of`, `adjoint_of` overloads with no transfer.
- CPU `var` bridges (one shared helper): `ScalarCl<var>(cpu_var)` uploads
  the value and adds the device adjoint to the CPU adjoint in the reverse
  pass; `to_host(ScalarCl<var>)` reads the value and adds the CPU adjoint
  into the device adjoint (by-value kernel arg, no readback).

### Autodiff integration

- `partials_propagator` specialization when the return type is
  `ScalarCl`: `build(ScalarCl<double>)` returns `ScalarCl<var>` (or the
  value if nothing is autodiff). Scalar edges hold device partials; a
  matrix partial assigned to a scalar edge is reduced on device. Reverse
  pass: matrix vars `x.adj() += res_adj * partial` (broadcast),
  `ScalarCl<var>` device multiply-add, CPU `var` explicit bridge.
- `adjoint_results_cl`: `ScalarCl<var>` results are reduced into the
  device adjoint (accumulate variant); the CPU `var` path is unchanged.
- `opencl_return_t<Ops...>` = `ScalarCl<var>` if any op is autodiff, else
  `ScalarCl<double>`. `return_type_t` is untouched.

## Commit plan (TDD per commit)

1. Traits (`is_scalar_cl`, `is_kernel_expression`).
2. `ScalarCl<double>`, `to_host`.
3. `scalar_buf_`, `as_operation_cl`, 1x1 evaluator, scalar-only functions.
4. `kernel_cl` hooks, reduction kernel, `-cl-std=CL1.2` compile test.
5. `ScalarCl<var>`, bridges, rebind tests.
6. `sum` family (9 functions, prim + rev).
7. `partials_propagator`, `adjoint_results`; `normal_lpdf` and
   `normal_id_glm_lpdf` end to end.
8. Remaining distributions, batched by family.
9. Scalar-mixed functions (`multiply`, `divide`, `add_diag` fix, `pow`,
   `fmax`, `fmin`, `fdim`, `hypot`, `rep_matrix`, ...).
10. Constraints (`ScalarCl<var>& lp` overloads).
11. GP covariance (`ScalarCl<double>` only).

Test helpers in `test/unit/math/opencl/util.hpp` learn `ScalarCl` before
commit 6; automatic `ScalarCl` argument variants are added once the
distributions migrate.

## Validation

- Focused tests per commit; OpenCL suite with
  `STAN_OPENCL=true python3 runTests.py test/unit/math/opencl`.
- `make test-math-dependencies`, header checks for touched headers, and a
  build without `STAN_OPENCL`.
- CPU-matching values and gradients; generated -> handwritten ->
  generated chains ordered through events.

## Out of scope

- Stanc changes.
- `ScalarCl<var>` GP hyperparameters (new derivative kernels).
- Transfer-count tracing.

## Implementation notes

### Types added during implementation

- Views `ScalarCl<const double&>` (read-only) and `ScalarCl<double&>`
  (writable) are what `ScalarCl<var>::val()` / `.adj()` return. They share
  the backing buffer and events and count as primitive device scalars
  (`is_prim_scalar_cl`), so they broadcast in expressions and can be passed
  to kernels.
- `ScalarCl<arena_matrix_cl<double>>` is `arena_t` of primitive device
  scalars. Reverse-pass callbacks created with `reverse_pass_callback` are
  never destroyed, so captured device scalars must live in the arena.

### Mechanisms

- `opencl::internal::as_operand(x)` turns a device scalar into its
  broadcasting `scalar_buf_` operation (other values pass through). Code that
  applies comparisons, checks or arithmetic to parameters that may be device
  scalars uses it so those pieces fuse into the surrounding kernel. Do not
  bind its result to a local that outlives a returned expression.
- `opencl::internal::as_kernel_buffer(x)` passes a device scalar's own buffer
  to a handwritten kernel and uploads host values.
- `opencl::make_partials_propagator` (explicitly qualified) builds the device
  propagator; `stan::math::make_partials_propagator` rejects OpenCL
  operands.
- `adjoint_results` reduces matrix-sized derivatives into `ScalarCl<var>`
  adjoints on the device.
- Handwritten kernels with scalar parameters are written once with
  `SCALAR_PARAM(name)` / `SCALAR_VALUE(name)` and compiled with the
  `scalar_params_by_value` or `scalar_params_buffer` prefix (GP kernels).
- Mixed operators `+ - * /` between device scalars and matrices, when an
  operand is autodiff, live in `stan::math::opencl` and forward to `add`,
  `subtract`, `multiply`, `divide`, `elt_divide`.

### Fixes to existing code found on the way

- Element-wise kernel generator operations on scalar-only operands reported
  impossible diagonals (`{-rows() + 1, cols() - 1}` with dynamic size), which
  marked results as diagonal-only (`dense_extreme_diagonals()`).
- `results() = expressions()` took the thread count from the first expression
  even if its size was dynamic.
- Sums of expressions containing reductions are evaluated first.
- `add_diag`'s reverse pass doubled the existing adjoint of a scalar `var`.
- OpenCL distributions returned the normalizing constant instead of 0 for an
  empty random variable with scalar parameters (`N == 0` vs the CPU's
  `size_zero`); GLMs keep `N == 0` (instances).
- `promote_scalar_t<double, matrix_cl<double>>` was `double`, so the
  `rep_*` OpenCL tests never ran on OpenCL.

### Remaining explicit host crossings

All are explicit (`opencl::to_host`, `from_matrix_cl` or an explicit
`ScalarCl<double>(x)` upload):

- CPU `var` operands: `ScalarCl<var>(var)` and CPU `var` edges of the
  partials propagator read one adjoint back in the reverse pass.
- Data-dependent early returns in ~27 distributions
  (`if (from_matrix_cl(flag).any()) return ...`) and the lazy argument checks
  in the GLMs (`if (!isfinite(to_host(logp)))`).
- Host-scalar branches of `bernoulli_lpmf` / `binomial_lpmf` compute on the
  host and upload the result.
- Norm check in prim `unit_vector_constrain`, probability sum check in
  `dirichlet_lpdf`, `sigma` checks in `gp_dot_prod_cov` with device scalars.
- `var& lp` / `double& lp` constraint overloads and `accumulator::add`.
- Every `check_cl` reads its flag back (pre-existing).

### Not done

- `ScalarCl<var>` GP hyperparameters (new derivative kernels, separate PR).
- Scalar-only calls of distributions (all arguments device scalars) are not
  supported; at least one argument must be an OpenCL matrix.
- Stanc-generated code using OpenCL distributions or reductions must convert
  results explicitly until the compiler is updated.
