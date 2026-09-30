# Reverse mode (`rev`)

Background: `common_pitfalls.md` and `reverse_mode_types.md` in
`doxygen/contributor_help_pages/`.
See `stan/math/rev/fun/multiply.hpp` for an example of what the code here should look like.

- **Memory.**
  - Everything a callback captures must have memory that only lives in the stack arena and is trivially destructable: `arena_t<T>`, `to_arena(x)`, `arena_matrix`, `make_zeroed_arena`, or a `var`.
  - Arena memory never runs destructors. Capturing a `std::vector`, `Eigen::Matrix` or an Eigen expression by value leaks or dangles.
- **Eigen and `auto`.**
  - Use `.val()`/`.adj()` for value/adjoint views of `var` matrices. Use `.val_op()`/`.adj_op()` when a view is needed inside a matrix product.
- **Two matrix representations.**
  - `Eigen::Matrix<var, ...>` is selected with `require_eigen_vt<is_var, T>`.
  - `var_value<Eigen::Matrix<double, ...>>` is selected with `require_var_matrix_t<T>`.
  - A new rev matrix function usually supports both; test with `expect_ad_matvar`.
  - Convert between them with `to_var_value` / `from_var_value`, and pick the return type with `conditional_var_value_t`.
- **Nesting and drivers.**
  - Inner gradients use a `nested_rev_autodiff` scope. Don't call `start_nested` by hand.
