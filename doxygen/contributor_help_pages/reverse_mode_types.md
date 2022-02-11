# Reverse Mode Types {#reverse_mode_types}

### Array of Structures vs. Structure of Array types (AOS vs. SOA, matvar vs. varmat)

There are two types of reverse mode autodiff in Math. They are colloquially
referred to as matvar and varmat, though it is more accurate to
refer to them as Array of Structure (AOS) and Structure of Array (SOA)
autodiff types.

The first autodiff type in Stan was `var`. A `var` is the autodiff version
of a `double`. Constructing matrix autodiff types with `var` is as simple
as building an Eigen matrix and filling it with vars, so
`Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>`. This is an AOS
autodiff style.

The SOA autodiff style defines vector and matrices themselves as autodiff
variables. The SOA matrix is `var_value<Eigen::MatrixXd>`. Defining an
entire matrix as an autodiff type makes it possible to do a number of
optimizations that are not possible when the autodiff types are just
variables stored inside another structure.



### Tests (SOA types)

SOA tests are separate from the regular AOS tests. To test SOA types, add for
every `expect_ad` call an `expect_ad_matvar` call. `expect_ad` builds on the
`prim` tests by using finite derivatives to check all of the higher order
autodiff. `expect_ad_matvar` builds on the `expect_ad` tests by checking the
SOA autodiff agains the AOS autodiff.
