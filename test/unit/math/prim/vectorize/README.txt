Instructions for implementing and testing a vectorized function
----------------------------------------------------------------

* Implement the function itself with vectorization in a single file.

  For an example, see: test/unit/math/prim/mat/vectorize/foo_fun.hpp

    *  Define a struct definining the function for scalars which works
       for all of the scalar types (prim, rev, fwd, and mix);
       see struct foo_fun in foo_fun.hpp

    *  Define the actual function in terms of the struct.
       see definition of foo() in foo_fun.hpp

* For testing, we need to define a test class that implements
  the function independently and supplies sequences of legal and
  illegal double and int inputs.

  For an example, see the definition of struct foo_base_test in
    test/unit/math/prim/mat/vectorize/foo_base_test.hpp

  * Need all of the functions to be defined

      * T appply(const T&);

      * double apply_base(int)

*** WARNING ***
The tests presuppose that all of the scalar possibilities ***have already been tested***.  These tests only test that nothing gets messed up by the vectorization.
