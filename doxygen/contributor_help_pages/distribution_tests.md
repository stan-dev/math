# Distribution Tests {#dist_tests}

We've written a framework for testing our univariate distribution functions. It's not well doc'ed and the framework is pretty complicated, but it is used to make sure our distributions are ok. (Well -- ok enough.)

This is an attempt to document how to use the framework and how we ended up with this framework.

## How

In this section, I'll describe how to run the tests and what's happening.

### How to Run All the Distribution Tests (for the impatient)

From within the Math library:

```
./runTests.py test/prob

```

For development on windows, add `GENERATE_DISTRIBUTION_TESTS=true` to the make/local file.

This will take hours, perhaps many hours. Most of the time is spent compiling. You might want to add the parallel flag to the python script.


### Steps Involved in Running All the Tests

Here, I'm just going to describe the steps that are taken to run the tests. Details on how to write a test and what's inside the framework are further down. These are the steps taken when calling `./runTests.py test/prob` -- just broken out step by step.

1. `make generate-tests` is called.
    1. This first builds the `test/prob/generate_tests` executable from `test/prob/generate_tests.cpp`.
    2. For each test file inside `test/prob/*/*`, it will call the executable with the test file as the first argument and the number of template instantiations per file within the second argument. For example, for testing the `bernoulli_log()` function, make will call: `test/prob/generate_tests test/prob/bernoulli/bernoulli_test.hpp 100`
    The call to the executable will generate 5 different test files, all within the `test/prob/bernoulli/` folder:
        - bernoulli\_00000\_generated\_fd\_test.cpp
        - bernoulli\_00000\_generated\_ffd\_test.cpp
        - bernoulli\_00000\_generated\_ffv\_test.cpp
        - bernoulli\_00000\_generated\_fv\_test.cpp
        - bernoulli\_00000\_generated\_v\_test.cpp
    These generated files are actual unit tests and can be run just like any other test in the math library. The `{fd, ffd, ffv, fv, v}` are the types of instantiations within the file. Those types map to: `fvar<double>`, `fvar<fvar<double>>`, `fvar<fvar<var>>`, `fvar<var>`, `var`.
    In disributions with many arguments, there will be more than 1 file per types of instantiations.
2. Each of the test files within `test/prob/*/*` are compiled.
3. Each of the test files within `test/prob/*/*` is run.


### How to Run a Single Distribution Test?

The easiest way is to:

1. Generate all the tests:
    `make generate-tests`
2. Run the test. For example, it were the bernoulli test for reverse mode:
    `./runTests.py test/prob/bernoulli/bernoulli_00000_generated_v_test.cpp`

## Why?

Before getting into how to write a distribution test, I'll walk through some of the reasons why we built a framework for testing distribution functions.

As with almost everything complicated in Stan, the testing framework was built in part out of necessity and part as a reaction. The framework was built after we started vectorizing functions. In earlier versions of Stan, we didn't have vectorized functions. For a 3 argument function, like `normal_log(T1 y, T2 mu, T3 sigma)`, we had 8 (2^3; each template argument can be `double` or `stan::nath::var`) instantiations that we wrote out by hand. That's not individual tests. At that time, we probably had 2 tests per instantiation: one for good arguments which would have tested the result, the gradients, and the propto flag, and a second test for invalid arguments which would have tested for exceptions. So, with 8 instantiations, that was 16 hand-written tests. That was managagble, but when we started to vectorize, we realized that it wasn't. With vectorization, that same function now allows for 512 different instantiations. We now are instantiating with different containers. For `double`, it's: `double`, `std::vector<double>`, `Eigen::Matrix<double, -1, 1>`, `Eigen::Matrix<double, 1, -1>`. And we have the instantiations with `double` replaced with `stan::math::var`.

I mentioned a part of this framework was created as a reaction. As we started vectorizing the distributions efficiently, we had buggy gradients in just some of the instantiations. Early on, we were also not as good with templated C++, so not every instantiation compiled. As a precaution, we started testing every instantiation and testing gradients across all of the possible instantiations. In the long run, it's actaully saved us quite a bit of time.


### Original Goals

Some of the original goals were to:

- make sure every instantiation (from the Stan language) was compilable by the C++ compiler
- verify that the distribution functions were calculating the correct values
- verify that the gradients were correct for each of the instantiations; this wasn't hard, but tricky to get right and would have been error prone if it had to be repeated
- verify that invalid arguments were caught
- verify that propto worked across the different instantiations
- verify that passing 0-length vectors worked as expected
- make it easy to test new distributions
- keep using Google test
- run tests on Windows

We ended up with our current framework mostly due to dealing with these constraints.


### Updated Goals

Recently, we added a few more goals:

- test cdfs, ccdfs, cdf_logs, etc.
- test forward mode instantiations

This lead to an expansion of the original framework, which may have added more complication than necessary.


## Writing a Distribution Test Header

Writing a single distribution test isn't overly difficult. It's a bit tedious, but it beats writing every possible instantiation out by hand.

I'll start out by outlining a distribution test. To write a test for cdfs, it's similar. The testing framework file that is associated with the distribution is located at: `test/prob/test_fixture_distr.hpp`. If you want to know the details of the tests, look in here.

Before diving into details, the overall structure of the test file is:

- Annotation line indicating the types of arguments.
- A class definition that subclasses `AgradDistributionTest` (the name is a relic of old-Stan).
- Implementation of five methods in the class:
    - `void valid_values(vector<vector<double> >& parameters, vector<double>& log_prob)`
    - `void invalid_values(vector<size_t>& index, vector<double>& value)`
    -
```
  template <class T_n, class T_prob, typename T2,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_n, T_prob>
  log_prob(const T_n& n, const T_prob& theta, const T2&,
           const T3&, const T4&, const T5&)
```
    -
```
  template <bool propto,
            class T_n, class T_prob, typename T2,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_n, T_prob>
  log_prob(const T_n& n, const T_prob& theta, const T2&,
           const T3&, const T4&, const T5&)
```
    -
```
  template <class T_n, class T_prob, typename T2,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_n, T_prob>
  log_prob_function(const T_n& n, const T_prob& theta, const T2&,
                    const T3&, const T4&, const T5&)
```
- optional: additional tests.



Details:

1. Start by creating a test file. This will be called `distribution_test.hpp` where you'll name the file after the distribution. Our convention is to place this file here: `test/prob/distribution/distribution_test.hpp` where `distribution` is replaced with the distribution name.
2. The first line of this file is an annotation that indicates how to generate the individual test files. Here's an example from the bernoulli test:
```
// Arguments: Ints, Doubles
```
    That's a comment line that starts with `//`, then a space, then the keyword `Arguments`, a colon `:`, then a comma-then-space separated list of argument types. Each of the distribution tests should start with exactly:
```
// Arguments:
```
    The types of arguments that are valid are: `Int`, `Double`, `Ints`, and `Doubles`. The first two indicate that the argument is not vectorized and only takes single, scalar values. The last two plural versions indicate that the argument types are vectorized.
3. The next lines are includes. Most of the files will have:
```
\#include <stan/math/prim.hpp>
\#include <stdexcept>

```
4.




## The Testing Framework

### test_fixture_distr.hpp

### test_fixture_cdf.hpp

### test_fixture_cdf_log.hpp

### test_fixture_ccdf_log.hpp
