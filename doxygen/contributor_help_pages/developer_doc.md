# Developer Guide {#developer_guide}

Thanks for checking out these docs. This is where you'll find information on contributing to the [Math library](https://github.com/stan-dev/math), how the source code is laid out, and other technical details that could help. This wiki focuses on things relevant to the Math library. There's also a separate [Stan wiki](https://github.com/stan-dev/stan/wiki) for things related to the language, services, algorithms.

For the most up to date guide on contributing to Stan Math see the [Getting Started Guide](@ref getting_started), [Common Pitfalls](@ref common_pitfalls), and [Adding New Distributions Guide](@ref new_distribution).

This wiki is a work in progress. If you have suggestions, please update the wiki or mention it on our forums on [this thread](https://discourse.mc-stan.org/t/what-doc-would-help-developers-of-the-math-library/).

-----------------

# Contents

- [Overview](#overview)
- [Licensing](#licensing)
- [Contributing](#contributing)
- [Code Review Guidelines](#code-review-guidelines)
- [Building and Running Tests (makefile documentation)](#building-and-running-tests)
- [Questions](#questions)

-----------------

# Overview {#overview}


The Stan Math library, referred to as Math or the Math library, is a C++ library for automatic differentiation. It's designed to be usable, extensive and extensible, efficient, scalable, stable, portable, and redistributable in order to facilitate the construction and utilization of algorithms that utilize derivatives.

The Math library implements:

- [reverse mode](https://en.wikipedia.org/wiki/Automatic_differentiation#Reverse_accumulation) automatic differentiation for computing gradients. This is fully tested and is utilized by [Stan](https://github.com/stan-dev/stan).
- [forward mode](https://en.wikipedia.org/wiki/Automatic_differentiation#Forward_accumulation) automatic differentiation for computing directional derivatives. This is not fully tested, but is close to complete.
- mixed mode automatic differentiation for computing higher order derivatives. Once forward mode is fully tested, this should work.

Some key features of the Math library's reverse mode automatic differentiation:
- object oriented design with overloading of operators
- arena-based memory management

For implementation details of the Math library's automatic differentiation, please read the arXiv paper "[The Stan Math Library: Reverse-Mode Automatic Differentiation in C++](https://arxiv.org/abs/1509.07164)."


# Licensing {#licensing}

We're committed to having a permissive open-source license. The Math library is [licensed with the BSD 3-Clause License](https://github.com/stan-dev/math/blob/develop/LICENSE%2Emd) and we only accept changes to the code base that compatible with this license.

# Contributing {#contributing}

Thanks for reading! We love contributions from everyone in the form of good discussion, issues, and pull requests.

## Issues

We reserve [issues](https://github.com/stan-dev/math/issues) for bugs and feature requests that are defined well enough for a developer to tackle. If you have general questions about the Math library, please see the [Discussion](#discussion) section.

### Bug Reports

Ideally, bug reports will include these pieces of information:

1. a description of the problem
2. a reproducible example
3. the expected outcome if the bug were fixed.

If there's an error and you can produce any of these pieces, it's much appreciated!


### Feature Requests

We track the development of new features using the same [issue tracker](https://github.com/stan-dev/math/issues). Ideally, feature requests will have:

1. a description of the feature
2. an example
3. the expected outcome if the feature existed.

Open feature requests should be the ones we want to implement in the Math library. We'll try to close vague feature requests that don't have enough information and move the discussion to the forums.

## Pull Requests

All changes to the Math library are handled through [pull requests](https://github.com/stan-dev/math/pulls). Each pull request should correspond to an issue. We follow a [modified GitFlow branching model](https://github.com/stan-dev/stan/wiki/Developer-process-overview#2-create-a-branch-for-the-issue) for development.

When a contributor creates a pull request for inclusion to the Math library, here are some of the things we expect:

1. the contribution maintains the Math library's open-source [license](https://github.com/stan-dev/stan/wiki/Stan-Licensing): 3-clause BSD
2. the code base remains stable after merging the pull request; we expect the `develop` branch to always be in a good state
3. the changes are maintainable. In code review, we look at the design of the proposed code. We also expect documentation. It should look like idiomatic C++.
4. the changes are tested. For bugs, we expect at least one test that fails before the patch and is fixed after the patch. For new features, we expect at least one test that shows expected behavior and one test that shows the behavior when there's an error.
5. the changes adhere to the Math library's [C++ standards](https://github.com/stan-dev/stan/wiki/Coding-Style-and-Idioms). Consistency really helps.

Pull requests are code reviewed after they pass our continuous integration tests. We expect all the above before a pull request is merged. We are an open-source project and once code makes it into the repository, it's on the community to maintain.

It is the responsibility of the contributor submitting the pull request that the code meets these requirements. We're open-source. Once the code gets into the code base, the community of developers take ownership of it.


## Discussion

For general questions, please ask on the forums with the ["Developers" tag](http://discourse.mc-stan.org/c/stan-dev).


# Code Review Guidelines {#code-review-guidelines}

(These are also listed on the
[wiki](https://github.com/stan-dev/stan/wiki/Developer-process-overview#code-review-guidelines))

All pull requests must have these things:

1. clear licensing information. Check that the person that owns the copyright is listed.
2. the code base remains stable after the merge. Every pull request should pass all the continuous integration tests. There are a very few exceptions when they don't, but that's rare.
3. the changes proposed are maintainable by the community. We're an open-source project. Once the changes go into the code base, we're all responsible for the code. Please read the code with these things in mind:
    - overall readability. If it's hard to understand the code and it can be made simpler, please mention that.
    - code design. It should either be designed with well-known C++ patterns or there should be a good reason it isn't.
    - documentation. The documentation that goes with the changes should go in this pull request.
4. the changes are tested. Please verify that there's at least one new test that forces the code to execute. Please also verify there's at least one test that shows how to handle errors in the code. Having just these two tests (without full coverage) is usually enough to trap errors in the future. It is also enough to use these to refactor the code when the time comes.
5. the changes adhere to the Math library's [C++ standards](https://github.com/stan-dev/stan/wiki/Coding-Style-and-Idioms). This is a large code base. We want the quality to be consistent so it's easy to navigate and understand other parts of the code base.

Most of the above is subjective, so please use your best judgement.

Pull requests should be narrow in scope. Please don't let inadvertent changes get into the code base.

## Who Can Review a Pull Request

It would help if anyone that's knowledgable in the code base reviews code. This is also for those that aren't active developers, but can read C++ fluently!


## Who Can Merge a Pull Request
Members of the Stan development team that contribute regularly to the Math library are allowed to merge the pull requests.

If all tests pass and the pull request has been reviewed, if no one merges, \@syclik will do so (not on any schedule, but at least once a week).


## How to Review a Pull Request
Thanks! Here are some guidelines on how to review. At a high level, you're making sure the code submitted is indeed what the submitter says they have submitted. If an issue exists and is good, then it should be fixed. There are lots of ways to solve a problem -- we're not always looking for the perfect solution. We're looking for something that can get into the code base without putting a maintenance burden on other developers.

Feel free to review a narrow part of the pull request. By that, I mean one of these things: build process, numerical computation, template metaprogramming. If you're reviewing a small portion of the pull request, please mention it in the comment and only mark it as a **comment** (not as an "Approve"). In most cases, this won't matter, but in the ones that span across different concerns, one of the core developers will coordinate the review and approve.

1. Before starting, review the linked Math issue. That issue should describe what should happen. The pull request will be how that happens.
2. Wait for the continuous integration to complete and pass. If the tests don't pass, the rest is optional.
3. Review the pull request text:
    1. If any of the checklist is missing, please reject the pull request by leaving a comment for the submitter. If there's no response after some time, please close the pull request.
    2. If the checklist is complete, read the description of the pull request. If the description of the fix doesn't seem right, comment on the pull request.
    3. If the description of the tests don't seem like they are enough to keep the technical debt down, comment on the pull request asking for more tests. If you can, ask for specific tests.
4. Review the contents of the pull request. Dig into the code:
    1. Verify that what was described is indeed what's happening in the code.
    2. Spot check the tests to verify that they exist.
    3. Verify documentation exists.
    4. Check that the code is really readable. When you sign off, you're suggesting that someone else will be able to figure out the code to fix or modify it.
    3. Add line-specific comments when you can.
    4. Please check that there aren't other fixes or other code that sneaks into the pull request. Submitters shouldn't be sneaking in patches under the guise of a pull request.
5. "Review Changes" -- hit the green button and submit the review. We're playing within GitHub's rules, so we need to coordinate together. GitHub will allow the pull request to be merged when all continuous integration tests pass, there are no reviews that request changes, and at least one reviewer has approved. Pull requests are tricky because they can actually spawn across different concerns. Here's how we'll coordinate:
    1. **Request changes**. If there are changes that need to be made, finalize your review with "Request changes."
    2. **Comment**.  If you're only reviewing a portion of the pull request, e.g. numeric computation, finalize your review with "Comment" if everything is good.
    3. **Approve**. Finalize the pull request with "Approve" when the pull request is ready to be merged. If you're looking at the whole pull request, then feel free to approve when it's ready. If you're only looking at a part of the pull request, please do not approve. One of the core Math developers will coordinate the different reviewers and approve based on the feedback.


# Building and Running Tests {#building-and-running-tests}

**Nota bene:** these build instructions are not in the released version yet. This is for the `develop` branch (post v2.18.0. Prior to this, the build instructions are similar, but not identical.

## Overview

The Math library is designed as a header-only library. Everything is included in a single translation unit.

Some of the external libraries require libraries to be linked in. These are:

- SUNDIALS for ordinary differential equations solving
- Boost MPI
- OpenCL

Within the Math library, the only build targets are the tests (there are no other executables to build directly from Math). The makefile can be used to build other C++ executables, with all the correct compiler flags for Math.

## Makefile Variables

To customize how the tests are built in C++, variables can be set in a file called `make/local`. This is the preferred way of customizing the build. You can also set it more permanently by setting variables in `~/.config/stan/make.local`.

Note: variables in `make/local` override variables in `~/.config/stan/make/local`.

There are a lot of make variables that can be set. In general, `CXXFLAGS_*` is for C++ compiler flags, `CPPFLAGS_*` is for C preprocessor flags, `LDFLAGS_*` is for linker flags, and `LDLBIS_*` is for libraries that need to be linked in.

These are the more common make flags that could be set:
- `CXX`: C++ compiler
- `CXXFLAGS_OS`: compiler flags specific to the operating system.
- `CPPFLAGS_OS`: C preprocessor flags specific to the operating system.
- `O`: optimization level. Defaults to `3`.
- `INC_FIRST`: this is a C++ compiler option. If you need to include any headers before Math's headers, this is where to specify it. Default is empty.
- `LDFLAGS_OS`: linker flags for the operating system
- `LDLIBS_OS`: link libraries for the operating system


These are the rest of the variables that can be set:

- C++ compiler flags
    - `CXXFLAGS_LANG`: sets the language. Currently defaults to `-std=c++1y`
    - `CXXFLAGS_WARNINGS`: compiler options to squash compiler warnings
    - `CXXFLAGS_BOOST`: Boost-specific compiler flags
    - `CXXFLAGS_EIGEN`: Eigen-specific compiler flags
    - `CXXFLAGS_OPENCL`: OpenCL-specific compiler flags
    - `CXXFLAGS_MPI`: MPI-specific compiler flags
- C preprocessor flags
    - `CPPFLAGS_LANG`:
    - `CPPFLAGS_WARNINGS`
    - `CPPFLAGS_BOOST`
    - `CPPFLAGS_EIGEN`
    - `CPPFLAGS_OPENCL`
    - `CPPFLAGS_MPI`
- Linker flags
    - `LDFLAGS_LANG`
    - `LDFLAGS_WARNINGS`
    - `LDFLAGS_BOOST`
    - `LDFLAGS_EIGEN`
    - `LDFLAGS_OPENCL`
    - `LDFLAGS_MPI`
- Link libraries
    - `LDLIBS_LANG`
    - `LDLIBS_WARNINGS`
    - `LDLIBS_BOOST`
    - `LDLIBS_EIGEN`
    - `LDLIBS_OPENCL`
    - `LDLIBS_MPI`
- Misc
    - `OS`: operating system. Defaults to `Darwin`, `Linux`, or `WindowsNT`.
    - `MATH`: the location of the math library. Defaults to empty variable (used for Stan and CmdStan).
    - `EIGEN`: location of the Eigen headers
    - `BOOST`: location of the Boost headers
    - `SUNDIALS`: location of the Sundials headers
    - `INC`: the includes for the C++ compiler
    - `INC_SUNDIALS`: the Sundials includes
    - `INC_GTEST`: the Google test includes
    - `EXE`: the executable file extension. On Windows, defaults to `.EXE`. On other operating system, defaults to an empty variable.

For debugging, there's a useful target for printing Stan variables. It's `print-*` where `*` is the variable. For example, on a Mac:

```
> make print-CXX
CXX = clang++
```

## Building Tests

The main function of the makefiles is to build test executables. The typical way to build test executables is to use the `runTests.py` python script, but we can use make to directly generate a test executable.

The name of the test executable is the name of the test with the file extension removed for Linux and Mac or replacing the file extension with `.exe` for Windows. For example, to run the test `test/unit/math_include_test.cpp`, we build on Linux and Mac with:

```
> make test/unit/math_include_test
```

or on Windows with:

```
> make test/unit/math_include_test.exe
```

We can then run the executable from the command line.


### Rebuilding Tests

Along with the test executable, the makefiles generate dependency files. The dependency files are generated by using the C++ compiler to list the headers that the test depends on. The dependency file is automatically generated by make and included by make so it knows what header files the test depends on. If none of those files changed, then the executable doesn't need to be rebuilt.

We can build the dependency file directly. For the `test/unit/math_include_test.cpp` example, we can build the dependency file using make:

```
> make test/unit/math_include_test.d
```

If you look at the first couple of lines, you'll see:

```
test/unit/math_include_test.o test/unit/math_include_test.d: \
  test/unit/math_include_test.cpp stan/math.hpp stan/math/rev/mat.hpp \
```

Make uses the timestamp of the files to determine whether it needs to be rebuilt. If any of the files listed after the targets (after the `:`) are updated, the executable will be rebuilt. If you've built all the unit tests and only change a single header file, building the unit tests again will selectively rebuild the tests that depends on that header.

## Running Tests

The easiest way to build and run tests is to use the `runTests.py` python script. To run unit tests:

```
> ./runTests.py test/unit
```

## Misc: Stan Testing

The [Stan testing](https://github.com/stan-dev/stan/wiki/Coding-Style-and-Idioms#unit-testing) process depends on the makefiles in Math.



# Where do I create a new issue

Stan's development is across multiple repositories. When in doubt of which one
to open an issue in, we ask that you make your best guess. We can always move it for you.

Language or compiler issues: stanc3

Algorithm issues: stan

Function or derivative issues: math

# Running tests from the math library

Running tests in the Stan Math Library require these libraries, which are all included in the repository:

- Boost
- Eigen
- Google Test
- CppLint (optional)

No additional configuration is necessary to start running with the default libraries.

If you want to use custom locations for the library locations, set these makefile variables:

- `EIGEN`
- `BOOST`
- `GTEST`
- `CPPLINT` (optional)

Example `~/.config/stan/make.local` file:
```
BOOST = ~/boost
```

# Running tests

To run tests, you will need a copy of the Math library, a C++ compiler, make, and python 2.x (for the test script).

To run the unit tests, type:
```
> ./runTests.py test/unit
```

To run the auto-generated distribution tests, type:
```
> ./runTests.py test/prob
```

To run the multiple translation unit tests, type:
```
> ./runTests.py test/unit/multiple_translation_units_test.cpp
```

If you see this message:

```
------------------------------------------------------------
make generate-tests -s
test/prob/generate_tests.cpp:9:10: fatal error: 'boost/algorithm/string.hpp' file not found
#include <boost/algorithm/string.hpp>
         ^
1 error generated.
make: *** [test/prob/generate_tests] Error 1
make generate-tests -s failed
```

the library paths have not been configured correctly.

To test headers,
```
> make test-headers
```

--------
# Questions {#questions}

This would have been an "faq," but nothing really gets asked frequently. We'll just collect thoughts on the Math library here.

## How do we keep consistent behavior in the Math library?
As of v2.18 (12/2018).

We try to keep the behavior of the functions in the Math library consistent for edge cases. Unfortunately, we implemented a lot of functions before we started realizing we should be consistent. That said, we prefer backwards compatibility to consistency -- we shouldn't change behavior until we get to a major version change. If we throw, we want well-formed messages because they do get passed back to the user from within Stan.

We follow these rules in order:
1. Existing functions should maintain their behavior.
2. Use the standard library's behavior. That means that we should be using the standard library functions for math if it exists.
3. Use Boost's math functions. If the error messages aren't nice, find some way to provide nice error messages.

## Why do we have `BOOST_DISABLE_ASSERTS` included in the C++ build?

Without including this, Boost will assert if certain inputs do not meet the preconditions of the function. Assertions are difficult to trap and recover from and we want to continue to have control over this behavior.

See [Discourse: Boost defines](https://discourse.mc-stan.org/t/boost-defines/10087) for more details.


--------
### Future topics to cover

- what goes in issues, pull requests, and the forums

- how to submit a pull request: requirements

- how to review a pull request

- A description of the map from Stan datatypes to C++ datatypes. (i.e. [stan -> c++] -> [vector[n] y -> Eigen::Matrix<double, Eigen::Dynamic, 1>)

- There is a way to visualize the objects put on the stack (stan::math::print_stack().) A description of this feature of Stan development would be cool, we obviously want the library as efficient as possible. If there’s any other features like this we haven’t seen yet, please let us know!

- More in depth description of all of the memory management features unique to the stan library (i.e. stan::math::recover_memory()).
- Description of how to use/print the gradient functions in Stan, and any other diagnostic tools for making sure our code/math is accurate ( some_var_.grad(vars, grad)) and then (for (auto g : grad) { std::cout << g << std::endl;}

- How does the architecture from rev-> prim work? I threw some print statements in gp_exp_quad_cov formerly (cov_exp_quad) to see how it works. rev fuctions are precompiled, and then call some of the prim functions on runtime? How does this work exactly? What kind of templating it happening to make sure rev is going to which prim? Can we have a quick visualization on a simple stan program to know what’s happening?
- where to find doxygen: https://mc-stan.org/math

- Intro to C++ concepts needed to develop for Stan.
- Description of directory structure of the Stan source
- Guidelines for what should be in the Stan code and what should be in user specific functions.
- Walkthrough of how how to implement, something, required documentation, pull request process etc.
