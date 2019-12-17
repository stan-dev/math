+++
title = "Changelog"
weight = 80
+++

---
## v2.1.2 11th December 2019 (Boost 1.72) [[release]](https://github.com/ned14/outcome/releases/tag/v2.1.2)

### Enhancements:

Improved compatibility with cmake tooling
: Standalone outcome is now `make install`-able, and cmake `find_package()` can find it.
Note that you must separately install and `find_package()` Outcome's dependency, quickcpplib,
else `find_package()` of Outcome will fail.

Non-permissive parsing is now default in Visual Studio
: The default targets in standalone Outcome's cmake now enable non-permissive parsing.
This was required partially because VS2019 16.3's quite buggy Concepts implementation is
unusuable in permissive parsing mode. Even then, lazy ADL two phase lookup is broken
in VS2019 16.3 with `/std:latest`, you may wish to use an earlier language standard.

**Breaking change!**
: The git submodule mechanism used by standalone Outcome of specifying dependent libraries
has been replaced with a cmake superbuild of dependencies mechanism instead. Upon cmake
configure, an internal copy of quickcpplib will be git cloned, built and installed into the
build directory from where an internal `find_package()` uses it. This breaks the use of
the unconfigured Outcome repo as an implementation of Outcome, one must now do one of:
 1. Add Outcome as subdirectory to cmake build.
 2. Use cmake superbuild (i.e. `ExternalProject_Add()`) to build and install Outcome into
 a local installation.
 3. Use one of the single header editions.

**Breaking change!**
: For standalone Outcome, the current compiler is now checked for whether it will compile
code containing C++ Concepts, and if it does, all cmake consumers of Outcome will enable
C++ Concepts. Set the cmake variable `BOOST_OUTCOME_C_CONCEPTS_FLAGS` to an empty string to prevent
auto detection and enabling of C++ Concepts support occurring.

`BOOST_OUTCOME_TRY` operation now hints to the compiler that operation will be successful
: [P1886 *Error speed benchmarking*](https://wg21.link/P1886) showed that there is
considerable gain in very small functions by hinting to the compiler whether the expression
is expected to be successful or not. `BOOST_OUTCOME_TRY` previously did not hint to the compiler
at all, but now it does. A new suite of macros `BOOST_OUTCOME_TRY_FAILURE_LIKELY` hint to the
compiler that failure is expected. If you wish to return to the previously unhinted
behaviour, define `BOOST_OUTCOME_TRY_LIKELY(expr)` to `(!!expr)`.

[#199](https://github.com/ned14/outcome/issues/199)
: Support for C++ Coroutines has been added. This comes in two parts, firstly there is
now an `BOOST_OUTCOME_CO_TRY()` operation suitable for performing the `TRY` operation from
within a C++ Coroutine. Secondly, in the header `outcome/coroutine_support.hpp` there are
implementations of `eager<OutcomeType>` and `lazy<OutcomeType>` which let you more
naturally and efficiently use `basic_result` or `basic_outcome` from within C++
Coroutines -- specifically, if the result or outcome will construct from an exception
pointer, exceptions thrown in the coroutine return an errored or excepted result with
the thrown exception instead of throwing the exception through the coroutine machinery
(which in current compilers, has a high likelihood of blowing up the program). Both
`eager<T>` and `lazy<T>` can accept any `T` as well. Both have been tested and found
working on VS2019 and clang 9.

[#210](https://github.com/ned14/outcome/issues/210)
: `make_error_code()` and `make_exception_ptr()` are now additionally considered for
compatible copy and move conversions for `basic_result<>`. This lets you construct
a `basic_result<T, E>` into a `basic_result<T, error_code>`, where `E` is a
custom type which has implemented the ADL discovered free function
`error_code make_error_code(E)`, but is otherwise unrelated to `error_code`.
The same availability applies for `exception_ptr` with `make_exception_ptr()` being
the ADL discovered free function. `basic_outcome<>` has less support for this than
`basic_result<>` in order to keep constructor count down, but it will accept via
this mechanism conversions from `basic_result<>` and `failure_type<>`.

### Bug fixes:

[#184](https://github.com/ned14/outcome/issues/207)
: The detection of `[[nodiscard]]` support in the compiler was very mildly broken.

---
## v2.1.1 19th August 2019 (Boost 1.71) [[release]](https://github.com/ned14/outcome/releases/tag/v2.1.1)

### Enhancements:

[#184](https://github.com/ned14/outcome/issues/184)
: As per request from Boost release managers, relocated `version.hpp` and
`revision.hpp` into detail, and added the Boost licence boilerplate to the top
of every source file which was missing one (I think). Also took the opportunity
to run the licence restamping script over all Outcome, so copyright dates are now
up to date.

[#185](https://github.com/ned14/outcome/issues/185)
: Add FAQ item explaining issue #185, and why we will do nothing to
fix it right now.

[#189](https://github.com/ned14/outcome/issues/189)
: Refactored the `BOOST_OUTCOME_TRY` implementation to use more clarified
customisation points capable of accepting very foreign inputs. Removed the
`std::experimental::expected<T, E>` specialisations, as those are no longer
necessary. Fixed the documentation for the customisation points which
previously claimed that they are ADL discovered, which they are not. Added
a recipe describing how to add in support for foreign input types.

[#183](https://github.com/ned14/outcome/issues/183)
: Added a separate `motivation/plug_error_code` specifically for Boost.

### Bug fixes:

-
: `BOOST_OUTCOME_VERSION_MINOR` hadn't been updated to 1.

[#181](https://github.com/ned14/outcome/issues/181)
: Fix issue #181 where Outcome didn't actually implement the strong swap guarantee,
despite being documented as doing so.

[#190](https://github.com/ned14/outcome/issues/190)
: Fix issue #190 in Boost edition where unit test suite was not runnable from
the Boost release distro.

[#182](https://github.com/ned14/outcome/issues/182)
: Fix issue #182 where `trait::is_exception_ptr_available<T>` was always true,
thus causing much weirdness, like not printing diagnostics and trying to feed
everything to `make_exception_ptr()`.

[#194](https://github.com/ned14/outcome/issues/192)
: Fix issue #192 where the `std::basic_outcome_failure_exception_from_error()`
was being defined twice for translation units which combine standalone and
Boost Outcome's.

---
## v2.1 12th Apr 2019 (Boost 1.70) [[release]](https://github.com/ned14/outcome/releases/tag/v2.1)

- [#180](https://github.com/ned14/outcome/issues/180)
    - `success()` and `failure()` now produce types marked `[[nodiscard]]`.

- `include/outcome/outcome.natvis` is now namespace permuted like the rest of
Outcome, so debugging Outcome based code in Visual Studio should look much
prettier than before.

- [#162](https://github.com/ned14/outcome/issues/162)
    - `.has_failure()` was returning false at times when it should have returned true.

- [#152](https://github.com/ned14/outcome/issues/152)
    - GCC 5 no longer can compile Outcome at all due to [https://stackoverflow.com/questions/45607450/gcc5-nested-variable-template-is-not-a-function-template](https://stackoverflow.com/questions/45607450/gcc5-nested-variable-template-is-not-a-function-template).
Added explicit version trap for GCC 5 to say it can not work. Note this is not a
breaking change, GCC 5 was never supported officially in any v2 Outcome.

- [#150](https://github.com/ned14/outcome/issues/150)
    - **BREAKING CHANGE** `result<T, E>`, `boost_result<T, E>` and `std_result<T, E>`
no longer implement hard UB on fetching a value from a valueless instance if `E` is
a UDT, they now fail to compile with a useful error message. If you wish hard UB,
use `unchecked<T, E>`, `boost_unchecked<T, E>` or `std_unchecked<T, E>` instead.

- [#140](https://github.com/ned14/outcome/issues/140)
    - Fixed a nasty corner case bug where value type's without a copy constructor
but with a move constructor would indicate via traits that copy construction
was available. Thanks to Microsoft's compiler team for reporting this issue.

- Added experimental `status_result` and `status_outcome` based on experimental
`status_code`.

- Boost edition is now 100% Boost, so defaults for `result` and `outcome` are
`boost::system::error_code::errc_t` and `boost::exception_ptr`. Moreover,
the test suite in the Boost edition now exclusively tests the Boost edition.
One can, of course, freely use the standalone edition with Boost, and the Boost
edition with `std` types.

- Renamed ADL discovered customisation point `throw_as_system_error_with_payload()`
to `outcome_throw_as_system_error_with_payload()`.

- [#135](https://github.com/ned14/outcome/issues/135)
    - Added much clearer compile failure when user tries `result<T, T>` or `outcome`
    where two or more types are identical. Thanks to Andrzej Krzemieński
    for suggesting a technique which combines SFINAE correctness with
    the remaining ability for `result<T, T>` etc to be a valid type, but
    not constructible.

- [#67](https://github.com/ned14/outcome/issues/67)
    - Fixed one of the oldest long open bugs in Outcome, that the noexcept
unit tests failed on OS X for an unknown reason.

- [#115](https://github.com/ned14/outcome/issues/115)
    - Outcome did not construct correctly from `failure_type`.

- Inexplicably outcome's error + exception constructor had been removed.
Nobody noticed during the Boost peer review, which is worrying seeing as that
constructor is needed for one of the main advertised features to Boost!

- [#107](https://github.com/ned14/outcome/issues/107) and [#116](https://github.com/ned14/outcome/issues/116)
    - `operator==` and `operator!=` now become disabled if the value, error and
    exception types do not implement the same operator.
    - Relatedly, both comparison operators simply didn't work right. Fixed.

- [#109](https://github.com/ned14/outcome/issues/109)
    - `swap()` now has correct `noexcept` calculation and now correctly orders
    the swaps to be whichever is the throwing swap first.

- Added reference dump of v2.1 ABI so we can check if ABI breakage detection
works in the next set of changes, plus Travis job to check ABI and API compatibility
per commit.

- [#124](https://github.com/ned14/outcome/issues/124)
    - `BOOST_OUTCOME_TRY` is now overloaded and selects `void` or `auto` edition
    according to input parameter count.

- [#120](https://github.com/ned14/outcome/issues/120)
    - Fix generation of double underscored temporary variables in
    `BOOST_OUTCOME_UNIQUE_NAME`, which is UB.

- [#110](https://github.com/ned14/outcome/issues/110)
    - Separated `result` from its hard coded dependency on the `<system_error>` header.
    - Renamed `result` and `outcome` to `basic_result` and `basic_outcome`.
    - Renamed `result.hpp` into `basic_result.hpp`.
    - Moved `<system_error>` and `<exception>` dependent code into new
    `std_result.hpp` and `std_outcome.hpp` header files.
    - Added `boost_result.hpp` and `boost_outcome.hpp` which use Boost.System
    and Boost.Exception (these are `result.hpp` and `outcome.hpp` in the Boost edition).

---
## v2.0 18th Jan 2018 [[release]](https://github.com/ned14/outcome/releases/tag/v2.0-boost-peer-review)

- Boost peer review edition. This is what was reviewed.
- Changelog from v1 can be found in the release notes for this release.
