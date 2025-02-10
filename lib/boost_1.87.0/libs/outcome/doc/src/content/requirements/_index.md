+++
title = "Prerequisites"
weight = 2
+++

Outcome is a header-only C++ 14 library known to work well on the latest
point releases of these compiler-platform combinations or better:

- clang 4.0.1 (LLVM) [FreeBSD, Linux, OS X]
- GCC 6.5 [Linux]
- Visual Studio 2017.9 [Windows]
- XCode 9 [MacOS]

For non-Windows non-POSIX platforms (typically embedded systems), Outcome
is usable in its Outcome.Experimental form with the macro `BOOST_OUTCOME_SYSTEM_ERROR2_NOT_POSIX`
defined.

It is worth turning on C++ 17 or C++ 20 if you can, as there are many usability and
performance improvements. Any Concepts TS or Coroutines TS implemented
by your compiler is automatically detected and used.


Known compiler issues (this was last updated April 2023):

- clang 3.5 - 3.9 can compile varying degrees of the test suite, the
problem is lack of complete and unbuggy C++ 14 language support.

- Older point releases of GCCs 7 and 8 have internal compiler error bugs
in their constexpr implementation which tend to be triggered by using
Outcome in constexpr. If you don't use Outcome in constexpr, you won't
see these problems. If you need your GCC to not ICE, upgrade to the
very latest point release, the constexpr ICE has been since fixed.

- Early editions of Visual Studio 2017 have many corner case problems.
From VS2017.9 onwards there remain a number of usually untroublesome corner
case issues, but use should be relatively unsurprising for most use cases.
Be aware that only from Visual Studio 2022 onwards are almost all corner
case problems fixed.

- Some point releases of GCC 10 with libstdc++ 10 can induce an infinite
template instantiation, which fails the build for some rare use cases. Earlier
or later GCCs or different point releases of the 10 series do not have
this issue.

---

"C++ 14" compilers which do not work, and will not work until their
maintainers fix them:

- GCC 5, due to a bug in nested template variables parsing which was fixed
in GCC 6. I appreciate that this upsets a lot of users. Please raise your
upset at https://gcc.gnu.org/bugzilla/. In the meantime, you can get fairly
far in Outcome with even clang 3.5.
- Any compiler which uses the libstdc++ version which comes with GCC 5, as it does
not implement enough of the C++ 14 standard library for Outcome to compile.
