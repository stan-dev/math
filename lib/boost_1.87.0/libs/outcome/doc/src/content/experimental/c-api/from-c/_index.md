+++
title = "C Results"
description = "Outcome's C Result support"
weight = 30
+++

The C macro API header `<boost/outcome/experimental/result.h>` has some macros for working with any kind of Result:

<dl>
<dt><code>BOOST_OUTCOME_C_DECLARE_RESULT(ident, T, E)</code>
<dd>Declares to C a <code>basic_result<T, E></code> type uniquely
identified by <code>ident</code>. <code>T</code> is available at the
member variable <code>.value</code>, and <code>E</code> is available
at the member variable <code>.error</code>. If you call this from within
C++, make SURE it is not within a <code>extern "C"</code> block!

<dt><code>BOOST_OUTCOME_C_RESULT(ident)</code>
<dd>A reference to a previously declared <code>result</code> type with
unique <code>ident</code>.

<dt><code>BOOST_OUTCOME_C_RESULT_HAS_VALUE(r)</code>
<dd>Evaluates to 1 (true) if the input <code>result</code> has a value.

<dt><code>BOOST_OUTCOME_C_RESULT_HAS_ERROR(r)</code>
<dd>Evaluates to 1 (true) if the input <code>result</code> has an error.

<dt><code>BOOST_OUTCOME_C_RESULT_ERROR_IS_ERRNO(r)</code>
<dd>Evaluates to 1 (true) if the input <code>result</code>'s error value
is a code in the POSIX <code>errno</code> domain.
</dl>

The above let you work, somewhat awkwardly, with any C-compatible
`basic_result<T, E>`. `basic_result<T, E>` is trivially copyable and
standard layout if its `T` and `E` are both so, and it has the C layout:

```c++
struct cxx_result_##ident
{
  union
  {
    T value;
    E error;
  };
  unsigned flags;
};
```

Note that this layout is different to that of [`BOOST_OUTCOME_C_DECLARE_STATUS_CODE`]({{% relref "../from-c" %}})
as the C++ `result` has a different layout if `E` is a status code.


