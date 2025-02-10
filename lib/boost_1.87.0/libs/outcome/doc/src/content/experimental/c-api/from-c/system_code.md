+++
title = "C system error results"
description = "Status code's `std::error` in C"
weight = 10
+++

In v2.2.11, C Result support went from second tier to first tier status, and
now you can create, query and manipulate a subset of Result types entirely from
within C by including `<boost/outcome/experimental/result.h>`.

The subset supported are those `result<T, E>` which are [a `status_result<T>`]({{% relref "../../status_result" %}})
i.e. the `E` is hardcoded to `experimental::error` which is the type erased runtime
polymorphic holder for any errored `status_code` whose payload is not bigger
than an `intptr_t`. This is the most useful subset of Outcome Experimental's
possible Result types, allowing arbitrary custom error coding schemes from
any unknown source to work seamlessly with all others, including errors from
the system or third party libraries.

The operations available to C are:

<dl>
<dt><code>BOOST_OUTCOME_C_DECLARE_RESULT_SYSTEM(ident, T)</code>
<dd>Declares to C a <code>status_result<T></code> type uniquely
identified by <code>ident</code>. <code>T</code> is available at the
member variable <code>.value</code>, and <code>struct cxx_status_code_system</code>
is available at the member variable <code>.error</code>. If in C++,
implements C extern functions for making successful and failure results
of this type. If you call this from within
C++, make SURE it is not within a <code>extern "C"</code> block!

<dt><code>BOOST_OUTCOME_C_RESULT_SYSTEM(ident)</code>
<dd>A reference to a previously declared <code>status_result</code> type with
unique <code>ident</code>.

<dt><code>BOOST_OUTCOME_C_MAKE_RESULT_SYSTEM_SUCCESS(ident, expr)</code> (needs C++ counterpart linked into final binary)
<dd>This invokes the aforementioned extern function which creates a <code>status_result</code>
with a successful value of type <code>T</code>.
<dt><code>BOOST_OUTCOME_C_MAKE_RESULT_SYSTEM_FAILURE_POSIX(ident, expr)</code> (needs C++ counterpart linked into final binary)
<dd>This invokes the aforementioned extern function which creates a <code>status_result</code>
with a failure of type <code>posix_code</code> representing a POSIX <code>errno</code>.
<dt><code>BOOST_OUTCOME_C_MAKE_RESULT_SYSTEM_FAILURE_SYSTEM(ident, expr)</code> (needs C++ counterpart linked into final binary)
<dd>This invokes the aforementioned extern function which creates a <code>status_result</code>
with a failure of type <code>posix_code</code> representing a POSIX <code>errno</code>
if on POSIX; if on Windows then a failure of type <code>win32_code</code>
representing a Win32 error code from a Windows API.

<br><br>
<dt><code>BOOST_OUTCOME_C_RESULT_HAS_VALUE(r)</code>
<dd>Evaluates to 1 (true) if the input <code>result</code> has a value.

<dt><code>BOOST_OUTCOME_C_RESULT_HAS_ERROR(r)</code>
<dd>Evaluates to 1 (true) if the input <code>result</code> has an error.

<dt><code>BOOST_OUTCOME_C_RESULT_ERROR_IS_ERRNO(r)</code>
<dd>Evaluates to 1 (true) if the input <code>result</code>'s error value
is a code in the POSIX <code>errno</code> domain.
<br><br>
<dt><code>BOOST_OUTCOME_C_RESULT_SYSTEM_TRY(expr)</code>
<dd>If the <code>status_result</code> returned by <code>expr</code> is
errored, exit the current function returning the result. This obviously
requires that the return type of the current function matches that of <code>
expr</code>.

<dt><code>BOOST_OUTCOME_C_RESULT_SYSTEM_TRY(cleanup, expr)</code>
<dd>Same as the above, but execute <code>cleanup</code> just before exiting the function
if returning failure.

<dt><code>BOOST_OUTCOME_C_RESULT_SYSTEM_TRY(var, cleanup, expr)</code>
<dd>Same as the above, but set <code>var</code> equal to the result's <code>.value</code> on success.

<dt><code>BOOST_OUTCOME_C_RESULT_SYSTEM_TRY(var, ident, cleanup, expr)</code>
<dd>Same as the above, but use <code>ident</code> as the return type instead. This allows
the return type of the calling function to differ from that of <code>expr</code>.
<br><br>
<dt><code>BOOST_OUTCOME_C_DECLARE_RESULT_SYSTEM_FROM_ENUM(ident, enum_name, uuid, {enum mapping-sequence, ...})</code>
<dd>This declares to C an extern function which creates a <code>status_result</code>
from a C enum. If in C++, it implements a <code>quick_status_code_from_enum</code> for
the C enum and the associated extern function, and you will need to supply <code>uuid</code>
and the appropriate enum value mapping sequence <a href="{{% relref "../../worked-example" %}}">
as per the <code>quick_status_code_from_enum</code> documentation</a>.
<dt><code>BOOST_OUTCOME_C_MAKE_RESULT_SYSTEM_FROM_ENUM(ident, enum_name, expr)</code> (needs C++ counterpart linked into final binary)
<dd>This invokes the aforementioned extern function which creates a <code>status_result</code>
from a C enum.
</dl>

The operations available to C++ are:

<dl>
<dt><code>BOOST_OUTCOME_C_TO_RESULT_SYSTEM_CODE(ident, status_code&lt;T&gt;)</code>
<dd>Returns a previously declared C Result from its matching C++ <code>status_code</code>.
NOTE that the destructor of the C++ status code is NOT called. If this is important
to your status code, it is 100% on you to ensure that your C Result reenters a C++
Result at the end of its lifetime.

<dt><code>to_result(any C Result)</code>
<dd>This is an overloaded C++ free function which returns the C++ status_code&lt;T&gt;
matching its input C Result.
</dl>

Using the above you can write C code using Outcome.Experimental's Result type
quite effectively. Let's look at an example of use next.
