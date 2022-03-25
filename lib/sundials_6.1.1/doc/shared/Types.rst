.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _Usage.CC.DataTypes:

Data Types
----------

The header file ``sundials_types.h`` contains the definition of the types:

* :c:type:`realtype` -- the floating-point type used by the SUNDIALS packages

* :c:type:`sunindextype` -- the integer type used for vector and matrix indices

* :c:type:`booleantype` -- the type used for logic operations within SUNDIALS

Floating point types
~~~~~~~~~~~~~~~~~~~~

.. c:type:: realtype

   The type ``realtype`` can be ``float``, ``double``, or ``long double``, with
   the default being ``double``. The user can change the precision of the
   arithmetic used in the SUNDIALS solvers at the configuration stage (see
   :cmakeop:`SUNDIALS_PRECISION`).

Additionally, based on the current precision, ``sundials_types.h`` defines
``BIG_REAL`` to be the largest value representable as a ``realtype``,
``SMALL_REAL`` to be the smallest value representable as a ``realtype``, and
``UNIT_ROUNDOFF`` to be the difference between :math:`1.0` and the minimum
``realtype`` greater than :math:`1.0`.

Within SUNDIALS, real constants are set by way of a macro called ``RCONST``. It
is this macro that needs the ability to branch on the definition of
``realtype``. In ANSI C, a floating-point constant with no suffix is stored as a
``double``. Placing the suffix "``F``" at the end of a floating point constant
makes it a ``float``, whereas using the suffix "``L``" makes it a ``long
double``. For example,

.. code-block:: c

   #define A 1.0
   #define B 1.0F
   #define C 1.0L

defines ``A`` to be a ``double`` constant equal to :math:`1.0`, ``B`` to be a
``float`` constant equal to :math:`1.0`, and ``C`` to be a ``long double``
constant equal to :math:`1.0`. The macro call ``RCONST(1.0)`` automatically
expands to ``1.0`` if ``realtype`` is ``double``, to ``1.0F`` if ``realtype`` is
``float``, or to ``1.0L`` if ``realtype`` is ``long double``. SUNDIALS uses the
``RCONST`` macro internally to declare all of its floating-point constants.

Additionally, SUNDIALS defines several macros for common mathematical functions
*e.g.*, ``fabs``, ``sqrt``, ``exp``, etc. in ``sundials_math.h``. The macros are
prefixed with ``SUNR`` and expand to the appropriate ``C`` function based on the
``realtype``. For example, the macro ``SUNRabs`` expands to the ``C`` function
``fabs`` when ``realtype`` is ``double``, ``fabsf`` when ``realtype`` is
``float``, and ``fabsl`` when ``realtype`` is ``long double``.

A user program which uses the type ``realtype``, the ``RCONST`` macro, and the
``SUNR`` mathematical function macros is precision-independent except for any
calls to precision-specific library functions. Our example programs use
``realtype``, ``RCONST``, and the ``SUNR`` macros. Users can, however, use the
type ``double``, ``float``, or ``long double`` in their code (assuming that this
usage is consistent with the typedef for ``realtype``) and call the appropriate
math library functions directly. Thus, a previously existing piece of C or C++
code can use SUNDIALS without modifying the code to use ``realtype``,
``RCONST``, or the ``SUNR`` macros so long as the SUNDIALS libraries are built
to use the corresponding precision (see :numref:`Installation.CMake.Options`).

Integer types used for indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. c:type:: sunindextype

   The type ``sunindextype`` is used for indexing array entries in SUNDIALS
   modules as well as for storing the total problem size (*e.g.*, vector
   lengths and matrix sizes). During configuration ``sunindextype`` may be
   selected to be either a 32- or 64-bit *signed* integer with the default being
   64-bit (see :cmakeop:`SUNDIALS_INDEX_SIZE`).

When using a 32-bit integer the total problem size is limited to
:math:`2^{31}-1` and with 64-bit integers the limit is :math:`2^{63}-1`. For
users with problem sizes that exceed the 64-bit limit an advanced configuration
option is available to specify the type used for ``sunindextype``
(see :cmakeop:`SUNDIALS_INDEX_TYPE`).

A user program which uses ``sunindextype`` to handle indices will work with both
index storage types except for any calls to index storage-specific external
libraries. Our ``C`` and ``C++`` example programs use ``sunindextype``. Users
can, however, use any compatible type (*e.g.*, ``int``, ``long int``,
``int32_t``, ``int64_t``, or ``long long int``) in their code, assuming that
this usage is consistent with the typedef for ``sunindextype`` on their
architecture. Thus, a previously existing piece of C or C++ code can use
SUNDIALS without modifying the code to use ``sunindextype``, so long as the
SUNDIALS libraries use the appropriate index storage type (for details see
:numref:`Installation.CMake.Options`).

Boolean type
~~~~~~~~~~~~

.. c:type:: booleantype

   As ANSI C89 (ISO C90) does not have a built-in boolean data type, SUNDIALS
   defines the type ``booleantype`` as an ``int``.

The advantage of using the name booleantype (instead of int) is an increase in
code readability. It also allows the programmer to make a distinction between
int and boolean data. Variables of type ``booleantype`` are intended to have
only the two values ``SUNFALSE`` (``0``) and ``SUNTRUE`` (``1``).
