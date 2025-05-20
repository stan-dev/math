+++
title = "Implicit construction"
weight = 10
+++

The preceding code had the compiler stamp out a custom status code domain
for a user supplied `enum`. You now get the following types:

{{% snippet "quick_status_code_from_enum.cpp" "implicit_construction" %}}

As you can see, this is less work than [plugging your custom enum
into `std::error_code`]({{% relref "/motivation/plug_error_code" %}}).
It also has C compatibility, and generates better codegen.
