// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/leaf/handle_errors.hpp>

namespace leaf = boost::leaf;

leaf::diagnostic_details f();
leaf::diagnostic_details a = f();
leaf::diagnostic_details b = a;
