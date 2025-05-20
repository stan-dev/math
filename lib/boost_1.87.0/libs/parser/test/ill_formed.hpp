// Copyright (C) 2020 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
// Copyright (C) 2019 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_PARSER_ILL_FORMED_HPP
#define BOOST_PARSER_ILL_FORMED_HPP

#include <boost/parser/detail/detection.hpp>


template<template<class...> class Template, typename... Args>
using ill_formed = std::integral_constant<
    bool,
    !boost::parser::detail::is_detected<Template, Args...>::value>;

template<template<class...> class Template, typename... Args>
using well_formed = std::integral_constant<
    bool,
    boost::parser::detail::is_detected<Template, Args...>::value>;

#endif
