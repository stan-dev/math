//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_LOCALE_TEST_POSIX_TOOLS_HPP
#define BOOST_LOCALE_TEST_POSIX_TOOLS_HPP

#include <clocale>
#include <string>
#if defined(__APPLE__) || defined(__FreeBSD__)
#    include <xlocale.h>
#endif

#ifndef LC_ALL_MASK
using locale_t = void*;
locale_t newlocale(int, const char*, locale_t)
{
    return nullptr;
}
void freelocale(locale_t) {}
#    define LC_ALL_MASK 0xFFFFFFFF
#endif

class locale_holder {
    locale_t l_;
    void reset(const locale_t l = nullptr)
    {
        if(l_)
            freelocale(l_);
        l_ = l;
    }

public:
    explicit locale_holder(locale_t l = nullptr) : l_(l) {}
    ~locale_holder() { reset(); }
    locale_holder(const locale_holder&) = delete;
    locale_holder& operator=(const locale_holder&) = delete;
    locale_holder& operator=(locale_t l)
    {
        reset(l);
        return *this;
    }
    operator locale_t() const { return l_; }
};

inline bool has_posix_locale(const std::string& name)
{
    locale_holder l(newlocale(LC_ALL_MASK, name.c_str(), nullptr));
    return !!l;
}

#endif
