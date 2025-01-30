//
// Copyright (c) 2022 alandefreitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//

#include <boost/url.hpp>

int main() {
    boost::urls::ipv4_address ip_addr({127, 0, 0, 1});
    if (ip_addr.to_bytes().size() != 4) {
        return 1;
    }
    return 0;
}
