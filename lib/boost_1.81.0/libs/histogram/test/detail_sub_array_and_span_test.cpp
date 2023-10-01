#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/detail/nonmember_container_access.hpp>
#include <boost/histogram/detail/reduce_command.hpp>
#include <boost/histogram/detail/span.hpp>
#include <boost/histogram/detail/sub_array.hpp>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include "throw_exception.hpp"

namespace boost {
namespace histogram {
namespace detail {
template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const sub_array<T, N>&) {
  return os;
}
std::ostream& operator<<(std::ostream& os, const reduce_command&) { return os; }
} // namespace detail
} // namespace histogram
} // namespace boost

using namespace boost::histogram::detail;

int main() {

  {
    sub_array<int, 3> a = {1, 2};

    BOOST_TEST_EQ(a.size(), 2);
    BOOST_TEST_EQ(a.max_size(), 3);
    BOOST_TEST_EQ(a.at(0), 1);
    BOOST_TEST_EQ(a.at(1), 2);

    sub_array<int, 3> b = {1, 2};
    BOOST_TEST_EQ(a, b);

    sub_array<int, 3> c = {1, 2, 3};
    BOOST_TEST_NE(a, c);

    sub_array<int, 3> d = {2, 2};
    BOOST_TEST_NE(a, d);

    auto sp = span<int>(a);
    BOOST_TEST_EQ(sp.size(), 2);
    BOOST_TEST_EQ(sp.front(), 1);
    BOOST_TEST_EQ(sp.back(), 2);

    const auto& ca = a;
    auto csp = span<const int>(ca);
    BOOST_TEST_EQ(csp.size(), 2);
    BOOST_TEST_EQ(csp.front(), 1);
    BOOST_TEST_EQ(csp.back(), 2);
  }

  {
    sub_array<int, 2> a(2, 1);
    sub_array<int, 2> b(1, 2);
    std::swap(a, b);

    BOOST_TEST_EQ(a.size(), 1);
    BOOST_TEST_EQ(b.size(), 2);
    BOOST_TEST_EQ(a[0], 2);
    BOOST_TEST_EQ(b[0], 1);
    BOOST_TEST_EQ(b[1], 1);
  }

  {
    const std::initializer_list<int> a = {1, 2};
    auto sp = span<const int>(a);
    BOOST_TEST_EQ(sp.size(), 2);
  }

  {
    const sub_array<reduce_command, 3> a(2);
    auto sp = span<const reduce_command>(a);
    BOOST_TEST_EQ(sp.size(), 2);
  }

  {
    const std::initializer_list<reduce_command> a = {reduce_command{}};
    auto sp = span<const reduce_command>(a);
    BOOST_TEST_EQ(sp.size(), 1);
  }

  return boost::report_errors();
}