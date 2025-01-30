// Copyright (C) 2005, 2006 Douglas Gregor.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// A test of the all_reduce() collective.
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <vector>
#include <algorithm>
#include <boost/serialization/string.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <numeric>

using boost::mpi::communicator;

// A simple point class that we can build, add, compare, and
// serialize.
struct point
{
  point() : x(0), y(0), z(0) { }
  point(int x, int y, int z) : x(x), y(y), z(z) { }

  int x;
  int y;
  int z;

 private:
  template<typename Archiver>
  void serialize(Archiver& ar, unsigned int /*version*/)
  {
    ar & x & y & z;
  }

  friend class boost::serialization::access;
};

std::ostream& operator<<(std::ostream& out, const point& p)
{
  return out << p.x << ' ' << p.y << ' ' << p.z;
}

bool operator==(const point& p1, const point& p2)
{
  return p1.x == p2.x && p1.y == p2.y && p1.z == p2.z;
}

bool operator!=(const point& p1, const point& p2)
{
  return !(p1 == p2);
}

point operator+(const point& p1, const point& p2)
{
  return point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}

// test lexical order
bool operator<(const point& p1, const point& p2)
{
  return (p1.x < p2.x 
          ? true
          : (p1.x > p2.x
             ? false 
             : p1.y < p2.y ));
}

namespace boost { namespace mpi {

  template <>
  struct is_mpi_datatype<point> : public mpl::true_ { };

} } // end namespace boost::mpi

template<typename Generator, typename Op>
bool
all_reduce_one_test(const communicator& comm, Generator generator,
                    const char* type_kind, Op op, const char* op_kind,
                    typename Generator::result_type init, bool in_place)
{
  typedef typename Generator::result_type value_type;
  value_type value = generator(comm.rank());

  using boost::mpi::all_reduce;
  using boost::mpi::inplace;

  if (comm.rank() == 0) {
    std::cout << "Reducing to " << op_kind << " of " << type_kind << "...";
    std::cout.flush();
  }

  value_type result_value;
  if (in_place) {
    all_reduce(comm, inplace(value), op);
    result_value = value;
  } else {
    result_value = all_reduce(comm, value, op);
  }
  
  // Compute expected result
  std::vector<value_type> generated_values;
  for (int p = 0; p < comm.size(); ++p)
    generated_values.push_back(generator(p));
  value_type expected_result = std::accumulate(generated_values.begin(),
                                               generated_values.end(),
                                               init, op);
  bool passed = result_value == expected_result;
  if (passed && comm.rank() == 0)
    std::cout << "OK." << std::endl;
  
  comm.barrier();
  return passed;
}

template<typename Generator, typename Op>
bool
all_reduce_array_test(const communicator& comm, Generator generator,
                      const char* type_kind, Op op, const char* op_kind,
                      typename Generator::result_type init, bool in_place)
{
  typedef typename Generator::result_type value_type;
  value_type value = generator(comm.rank());
  std::vector<value_type> send(10, value);

  using boost::mpi::all_reduce;
  using boost::mpi::inplace;

  if (comm.rank() == 0) {
      char const* place = in_place ? "in place" : "out of place";
      std::cout << "Reducing (" << place << ") array to " << op_kind << " of " << type_kind << "...";
      std::cout.flush();
  }
  std::vector<value_type> result;
  if (in_place) {
    all_reduce(comm, inplace(&(send[0])), send.size(), op);
    result.swap(send);
  } else {
    std::vector<value_type> recv(10, value_type());
    all_reduce(comm, &(send[0]), send.size(), &(recv[0]), op);
    result.swap(recv);
  }

  // Compute expected result
  std::vector<value_type> generated_values;
  for (int p = 0; p < comm.size(); ++p)
    generated_values.push_back(generator(p));
  value_type expected_result = std::accumulate(generated_values.begin(),
                                               generated_values.end(),
                                               init, op);
  
  bool passed = (std::equal_range(result.begin(), result.end(), expected_result)
                 == std::make_pair(result.begin(), result.end()));
  if (passed && comm.rank() == 0)
    std::cout << "OK." << std::endl;
  
  comm.barrier();
  return passed;
}

// Test the 4 families of all reduce: (value, array) X (in place, out of place)
// Return numbers of failed tests
template<typename Generator, typename Op>
int
all_reduce_test(const communicator& comm, Generator generator,
                const char* type_kind, Op op, const char* op_kind,
                typename Generator::result_type init)
{
  const bool in_place = true;
  const bool out_of_place = false;
  
  int failed = 0;
  failed += int(!all_reduce_one_test(comm, generator, type_kind, op, op_kind,  init, in_place));
  failed += int(!all_reduce_one_test(comm, generator, type_kind, op, op_kind,  init, out_of_place));
  failed += int(!all_reduce_array_test(comm, generator, type_kind, op, op_kind, init, in_place));
  failed += int(!all_reduce_array_test(comm, generator, type_kind, op, op_kind, init, out_of_place));
  return failed;
}

// Generates integers to test with all_reduce()
struct int_generator
{
  typedef int result_type;

  int_generator(int base = 1) : base(base) { }

  int operator()(int p) const { return base + p; }

 private:
  int base;
};

// Generate points to test with all_reduce()
struct point_generator
{
  typedef point result_type;

  point_generator(point origin) : origin(origin) { }

  point operator()(int p) const
  {
    return point(origin.x + 1, origin.y + 1, origin.z + 1);
  }

 private:
  point origin;
};

struct string_generator
{
  typedef std::string result_type;

  std::string operator()(int p) const
  {
    std::string result = boost::lexical_cast<std::string>(p);
    result += " rosebud";
    if (p != 1) result += 's';
    return result;
  }
};

struct secret_int_bit_and
{
  int operator()(int x, int y) const { return x & y; }
};

struct wrapped_int
{
  wrapped_int() : value(0) { }
  explicit wrapped_int(int value) : value(value) { }

  template<typename Archive>
  void serialize(Archive& ar, unsigned int /* version */)
  {
    ar & value;
  }

  int value;
};

wrapped_int operator+(const wrapped_int& x, const wrapped_int& y)
{
  return wrapped_int(x.value + y.value);
}

bool operator==(const wrapped_int& x, const wrapped_int& y)
{
  return x.value == y.value;
}

bool operator<(const wrapped_int& x, const wrapped_int& y)
{
  return x.value < y.value;
}

// Generates wrapped_its to test with all_reduce()
struct wrapped_int_generator
{
  typedef wrapped_int result_type;

  wrapped_int_generator(int base = 1) : base(base) { }

  wrapped_int operator()(int p) const { return wrapped_int(base + p); }

 private:
  int base;
};

namespace boost { namespace mpi {

// Make std::plus<wrapped_int> commutative.
template<>
struct is_commutative<std::plus<wrapped_int>, wrapped_int>
  : mpl::true_ { };

} } // end namespace boost::mpi

int
main() {
  using namespace boost::mpi;
  environment env;
  communicator comm;
  
  int failed = 0;
  // Built-in MPI datatypes with built-in MPI operations
  failed += all_reduce_test(comm, int_generator(), "integers", std::plus<int>(), "sum", 0);
  failed += all_reduce_test(comm, int_generator(), "integers", std::multiplies<int>(), "product", 1);
  failed += all_reduce_test(comm, int_generator(), "integers", maximum<int>(), "maximum", 0);
  failed += all_reduce_test(comm, int_generator(), "integers", minimum<int>(), "minimum", 2);
  
  // User-defined MPI datatypes with operations that have the
  // same name as built-in operations.
  failed += all_reduce_test(comm, point_generator(point(0,0,0)), "points", std::plus<point>(),
                            "sum", point());

  // Built-in MPI datatypes with user-defined operations
  failed += all_reduce_test(comm, int_generator(17), "integers", secret_int_bit_and(), 
                            "bitwise and", -1);
  
  // Arbitrary types with user-defined, commutative operations.
  failed += all_reduce_test(comm, wrapped_int_generator(17), "wrapped integers",
                            std::plus<wrapped_int>(), "sum", wrapped_int(0));

  // Arbitrary types with (non-commutative) user-defined operations
  failed += all_reduce_test(comm, string_generator(), "strings",
                            std::plus<std::string>(), "concatenation", std::string());
  if (comm.rank() == 0) {
    if (failed > 0) {
      std::cout << "Failed " << failed  << " tests.\n";
    } else {
      std::cout << "PASSED.\n";
    }
  }
  return failed;
}
