
// Copyright 2008-2009 Daniel James.
// Copyright 2024 Braden Ganetsky.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or move at http://www.boost.org/LICENSE_1_0.txt)

#if !defined(BOOST_UNORDERED_TEST_HELPERS_COUNT_HEAD)
#define BOOST_UNORDERED_TEST_HELPERS_COUNT_HEAD

#include <boost/core/lightweight_test.hpp>
#include <boost/container_hash/hash.hpp>

namespace test {
  struct object_count
  {
    int instances;
    int constructions;

    object_count() : instances(0), constructions(0) {}
    void reset() { *this = object_count(); }

    void construct()
    {
      ++instances;
      ++constructions;
    }

    void destruct()
    {
      if (instances == 0) {
        BOOST_ERROR("Unbalanced constructions.");
      } else {
        --instances;
      }
    }

    bool operator==(object_count const& x) const
    {
      return instances == x.instances && constructions == x.constructions;
    }

    bool operator!=(object_count const& x) const { return !(*this == x); }

    friend std::ostream& operator<<(std::ostream& out, object_count const& c)
    {
      out << "[instances: " << c.instances
          << ", constructions: " << c.constructions << "]";
      return out;
    }
  };

  // This won't be a problem as I'm only using a single compile unit
  // in each test (this is actually require by the minimal test
  // framework).
  //
  // boostinspect:nounnamed
  namespace {
    object_count global_object_count;
  }

  struct counted_object
  {
    counted_object() { global_object_count.construct(); }
    counted_object(counted_object const&) { global_object_count.construct(); }
    ~counted_object() { global_object_count.destruct(); }
  };

  struct check_instances
  {
    int instances_;
    int constructions_;

    check_instances()
        : instances_(global_object_count.instances),
          constructions_(global_object_count.constructions)
    {
    }
    ~check_instances()
    {
      BOOST_TEST(global_object_count.instances == instances_);
    }

    int instances() const { return global_object_count.instances - instances_; }
    int constructions() const
    {
      return global_object_count.constructions - constructions_;
    }
  };

  struct smf_count
  {
    int default_constructions = 0;
    int copy_constructions = 0;
    int move_constructions = 0;
    int copy_assignments = 0;
    int move_assignments = 0;
    int destructions = 0;

#if (BOOST_CXX_VERSION < 201402L) || (defined(_MSC_VER) && _MSC_VER < 1910)
    smf_count() = default;

    smf_count(int default_constructions_, int copy_constructions_,
      int move_constructions_, int copy_assignments_, int move_assignments_,
      int destructions_)
        : default_constructions(default_constructions_),
          copy_constructions(copy_constructions_),
          move_constructions(move_constructions_),
          copy_assignments(copy_assignments_),
          move_assignments(move_assignments_), destructions(destructions_)
    {
    }
#endif

    void reset() { *this = smf_count(); }

    void default_construct() { ++default_constructions; }
    void copy_construct() { ++copy_constructions; }
    void move_construct() { ++move_constructions; }
    void copy_assign() { ++copy_assignments; }
    void move_assign() { ++move_assignments; }
    void destruct() { ++destructions; }

    friend bool operator==(smf_count const& lhs, smf_count const& rhs)
    {
      return lhs.default_constructions == rhs.default_constructions &&
             lhs.copy_constructions == rhs.copy_constructions &&
             lhs.move_constructions == rhs.move_constructions &&
             lhs.copy_assignments == rhs.copy_assignments &&
             lhs.move_assignments == rhs.move_assignments &&
             lhs.destructions == rhs.destructions;
    }

    friend std::ostream& operator<<(std::ostream& out, smf_count const& c)
    {
      out << "[default_constructions: " << c.default_constructions
          << ", copy_constructions: " << c.copy_constructions
          << ", move_constructions: " << c.move_constructions
          << ", copy_assignments: " << c.copy_assignments
          << ", move_assignments: " << c.move_assignments
          << ", destructions: " << c.destructions << "]";
      return out;
    }
  };

  template <class Tag> class smf_counted_object
  {
  public:
    static smf_count count;
    static void reset_count() { count.reset(); }

    smf_counted_object(int index) : smf_counted_object() { index_ = index; }

    smf_counted_object() : index_(++running_index)
    {
      count.default_construct();
    }
    smf_counted_object(smf_counted_object const& rhs) : index_(rhs.index_)
    {
      count.copy_construct();
    }
    smf_counted_object(smf_counted_object&& rhs) noexcept : index_(rhs.index_)
    {
      count.move_construct();
    }
    smf_counted_object& operator=(smf_counted_object const& rhs)
    {
      count.copy_assign();
      index_ = rhs.index_;
      return *this;
    }
    smf_counted_object& operator=(smf_counted_object&& rhs) noexcept
    {
      count.move_assign();
      index_ = rhs.index_;
      return *this;
    }
    ~smf_counted_object() { count.destruct(); }

    friend bool operator==(
      smf_counted_object const& lhs, smf_counted_object const& rhs)
    {
      return lhs.index_ == rhs.index_;
    }

    friend std::size_t hash_value(smf_counted_object const& x)
    {
      return boost::hash<int>()(x.index_);
    }

    int index_;

  private:
    static int running_index;
  };
  template <class Tag> smf_count smf_counted_object<Tag>::count = {};
  template <class Tag> int smf_counted_object<Tag>::running_index = 0;
} // namespace test

#endif
