// Copyright 2005 Douglas Gregor.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// A test of the communicator that transmits skeletons and
// content for data types.
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/list.hpp>
#include <boost/mpi/skeleton_and_content.hpp>
#include <boost/mpi/nonblocking.hpp>
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/mpi/collectives/broadcast.hpp>

#include "mpi_test_utils.hpp"

using boost::mpi::communicator;

using boost::mpi::packed_skeleton_iarchive;
using boost::mpi::packed_skeleton_oarchive;

int
test_skeleton_and_content(const communicator& comm, int root,
                          bool manual_broadcast)
{
  using boost::mpi::skeleton;
  using boost::mpi::content;
  using boost::mpi::get_content;
  using boost::make_counting_iterator;
  using boost::mpi::broadcast;

  int failed = 0;
  int list_size = comm.size() + 7;
  if (comm.rank() == root) {
    // Fill in the seed data
    std::list<int> original_list;
    for (int i = 0; i < list_size; ++i)
      original_list.push_back(i);

    std::cout << "Broadcasting integer list skeleton from root " << root
              << "...";
    if (manual_broadcast) {
      // Broadcast the skeleton (manually)
      for (int p = 0; p < comm.size(); ++p)
        if (p != root) comm.send(p, 0, skeleton(original_list));
    } else {
      broadcast(comm, skeleton(original_list), root);
    }
    std::cout << "OK." << std::endl;

    // Broadcast the content (manually)
    std::cout << "Broadcasting integer list content from root " << root
              << "...";
    {
      content c = get_content(original_list);
      for (int p = 0; p < comm.size(); ++p)
        if (p != root) comm.send(p, 1, c);
    }
    std::cout << "OK." << std::endl;

    // Reverse the list, broadcast the content again
    std::reverse(original_list.begin(), original_list.end());
    std::cout << "Broadcasting reversed integer list content from root "
              << root << "...";
    {
      content c = get_content(original_list);
      for (int p = 0; p < comm.size(); ++p)
        if (p != root) comm.send(p, 2, c);
    }
    std::cout << "OK." << std::endl;
  } else {
    // Allocate some useless data, to try to get the addresses of the
    // list<int>'s used later to be different across processes.
    std::list<int> junk_list(comm.rank() * 3 + 1, 17);

    // Receive the skeleton to build up the transferred list
    std::list<int> transferred_list;
    if (manual_broadcast) {
      comm.recv(root, 0, skeleton(transferred_list));
    } else {
      broadcast(comm, skeleton(transferred_list), root);
    }
    BOOST_MPI_CHECK((int)transferred_list.size() == list_size, failed);
    
    // Receive the content and check it
    comm.recv(root, 1, get_content(transferred_list));
    BOOST_MPI_CHECK(std::equal(make_counting_iterator(0),
                           make_counting_iterator(list_size),
                         transferred_list.begin()), failed);
    
    // Receive the reversed content and check it
    comm.recv(root, 2, get_content(transferred_list));
    BOOST_MPI_CHECK(std::equal(make_counting_iterator(0),
                           make_counting_iterator(list_size),
                         transferred_list.rbegin()), failed);
  }

  (comm.barrier)();
  return failed;
}

int
test_skeleton_and_content_nonblocking(const communicator& comm, int root)
{
  using boost::mpi::skeleton;
  using boost::mpi::content;
  using boost::mpi::get_content;
  using boost::make_counting_iterator;
  using boost::mpi::broadcast;
  using boost::mpi::request;
  using boost::mpi::wait_all;
  int failed = 0;
  
  int list_size = comm.size() + 7;
  if (comm.rank() == root) {
    // Fill in the seed data
    std::list<int> original_list;
    for (int i = 0; i < list_size; ++i)
      original_list.push_back(i);

    std::cout << "Non-blocking broadcast of integer list skeleton from root " << root
              << "...";

    // Broadcast the skeleton (manually)
    {
      std::vector<request> reqs;
      for (int p = 0; p < comm.size(); ++p)
        if (p != root) 
          reqs.push_back(comm.isend(p, 0, skeleton(original_list)));
      wait_all(reqs.begin(), reqs.end());
    }
    std::cout << "OK." << std::endl;

    // Broadcast the content (manually)
    std::cout << "Non-blocking broadcast of integer list content from root " << root
              << "...";
    {
      content c = get_content(original_list);
      std::vector<request> reqs;
      for (int p = 0; p < comm.size(); ++p)
        if (p != root) reqs.push_back(comm.isend(p, 1, c));
      wait_all(reqs.begin(), reqs.end());
    }
    std::cout << "OK." << std::endl;

    // Reverse the list, broadcast the content again
    std::reverse(original_list.begin(), original_list.end());
    std::cout << "Non-blocking broadcast of reversed integer list content from root "
              << root << "...";
    {
      std::vector<request> reqs;
      content c = get_content(original_list);
      for (int p = 0; p < comm.size(); ++p)
        if (p != root) reqs.push_back(comm.isend(p, 2, c));
      wait_all(reqs.begin(), reqs.end());
    }
    std::cout << "OK." << std::endl;
  } else {
    // Allocate some useless data, to try to get the addresses of the
    // list<int>'s used later to be different across processes.
    std::list<int> junk_list(comm.rank() * 3 + 1, 17);

    // Receive the skeleton to build up the transferred list
    std::list<int> transferred_list;
    request req = comm.irecv(root, 0, skeleton(transferred_list));
    req.wait();
    BOOST_MPI_CHECK((int)transferred_list.size() == list_size, failed);
    
    // Receive the content and check it
    req = comm.irecv(root, 1, get_content(transferred_list));
    req.wait();
    BOOST_MPI_CHECK(std::equal(make_counting_iterator(0),
                           make_counting_iterator(list_size),
                         transferred_list.begin()), failed);
    
    // Receive the reversed content and check it
    req = comm.irecv(root, 2, get_content(transferred_list));
    req.wait();
    BOOST_MPI_CHECK(std::equal(make_counting_iterator(0),
                           make_counting_iterator(list_size),
                         transferred_list.rbegin()), failed);
  }

  (comm.barrier)();
  return failed;
}

int main()
{
  boost::mpi::environment env;
  communicator comm;
  int failed = 0;
  BOOST_MPI_CHECK(comm.size() > 1, failed);
  if (failed == 0) {
    BOOST_MPI_COUNT_FAILED(test_skeleton_and_content(comm, 0, true), failed);
    BOOST_MPI_COUNT_FAILED(test_skeleton_and_content(comm, 0, false), failed);
    BOOST_MPI_COUNT_FAILED(test_skeleton_and_content(comm, 1, true), failed);
    BOOST_MPI_COUNT_FAILED(test_skeleton_and_content(comm, 1, false), failed);
    BOOST_MPI_COUNT_FAILED(test_skeleton_and_content_nonblocking(comm, 0), failed);
    BOOST_MPI_COUNT_FAILED(test_skeleton_and_content_nonblocking(comm, 1), failed);
  }
  return failed;
}
