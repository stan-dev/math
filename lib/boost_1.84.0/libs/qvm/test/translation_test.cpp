// Copyright 2008-2022 Emil Dotchevski and Reverge Studios, Inc.

// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_QVM_TEST_SINGLE_HEADER
#   include BOOST_QVM_TEST_SINGLE_HEADER
#else
#   include <boost/qvm/vec_operations.hpp>
#   include <boost/qvm/vec.hpp>
#   include <boost/qvm/map_mat_vec.hpp>
#endif

#include "test_qvm_matrix.hpp"
#include "test_qvm_vector.hpp"
#include "gold.hpp"

namespace
    {
    template <int R, int C>
    void
    test()
        {
        using namespace boost::qvm;
        test_qvm::matrix<M1,R,C> x(42,1);
        test_qvm::vector<V1,C-1> y=translation(x);
        for( int i=0; i!=C-1; ++i )
            y.b[i]=x.a[i][C-1];
        BOOST_QVM_TEST_EQ(y.a,y.b);
        translation(x) *= 2;
        for( int i=0; i!=C-1; ++i )
            x.b[i][C-1] *= 2;
        BOOST_QVM_TEST_EQ(x.a,x.b);
        translation(x) + translation(x);
        -translation(x);
        }
    }

int
main()
    {
    test<3,3>();
    test<2,3>();
    test<4,4>();
    test<3,4>();
    test<5,5>();
    test<4,5>();
    return boost::report_errors();
    }
