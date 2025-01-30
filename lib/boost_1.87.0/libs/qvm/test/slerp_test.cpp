// Copyright 2008-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef NDEBUG
#undef NDEBUG
#endif

#ifdef BOOST_QVM_TEST_SINGLE_HEADER
#   include BOOST_QVM_TEST_SINGLE_HEADER
#else
#   include <boost/qvm/quat_operations.hpp>
#   include <boost/qvm/quat.hpp>
#endif

#include "test_qvm_quaternion.hpp"
#include "gold.hpp"

namespace
    {
    void
    test()
        {
        using namespace boost::qvm;
            {
            test_qvm::quaternion<Q1> a=rotx_quat(1.0f);
            test_qvm::quaternion<Q1> b=rotx_quat(.5f);
            test_qvm::quaternion<Q1> aa=slerp360(a,b,0);
            test_qvm::quaternion<Q1> bb=slerp360(a,b,1);
            BOOST_QVM_TEST_CLOSE(aa.a,a.a,0.0001f);
            BOOST_QVM_TEST_CLOSE(bb.a,b.a,0.0001f);
            }
        for( float a1=0; a1<6.28f; a1+=0.1f )
            {
            test_qvm::quaternion<Q1> const qx1=rotx_quat(a1);
            test_qvm::quaternion<Q1> const qy1=roty_quat(a1);
            test_qvm::quaternion<Q1> const qz1=rotz_quat(a1);
            for( float a2=0; a2<6.28f; a2+=0.1f )
                {
                test_qvm::quaternion<Q1> const qx2=rotx_quat(a2);
                test_qvm::quaternion<Q1> const qy2=roty_quat(a2);
                test_qvm::quaternion<Q1> const qz2=rotz_quat(a2);
                for( float t=0; t<1; t+=0.1f )
                    {
                        {
                        test_qvm::quaternion<Q1> const qx=rotx_quat(a1*(1-t)+a2*t);
                        test_qvm::quaternion<Q1> const qy=roty_quat(a1*(1-t)+a2*t);
                        test_qvm::quaternion<Q1> const qz=rotz_quat(a1*(1-t)+a2*t);
                        test_qvm::quaternion<Q1> const qsx=slerp360(qx1,qx2,t);
                        test_qvm::quaternion<Q1> const qsy=slerp360(qref(qy1),qy2,t);
                        test_qvm::quaternion<Q1> const qsz=slerp360(qz1,qref(qz2),t);
                        BOOST_QVM_TEST_CLOSE(qx.a,qsx.a,0.001f);
                        BOOST_QVM_TEST_CLOSE(qy.a,qsy.a,0.001f);
                        BOOST_QVM_TEST_CLOSE(qz.a,qsz.a,0.001f);
                        }
                        {
                        test_qvm::quaternion<Q1> const x1=slerp180(qx1,qx2,t);
                        test_qvm::quaternion<Q1> const x2=slerp360(dot(qx1,qx2)<0 ? -qx1 : qx1,qx2,t);
                        test_qvm::quaternion<Q1> const y1=slerp180(qy1,qy2,t);
                        test_qvm::quaternion<Q1> const y2=slerp360(dot(qy1,qy2)<0 ? -qy1 : qy1,qy2,t);
                        test_qvm::quaternion<Q1> const z1=slerp180(qz1,qz2,t);
                        test_qvm::quaternion<Q1> const z2=slerp360(dot(qz1,qz2)<0 ? -qz1 : qz1,qz2,t);
                        BOOST_QVM_TEST_CLOSE(x1.a,x2.a,0.001f);
                        BOOST_QVM_TEST_CLOSE(y1.a,y2.a,0.001f);
                        BOOST_QVM_TEST_CLOSE(z1.a,z2.a,0.001f);
                        }
                    }
                }
            }
        }
    }

int
main()
    {
    test();
    return boost::report_errors();
    }
