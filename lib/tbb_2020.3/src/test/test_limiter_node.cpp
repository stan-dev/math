/*
    Copyright (c) 2005-2020 Intel Corporation

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/

#include "harness.h"

#if __TBB_CPF_BUILD
#define TBB_DEPRECATED_FLOW_NODE_EXTRACTION 1
#include "harness_graph.h"
#endif

#include "tbb/flow_graph.h"
#include "tbb/atomic.h"
#include "tbb/task_scheduler_init.h"
#include "test_follows_and_precedes_api.h"

const int L = 10;
const int N = 1000;

using tbb::flow::internal::SUCCESSFULLY_ENQUEUED;

template< typename T >
struct serial_receiver : public tbb::flow::receiver<T>, NoAssign {
   T next_value;
   tbb::flow::graph& my_graph;

   serial_receiver(tbb::flow::graph& g) : next_value(T(0)), my_graph(g) {}

   tbb::task *try_put_task( const T &v ) __TBB_override {
       ASSERT( next_value++  == v, NULL );
       return const_cast<tbb::task *>(SUCCESSFULLY_ENQUEUED);
   }

    tbb::flow::graph& graph_reference() const __TBB_override {
        return my_graph;
    }

#if TBB_DEPRECATED_FLOW_NODE_EXTRACTION
    typedef typename tbb::flow::receiver<T>::built_predecessors_type built_predecessors_type;
    typedef typename tbb::flow::receiver<T>::predecessor_list_type predecessor_list_type;
    typedef typename tbb::flow::receiver<T>::predecessor_type predecessor_type;
    built_predecessors_type bpt;
    built_predecessors_type &built_predecessors() __TBB_override { return bpt; }
    void internal_add_built_predecessor( predecessor_type & ) __TBB_override { }
    void internal_delete_built_predecessor( predecessor_type & ) __TBB_override { }
    void copy_predecessors( predecessor_list_type & ) __TBB_override { }
    size_t predecessor_count() __TBB_override { return 0; }
#endif

   void reset_receiver(tbb::flow::reset_flags /*f*/) __TBB_override {next_value = T(0);}
};

template< typename T >
struct parallel_receiver : public tbb::flow::receiver<T>, NoAssign {

    tbb::atomic<int> my_count;
    tbb::flow::graph& my_graph;

    parallel_receiver(tbb::flow::graph& g) : my_graph(g) { my_count = 0; }

    tbb::task *try_put_task( const T &/*v*/ ) __TBB_override {
       ++my_count;
       return const_cast<tbb::task *>(tbb::flow::internal::SUCCESSFULLY_ENQUEUED);
    }

    tbb::flow::graph& graph_reference() const __TBB_override {
        return my_graph;
    }

#if TBB_DEPRECATED_FLOW_NODE_EXTRACTION
    typedef typename tbb::flow::receiver<T>::built_predecessors_type built_predecessors_type;
    typedef typename tbb::flow::receiver<T>::predecessor_list_type predecessor_list_type;
    typedef typename tbb::flow::receiver<T>::predecessor_type predecessor_type;
    built_predecessors_type bpt;
    built_predecessors_type &built_predecessors() __TBB_override { return bpt; }
    void internal_add_built_predecessor( predecessor_type & ) __TBB_override { }
    void internal_delete_built_predecessor( predecessor_type & ) __TBB_override { }
    void copy_predecessors( predecessor_list_type & ) __TBB_override { }
    size_t predecessor_count( ) __TBB_override { return 0; }
#endif
    void reset_receiver(tbb::flow::reset_flags /*f*/) __TBB_override {my_count = 0;}
};

template< typename T >
struct empty_sender : public tbb::flow::sender<T> {
        typedef typename tbb::flow::sender<T>::successor_type successor_type;

        bool register_successor( successor_type & ) __TBB_override { return false; }
        bool remove_successor( successor_type & ) __TBB_override { return false; }
#if TBB_DEPRECATED_FLOW_NODE_EXTRACTION
        typedef typename tbb::flow::sender<T>::built_successors_type built_successors_type;
        typedef typename tbb::flow::sender<T>::successor_list_type successor_list_type;
        built_successors_type bst;
        built_successors_type &built_successors() __TBB_override { return bst; }
        void    internal_add_built_successor( successor_type & ) __TBB_override { }
        void internal_delete_built_successor( successor_type & ) __TBB_override { }
        void copy_successors( successor_list_type & ) __TBB_override { }
        size_t successor_count() __TBB_override { return 0; }
#endif
};


template< typename T >
struct put_body : NoAssign {

    tbb::flow::limiter_node<T> &my_lim;
    tbb::atomic<int> &my_accept_count;

    put_body( tbb::flow::limiter_node<T> &lim, tbb::atomic<int> &accept_count ) :
        my_lim(lim), my_accept_count(accept_count) {}

    void operator()( int ) const {
        for ( int i = 0; i < L; ++i ) {
            bool msg = my_lim.try_put( T(i) );
            if ( msg == true )
               ++my_accept_count;
        }
    }
};

template< typename T >
struct put_dec_body : NoAssign {

    tbb::flow::limiter_node<T> &my_lim;
    tbb::atomic<int> &my_accept_count;

    put_dec_body( tbb::flow::limiter_node<T> &lim, tbb::atomic<int> &accept_count ) :
        my_lim(lim), my_accept_count(accept_count) {}

    void operator()( int ) const {
        int local_accept_count = 0;
        while ( local_accept_count < N ) {
            bool msg = my_lim.try_put( T(local_accept_count) );
            if ( msg == true ) {
                ++local_accept_count;
                ++my_accept_count;
                my_lim.decrement.try_put( tbb::flow::continue_msg() );
            }
        }
    }

};

template< typename T >
void test_puts_with_decrements( int num_threads, tbb::flow::limiter_node< T >& lim , tbb::flow::graph& g) {
    parallel_receiver<T> r(g);
    empty_sender< tbb::flow::continue_msg > s;
    tbb::atomic<int> accept_count;
    accept_count = 0;
    tbb::flow::make_edge( lim, r );
    tbb::flow::make_edge(s, lim.decrement);
#if TBB_DEPRECATED_FLOW_NODE_EXTRACTION
    ASSERT(lim.decrement.predecessor_count() == 1, NULL);
    ASSERT(lim.successor_count() == 1, NULL);
    ASSERT(lim.predecessor_count() == 0, NULL);
    typename tbb::flow::interface11::internal::decrementer
        <tbb::flow::limiter_node<T>, tbb::flow::continue_msg>::predecessor_list_type dec_preds;
    lim.decrement.copy_predecessors(dec_preds);
    ASSERT(dec_preds.size() == 1, NULL);
#endif
    // test puts with decrements
    NativeParallelFor( num_threads, put_dec_body<T>(lim, accept_count) );
    int c = accept_count;
    ASSERT( c == N*num_threads, NULL );
    ASSERT( r.my_count == N*num_threads, NULL );
}

//
// Tests
//
// limiter only forwards below the limit, multiple parallel senders / single receiver
// multiple parallel senders that put to decrement at each accept, limiter accepts new messages
//
//
template< typename T >
int test_parallel(int num_threads) {

   // test puts with no decrements
   for ( int i = 0; i < L; ++i ) {
       tbb::flow::graph g;
       tbb::flow::limiter_node< T > lim(g, i);
       parallel_receiver<T> r(g);
       tbb::atomic<int> accept_count;
       accept_count = 0;
       tbb::flow::make_edge( lim, r );
       // test puts with no decrements
       NativeParallelFor( num_threads, put_body<T>(lim, accept_count) );
       g.wait_for_all();
       int c = accept_count;
       ASSERT( c == i, NULL );
   }

   // test puts with decrements
   for ( int i = 1; i < L; ++i ) {
       tbb::flow::graph g;
       tbb::flow::limiter_node< T > lim(g, i);
       test_puts_with_decrements(num_threads, lim, g);
       tbb::flow::limiter_node< T > lim_copy( lim );
       test_puts_with_decrements(num_threads, lim_copy, g);
   }

   return 0;
}

//
// Tests
//
// limiter only forwards below the limit, single sender / single receiver
// at reject, a put to decrement, will cause next message to be accepted
//
template< typename T >
int test_serial() {

   // test puts with no decrements
   for ( int i = 0; i < L; ++i ) {
       tbb::flow::graph g;
       tbb::flow::limiter_node< T > lim(g, i);
       serial_receiver<T> r(g);
       tbb::flow::make_edge( lim, r );
       for ( int j = 0; j < L; ++j ) {
           bool msg = lim.try_put( T(j) );
           ASSERT( ( j < i && msg == true ) || ( j >= i && msg == false ), NULL );
       }
       g.wait_for_all();
   }

   // test puts with decrements
   for ( int i = 1; i < L; ++i ) {
       tbb::flow::graph g;
       tbb::flow::limiter_node< T > lim(g, i);
       serial_receiver<T> r(g);
       empty_sender< tbb::flow::continue_msg > s;
       tbb::flow::make_edge( lim, r );
       tbb::flow::make_edge(s, lim.decrement);
       for ( int j = 0; j < N; ++j ) {
           bool msg = lim.try_put( T(j) );
           ASSERT( ( j < i && msg == true ) || ( j >= i && msg == false ), NULL );
           if ( msg == false ) {
               lim.decrement.try_put( tbb::flow::continue_msg() );
               msg = lim.try_put( T(j) );
               ASSERT( msg == true, NULL );
           }
       }
   }
   return 0;
}

// reported bug in limiter (http://software.intel.com/en-us/comment/1752355)
#define DECREMENT_OUTPUT 1  // the port number of the decrement output of the multifunction_node
#define LIMITER_OUTPUT 0    // port number of the integer output

typedef tbb::flow::multifunction_node<int, tbb::flow::tuple<int,tbb::flow::continue_msg> > mfnode_type;

tbb::atomic<size_t> emit_count;
tbb::atomic<size_t> emit_sum;
tbb::atomic<size_t> receive_count;
tbb::atomic<size_t> receive_sum;

struct mfnode_body {
    int max_cnt;
    tbb::atomic<int>* my_cnt;
    mfnode_body(const int& _max, tbb::atomic<int> &_my) : max_cnt(_max), my_cnt(&_my)  { }
    void operator()(const int &/*in*/, mfnode_type::output_ports_type &out) {
        int lcnt = ++(*my_cnt);
        if(lcnt > max_cnt) {
            return;
        }
        // put one continue_msg to the decrement of the limiter.
        if(!tbb::flow::get<DECREMENT_OUTPUT>(out).try_put(tbb::flow::continue_msg())) {
            ASSERT(false,"Unexpected rejection of decrement");
        }
        {
            // put messages to the input of the limiter_node until it rejects.
            while( tbb::flow::get<LIMITER_OUTPUT>(out).try_put(lcnt) ) {
                emit_sum += lcnt;
                ++emit_count;
            }
        }
    }
};

struct fn_body {
    int operator()(const int &in) {
        receive_sum += in;
        ++receive_count;
        return in;
    }
};

//                   +------------+
//    +---------+    |            v
//    | mf_node |0---+       +----------+          +----------+
// +->|         |1---------->| lim_node |--------->| fn_node  |--+
// |  +---------+            +----------+          +----------+  |
// |                                                             |
// |                                                             |
// +-------------------------------------------------------------+
//
void
test_multifunction_to_limiter(int _max, int _nparallel) {
    tbb::flow::graph g;
    emit_count = 0;
    emit_sum = 0;
    receive_count = 0;
    receive_sum = 0;
    tbb::atomic<int> local_cnt;
    local_cnt = 0;
    mfnode_type mf_node(g, tbb::flow::unlimited, mfnode_body(_max, local_cnt));
    tbb::flow::function_node<int, int> fn_node(g, tbb::flow::unlimited, fn_body());
    tbb::flow::limiter_node<int> lim_node(g, _nparallel);
    tbb::flow::make_edge(tbb::flow::output_port<LIMITER_OUTPUT>(mf_node), lim_node);
    tbb::flow::make_edge(tbb::flow::output_port<DECREMENT_OUTPUT>(mf_node), lim_node.decrement);
    tbb::flow::make_edge(lim_node, fn_node);
    tbb::flow::make_edge(fn_node, mf_node);
#if TBB_DEPRECATED_FLOW_NODE_EXTRACTION
    REMARK("pred cnt == %d\n",(int)(lim_node.predecessor_count()));
    REMARK("succ cnt == %d\n",(int)(lim_node.successor_count()));
    tbb::flow::limiter_node<int>::successor_list_type my_succs;
    lim_node.copy_successors(my_succs);
    REMARK("succ cnt from vector  == %d\n",(int)(my_succs.size()));
    tbb::flow::limiter_node<int>::predecessor_list_type my_preds;
    lim_node.copy_predecessors(my_preds);
    REMARK("pred cnt from vector  == %d\n",(int)(my_preds.size()));
#endif
    mf_node.try_put(1);
    g.wait_for_all();
    ASSERT(emit_count == receive_count, "counts do not match");
    ASSERT(emit_sum == receive_sum, "sums do not match");

    // reset, test again
    g.reset();
    emit_count = 0;
    emit_sum = 0;
    receive_count = 0;
    receive_sum = 0;
    local_cnt = 0;;
    mf_node.try_put(1);
    g.wait_for_all();
    ASSERT(emit_count == receive_count, "counts do not match");
    ASSERT(emit_sum == receive_sum, "sums do not match");
}


void
test_continue_msg_reception() {
    tbb::flow::graph g;
    tbb::flow::limiter_node<int> ln(g,2);
    tbb::flow::queue_node<int>   qn(g);
    tbb::flow::make_edge(ln, qn);
    ln.decrement.try_put(tbb::flow::continue_msg());
    ln.try_put(42);
    g.wait_for_all();
    int outint;
    ASSERT(qn.try_get(outint) && outint == 42, "initial put to decrement stops node");
}

#if TBB_DEPRECATED_LIMITER_NODE_CONSTRUCTOR
using namespace tbb::flow;
void run_and_check_result(graph& g, limiter_node<int>& limit, queue_node<int>& queue, broadcast_node<continue_msg>& broad) {
    ASSERT( limit.try_put(1), NULL );
    ASSERT( limit.try_put(2), NULL );
    ASSERT( !limit.try_put(3), NULL );
    ASSERT( broad.try_put(continue_msg()), NULL );
    ASSERT( limit.decrement.try_put(continue_msg()), NULL );
    ASSERT( limit.try_put(4), NULL );
    ASSERT( !limit.try_put(5), NULL );
    g.wait_for_all();

    int list[] = {1, 2, 4};
    int var = 0;
    for (size_t i = 0; i < sizeof(list)/sizeof(list[0]); i++) {
        queue.try_get(var);
        ASSERT(var==list[i], "some data dropped, input does not match output");
    }
}

void test_num_decrement_predecessors() {
    graph g;
    queue_node<int> output_queue(g);
    limiter_node<int> limit1(g, 2, /*number_of_predecessors*/1);
    limiter_node<int, continue_msg> limit2(g, 2, /*number_of_predecessors*/1);
    broadcast_node<continue_msg> broadcast(g);

    make_edge(limit1, output_queue);
    make_edge(limit2, output_queue);

    make_edge(broadcast, limit1.decrement);
    make_edge(broadcast, limit2.decrement);

    run_and_check_result(g, limit1, output_queue, broadcast);
    run_and_check_result(g, limit2, output_queue, broadcast);
}
#else // TBB_DEPRECATED_LIMITER_NODE_CONSTRUCTOR
//
// This test ascertains that if a message is not successfully put
// to a successor, the message is not dropped but released.
//

void test_reserve_release_messages() {
    using namespace tbb::flow;
    graph g;

    //making two queue_nodes: one broadcast_node and one limiter_node
    queue_node<int> input_queue(g);
    queue_node<int> output_queue(g);
    broadcast_node<int> broad(g);
    limiter_node<int, int> limit(g,2); //threshold of 2

    //edges
    make_edge(input_queue, limit);
    make_edge(limit, output_queue);
    make_edge(broad,limit.decrement);

    int list[4] = {19, 33, 72, 98}; //list to be put to the input queue

    input_queue.try_put(list[0]); // succeeds
    input_queue.try_put(list[1]); // succeeds
    input_queue.try_put(list[2]); // fails, stored in upstream buffer
    g.wait_for_all();

    remove_edge(limit, output_queue); //remove successor

    //sending message to the decrement port of the limiter
    broad.try_put(1); //failed message retrieved.
    g.wait_for_all();

    make_edge(limit, output_queue); //putting the successor back

    broad.try_put(1);  //drop the count

    input_queue.try_put(list[3]);  //success
    g.wait_for_all();

    int var=0;

    for (int i=0; i<4; i++) {
        output_queue.try_get(var);
        ASSERT(var==list[i], "some data dropped, input does not match output");
        g.wait_for_all();
    }
}

void test_decrementer() {
    const int threshold = 5;
    tbb::flow::graph g;
    tbb::flow::limiter_node<int, int> limit(g, threshold);
    tbb::flow::queue_node<int> queue(g);
    make_edge(limit, queue);
    int m = 0;
    ASSERT( limit.try_put( m++ ), "Newly constructed limiter node does not accept message." );
    ASSERT( limit.decrement.try_put( -threshold ), // close limiter's gate
            "Limiter node decrementer's port does not accept message." );
    ASSERT( !limit.try_put( m++ ), "Closed limiter node's accepts message." );
    ASSERT( limit.decrement.try_put( threshold + 5 ),  // open limiter's gate
            "Limiter node decrementer's port does not accept message." );
    for( int i = 0; i < threshold; ++i )
        ASSERT( limit.try_put( m++ ), "Limiter node does not accept message while open." );
    ASSERT( !limit.try_put( m ), "Limiter node's gate is not closed." );
    g.wait_for_all();
    int expected[] = {0, 2, 3, 4, 5, 6};
    int actual = -1; m = 0;
    while( queue.try_get(actual) )
        ASSERT( actual == expected[m++], NULL );
    ASSERT( sizeof(expected) / sizeof(expected[0]) == m, "Not all messages have been processed." );
    g.wait_for_all();

    const size_t threshold2 = size_t(-1);
    tbb::flow::limiter_node<int, long long> limit2(g, threshold2);
    make_edge(limit2, queue);
    ASSERT( limit2.try_put( 1 ), "Newly constructed limiter node does not accept message." );
    long long decrement_value = (long long)( size_t(-1)/2 );
    ASSERT( limit2.decrement.try_put( -decrement_value ),
            "Limiter node decrementer's port does not accept message" );
    ASSERT( limit2.try_put( 2 ), "Limiter's gate should not be closed yet." );
    ASSERT( limit2.decrement.try_put( -decrement_value ),
            "Limiter node decrementer's port does not accept message" );
    ASSERT( !limit2.try_put( 3 ), "Overflow happened for internal counter." );
    int expected2[] = {1, 2};
    actual = -1; m = 0;
    while( queue.try_get(actual) )
        ASSERT( actual == expected2[m++], NULL );
    ASSERT( sizeof(expected2) / sizeof(expected2[0]) == m, "Not all messages have been processed." );
    g.wait_for_all();
}
#endif // TBB_DEPRECATED_LIMITER_NODE_CONSTRUCTOR

#if TBB_DEPRECATED_FLOW_NODE_EXTRACTION
void test_extract() {
    tbb::flow::graph g;
    int j;
    tbb::flow::limiter_node<int> node0(g, /*threshold*/1);
    tbb::flow::queue_node<int> q0(g);
    tbb::flow::queue_node<int> q1(g);
    tbb::flow::queue_node<int> q2(g);
    tbb::flow::broadcast_node<tbb::flow::continue_msg> b0(g);
    tbb::flow::broadcast_node<tbb::flow::continue_msg> b1(g);

    for( int i = 0; i < 2; ++i ) {
        REMARK("At pass %d\n", i);
        ASSERT(node0.predecessor_count() == 0, "incorrect predecessor count at start");
        ASSERT(node0.successor_count() == 0, "incorrect successor count at start");
        ASSERT(node0.decrement.predecessor_count() == 0, "incorrect decrement pred count at start");

        tbb::flow::make_edge(q0, node0);
        tbb::flow::make_edge(q1, node0);
        tbb::flow::make_edge(node0, q2);
        tbb::flow::make_edge(b0, node0.decrement);
        tbb::flow::make_edge(b1, node0.decrement);
        g.wait_for_all();

        /*    b0   b1              */
        /*      \  |               */
        /*  q0\  \ |               */
        /*     \  \|               */
        /*      +-node0---q2       */
        /*     /                   */
        /*  q1/                    */

        q0.try_put(i);
        g.wait_for_all();
        ASSERT(node0.predecessor_count() == 2, "incorrect predecessor count after construction");
        ASSERT(node0.successor_count() == 1, "incorrect successor count after construction");
        ASSERT(node0.decrement.predecessor_count() == 2, "incorrect decrement pred count after construction");
        ASSERT(q2.try_get(j), "fetch of value forwarded to output queue failed");
        ASSERT(j == i, "improper value forwarded to output queue");
        q0.try_put(2*i);
        g.wait_for_all();
        ASSERT(!q2.try_get(j), "limiter_node forwarded item improperly");
        b0.try_put(tbb::flow::continue_msg());
        g.wait_for_all();
        ASSERT(!q2.try_get(j), "limiter_node forwarded item improperly");
        b0.try_put(tbb::flow::continue_msg());
        g.wait_for_all();
        ASSERT(q2.try_get(j) && j == 2*i, "limiter_node failed to forward item");

        tbb::flow::limiter_node<int>::successor_list_type sv;
        tbb::flow::limiter_node<int>::predecessor_list_type pv;
        tbb::flow::continue_receiver::predecessor_list_type dv;
        tbb::flow::limiter_node<int>::successor_list_type sv1;
        tbb::flow::limiter_node<int>::predecessor_list_type pv1;
        tbb::flow::continue_receiver::predecessor_list_type dv1;

        node0.copy_predecessors(pv);
        node0.copy_successors(sv);
        node0.decrement.copy_predecessors(dv);
        pv1.push_back(&(q0));
        pv1.push_back(&(q1));
        sv1.push_back(&(q2));
        dv1.push_back(&(b0));
        dv1.push_back(&(b1));

        ASSERT(pv.size() == 2, "improper size for predecessors");
        ASSERT(sv.size() == 1, "improper size for successors");
        ASSERT(lists_match(pv,pv1), "predecessor lists do not match");
        ASSERT(lists_match(sv,sv1), "successor lists do not match");
        ASSERT(lists_match(dv,dv1), "successor lists do not match");

        if(i == 0) {
            node0.extract();
            ASSERT(node0.predecessor_count() == 0, "incorrect predecessor count after extraction");
            ASSERT(node0.successor_count() == 0, "incorrect successor count after extraction");
            ASSERT(node0.decrement.predecessor_count() == 0, "incorrect decrement pred count after extraction");
        }
        else {
            q0.extract();
            b0.extract();
            q2.extract();

            ASSERT(node0.predecessor_count() == 1, "incorrect predecessor count after extract second iter");
            ASSERT(node0.successor_count() == 0, "incorrect successor count after extract second iter");
            ASSERT(node0.decrement.predecessor_count() == 1, "incorrect decrement pred count after extract second iter");

            node0.copy_predecessors(pv);
            node0.copy_successors(sv);
            node0.decrement.copy_predecessors(dv);
            pv1.clear();
            sv1.clear();
            dv1.clear();
            pv1.push_back(&(q1));
            dv1.push_back(&(b1));

            ASSERT(lists_match(pv,pv1), "predecessor lists do not match second iter");
            ASSERT(lists_match(sv,sv1), "successor lists do not match second iter");
            ASSERT(lists_match(dv,dv1), "successor lists do not match second iter");

            q1.extract();
            b1.extract();
        }
        ASSERT(node0.predecessor_count() == 0, "incorrect predecessor count after extract");
        ASSERT(node0.successor_count() == 0, "incorrect successor count after extract");
        ASSERT(node0.decrement.predecessor_count() == 0, "incorrect decrement pred count after extract");

    }

}
#endif  // TBB_DEPRECATED_FLOW_NODE_EXTRACTION

#if __TBB_PREVIEW_FLOW_GRAPH_NODE_SET
#include <array>
#include <vector>
void test_follows_and_precedes_api() {
    using msg_t = tbb::flow::continue_msg;

    std::array<msg_t, 3> messages_for_follows= { {msg_t(), msg_t(), msg_t()} };
    std::vector<msg_t> messages_for_precedes = {msg_t()};

    follows_and_precedes_testing::test_follows
        <msg_t, tbb::flow::limiter_node<msg_t, msg_t>>(messages_for_follows, 1000);
    follows_and_precedes_testing::test_precedes
        <msg_t, tbb::flow::limiter_node<msg_t, msg_t>>(messages_for_precedes, 1000);

}
#endif

#if __TBB_CPP17_DEDUCTION_GUIDES_PRESENT
void test_deduction_guides() {
    using namespace tbb::flow;

    graph g;
    broadcast_node<int> br(g);
    limiter_node<int> l0(g, 100);

#if __TBB_PREVIEW_FLOW_GRAPH_NODE_SET
    limiter_node l1(follows(br), 100);
    static_assert(std::is_same_v<decltype(l1), limiter_node<int>>);

    limiter_node l2(precedes(br), 100);
    static_assert(std::is_same_v<decltype(l2), limiter_node<int>>);
#endif

    limiter_node l3(l0);
    static_assert(std::is_same_v<decltype(l3), limiter_node<int>>);
}
#endif

int TestMain() {
    for (int i = 1; i <= 8; ++i) {
        tbb::task_scheduler_init init(i);
        test_serial<int>();
        test_parallel<int>(i);
    }
    test_continue_msg_reception();
    test_multifunction_to_limiter(30,3);
    test_multifunction_to_limiter(300,13);
    test_multifunction_to_limiter(3000,1);
#if TBB_DEPRECATED_LIMITER_NODE_CONSTRUCTOR
    test_num_decrement_predecessors();
#else
    test_reserve_release_messages();
    test_decrementer();
#endif
#if TBB_DEPRECATED_FLOW_NODE_EXTRACTION
    test_extract();
#endif
#if __TBB_PREVIEW_FLOW_GRAPH_NODE_SET
    test_follows_and_precedes_api();
#endif
#if __TBB_CPP17_DEDUCTION_GUIDES_PRESENT
    test_deduction_guides();
#endif
   return Harness::Done;
}
