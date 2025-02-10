#include <cstdlib>
#include <boost/json/value.hpp>
#include <boost/json/string.hpp>
#include <boost/json/monotonic_resource.hpp>
#include <boost/json/static_resource.hpp>


using namespace boost::json;


int main()
{
    value jv;
    // TEST_EXPR( 'jv', 'null' )

    jv = true;
    // TEST_EXPR( 'jv', 'true' )

    jv = false;
    // TEST_EXPR( 'jv', 'false' )

    jv = 1;
    // TEST_EXPR( 'jv', '1' )

    jv = 1u;
    // TEST_EXPR( 'jv', '1' )

    jv = 1.5;
    // TEST_EXPR( 'jv', '1.5' )

    string js;
    // TEST_EXPR( 'js', '""' )

    js = "1";
    // TEST_EXPR( 'js', '"1"' )

    js = "this is a very long string, unusually long even, definitely not short";
    // TEST_EXPR( 'js', '"this is a very long string, unusually long even, definitely not short"' )

    array ja;
    // TEST_EXPR( 'ja', 'array [size=0, capacity=0]' )

    ja.push_back("a");
    // TEST_EXPR( 'ja', 'array [size=1, capacity=1] = {"a"}' )

    ja.push_back(true);
    // TEST_EXPR( 'ja', 'array [size=2, capacity=2] = {"a", true}' )

    ja.insert(ja.end(), {1, 2, 3, 4});
    // TEST_EXPR( 'ja', 'array [size=6, capacity=6] = {"a", true, 1, 2, 3, 4}' )

    ja.push_back(5);
    // TEST_EXPR( 'ja', 'array [size=7, capacity=9] = {"a", true, 1, 2, 3, 4, 5}' )

    ja[ja.size() - 1] = array{1,2,3};
    // TEST_EXPR( 'ja', 'array [size=7, capacity=9] = {"a", true, 1, 2, 3, 4, array [size=3, capacity=3] = {1, 2, 3}}' )

    object jo;
    // TEST_EXPR( 'jo', 'object [size=0, capacity=0]' )

    jo["a"] = "b";
    // TEST_EXPR( 'jo', 'object [size=1, capacity=1] = {["a"] = "b"}' )

    jo["b"] = "c";
    // TEST_EXPR( 'jo', 'object [size=2, capacity=2] = {["a"] = "b", ["b"] = "c"}' )

    jo.insert({ {"c", "d"}, {"d", "e"} });
    // TEST_EXPR( 'jo', 'object [size=4, capacity=4] = {["a"] = "b", ["b"] = "c", ["c"] = "d", ["d"] = "e"}' )

    jo["e"] = "f";
    // TEST_EXPR( 'jo', 'object [size=5, capacity=6] = {["a"] = "b", ["b"] = "c", ["c"] = "d", ["d"] = "e", ["e"] = "f"}' )

    key_value_pair kv = *jo.begin();
    (void)kv;
    // TEST_EXPR( 'kv', '["a"] = "b"' )

    storage_ptr sp = jv.storage();
    // TEST_EXPR( 'sp', 'storage_ptr [resource=default]' )

    unsigned char buf[1024];
    {
        static_resource sr(buf);
        // TEST_EXPR( 'sr', 'static_resource [buffer={0}, head={0}, free=1024]', '/a &buf' )

        sr.allocate(200);
        unsigned char* new_head = buf + 200;
        (void)new_head;
        // TEST_EXPR( 'sr', 'static_resource [buffer={0}, head={1}, free=824]', '/a &buf', '/a new_head' )

        sp = &sr;
        // TEST_EXPR( 'sp', 'storage_ptr [trivial, resource=static_resource [buffer={0}, head={1}, free=824]]', '/a &buf', '/a new_head' )

        sr.release();
    }

    sp = make_shared_resource<static_resource>(buf);
    // TEST_EXPR( 'sp', 'storage_ptr [trivial, shared, refs=1, resource=static_resource [buffer={0}, head={0}, free=1024]]', '/a &buf' )
    {
        auto sp2 = sp;
        // TEST_EXPR( 'sp', 'storage_ptr [trivial, shared, refs=2, resource=static_resource [buffer={0}, head={0}, free=1024]]', '/a &buf' )
        (void)sp;
    }

    {
        monotonic_resource mr;
        // TEST_EXPR( 'mr', 'monotonic_resource [buffer=0x0, block=0x0, head=0x0, free=0]' )

        storage_ptr sp2 = &mr;
        // TEST_EXPR( 'sp2', 'storage_ptr [trivial, resource=monotonic_resource [buffer=0x0, block=0x0, head=0x0, free=0]]' )
        (void)sp2;
    }

    monotonic_resource mr(buf, 10, sp);
    // TEST_EXPR( 'mr', 'monotonic_resource [buffer={0}, block={0}, head={0}, free=10, upstream=storage_ptr [trivial, shared, refs=2, resource=static_resource [buffer={0}, head={0}, free=1024]]]', '/a &buf' )

    mr.allocate(4);
    unsigned char* new_head = buf + 4;
    (void)new_head;
    // TEST_EXPR( 'mr', 'monotonic_resource [buffer={0}, block={0}, head={1}, free=6, upstream=storage_ptr [trivial, shared, refs=2, resource=static_resource [buffer={0}, head={0}, free=1024]]]', '/a &buf', '/a new_head' )

    return EXIT_SUCCESS;
}
