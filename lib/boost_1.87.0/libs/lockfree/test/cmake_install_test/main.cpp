#include <boost/lockfree/queue.hpp>

int main()
{
    boost::lockfree::queue< int > q( 64 );
}
