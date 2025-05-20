#ifndef BOOST_UNORDERED_HEADER
#error "this test requires a class template be passed as a macro"
#endif

#define STRINGIZE(text) STRINGIZE_I(text)
#define STRINGIZE_I(...) #__VA_ARGS__

#define HEADER STRINGIZE(BOOST_UNORDERED_HEADER)
#include HEADER
