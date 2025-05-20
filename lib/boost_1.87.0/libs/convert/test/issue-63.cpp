#include <boost/convert.hpp>
#include <boost/convert/charconv.hpp>
#include <iostream>

int main()
{
    auto cnv4 = boost::cnv::charconv();
    std::string PuNum = "2";
    auto optPuNum = boost::convert<size_t>(PuNum, cnv4);

    if (optPuNum)
        std::cout << *optPuNum << "\n";
    else
        std::cout << "nullopt\n";

    return 0;
}
