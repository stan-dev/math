for f in "*.hpp"
do
    if grep -q "/meta/" $f
    then
       echo $f
       sed -i "/\/meta\//d" $f
       sed -i "4i#include <stan/math/prim/meta.hpp>" $f
    fi
done
