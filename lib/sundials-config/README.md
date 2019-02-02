# Sundials for Stan Notes, Feb 2nd 2019

The official sundials source code contains printf, sprintf and fprintf
calls which we have to make optional when compiling the library. This
is needed due to restrictions imposed by CRAN which disallows these
function calls to even appear in any binary part of a R package. To
achieve this we modify the vanilla sundials source files such adequate
replacements are done. The modified sundials source has all of these
calls swapped by a macro call which is defined in
sundials-config/include/stan_sundials_printf_override.hpp.

To update to a new sundials version the respective replacements are
handled by the prepare.sh script. To use this script, please do

- Download the sundials tar.gz file and place it in stan-math/lib
- Edit the prepare.sh script such that
   - The version is updated in the SUNDIALS_VERSION variable
   - This should be all what needs to change in this script, but if sundials conventions change, then you need to adapt that
- Change into the stan-math/lib directory and call the ./sundials-config/prepare.sh script
- This will extract the new source code and apply all changes
- In case of a large version change for sundials, then the you need to recreate a new
  sundials-config/include/sundials/sundials_config.h file
- To do so please follow the build instructions for sundials to setup their cmake system with the default settings. This step is in most cases not needed.
- Finally remove with git the old sundials and check in the new sundials version
- Update the sundials version in make/libraries and README.md

Done.
