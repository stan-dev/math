# ![Boost.Geometry](../../../doc/other/logo/logo_bkg.png)

# wxWidgets

## Introduction

[wxWidgets](https://www.wxwidgets.org/) is a stable and powerful open source framework for developing native cross-platform GUI applications in C++.

## Building wxWidgets

There are several possibilities. This documentation uses the CMake approach.

* Retrieve wxWidgets from github
* Be sure to also retrieve the git submodules
* Build and install with cmake, such that it can be found from anywhere.   

```
cd ~git
git clone --recurse-submodules git@github.com:wxWidgets/wxWidgets.git
cd wxWidgets
mkdir my_build_folder
cd my_build_folder
cmake ..
cmake --build .
sudo cmake --build . --target install
```

It is (on macOs) now installed in `/usr/local/lib/`

## Building this example

Assuming you want to build it with CMake

```
cd example/with_external_libs/wxwidgets
mkdir my_build_folder
cd my_build_folder
cmake ..
cmake --build .
```

## Running this example

You can pass an Ascii file with WKT polygons as the first command line argument. There are several
packed with Boost.Geometry as examples and as test data.

For example: `./wx_widgets_world_mapper ../../../data/world.wkt`
