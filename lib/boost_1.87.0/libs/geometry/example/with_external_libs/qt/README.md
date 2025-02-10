# ![Boost.Geometry](../../../doc/other/logo/logo_bkg.png)

# Qt

## Introduction

[Qt](https://www.qt.io/product/framework) is a stable and powerful open source framework for developing native cross-platform GUI applications in C++.

## Installing Qt

###  Mac

On a Mac, installing Qt for this purpose is trivial:

`brew install qt`

And then you can use the standard CMake workflow.

During that, either ignore this warning: `Could NOT find WrapVulkanHeaders (missing: Vulkan_INCLUDE_DIR)`

or install vulkan-tools (`brew install vulkan-tools`), though for 2D applications, it might not be necessary.

### Linux Ubuntu

Since Qt 6.3 there is `qt_standard_project_setup()` used in the `CMakeLists.txt`.

If you have an older Qt version (on Ubuntu 22 you might have `QT6.2.4`), replace `qt_standard_project_setup` with:

```
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
```

And then you can use the standard CMake workflow.

If you want to use Qt with Vulkan, or get compiler errors or warnings because of Vulkan missing,
and you want to solve that, then install Vulkan as described [here](https://vulkan.lunarg.com/doc/sdk/1.3.239.0/linux/getting_started_ubuntu.html).

### Other platforms

Either you have Qt installed already, or you install it according to the [Qt documentation](https://doc.qt.io/qt-6/get-and-install-qt.html),
and change the `CMakeLists.txt` if necessary.

## Building this example

Assuming you want to build it with CMake, an example workflow is:

```
cd example/with_external_libs/qt
mkdir my_build_folder
cd my_build_folder
cmake .. -G Ninja
ninja
```

## Running this example

You can pass an Ascii file with WKT polygons as the first command line argument. There are several
packed with Boost.Geometry as examples and as test data.

For example: `././qt_world_mapper.app/Contents/MacOS/qt_world_mapper ../../../data/world.wkt`
