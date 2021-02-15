# Release Notes

All notable changes to [Boost.GIL](https://github.com/boostorg/gil/) project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

## [1.72.0] - 2019-12-11

### Added
- GSoC 2019: Lanczos resampling for image down scaling ([PR #309](https://github.com/boostorg/gil/pull/309)).
- GSoC 2019: Methods for binary thresholding, inverted binary thresholding and truncation thresholding ([PR #313](https://github.com/boostorg/gil/pull/313)).
- GSoC 2019: Otsu thresholding method ([PR #314](https://github.com/boostorg/gil/pull/314)).
- GSoC 2019: Adaptive thresholding using mean or gaussian-weighted sum of the neighbourhood area ([PR #315](https://github.com/boostorg/gil/pull/315)).
- GSoC 2019: Harris response calculation (corner detector without non-maximum filtering) ([PR #350](https://github.com/boostorg/gil/pull/350)).
- GSoC 2019: Hessian corner detector ([PR #364](https://github.com/boostorg/gil/pull/364)).
- GSoC 2019: Types for defining 2D kernel, `kernel_2d` and `kernel_2d_fixed`, in Numeric extension ([PR #361](https://github.com/boostorg/gil/pull/361)).
- GSoC 2019: Implementation of 2D convolution as new function `convolve_2d` ([PR #367](https://github.com/boostorg/gil/pull/367)).
- GSoC 2019: Box filtering using the average filter ([PR #383](https://github.com/boostorg/gil/pull/383)).
- GSoC 2019: Blur function based on normalized mean filter ([PR #383](https://github.com/boostorg/gil/pull/383)).
- GSoC 2019: Sobel and Scharr operators ([PR #392](https://github.com/boostorg/gil/pull/392)).
- GSoC 2019: Median filter to remove noise from image ([PR #393](https://github.com/boostorg/gil/pull/393)).
- Continued adding new test cases and significantly improved overall test coverage.
- Documented purpose of `cached_location_t` ([PR #287](https://github.com/boostorg/gil/pull/287)).
- Function `convolve_1d` in Numeric extension for convenient use of `convolve_rows` and `convolve_cols` ([PR #347](https://github.com/boostorg/gil/pull/347) and [PR #367](https://github.com/boostorg/gil/pull/367)).
- Function `extend_boundary` in Numeric extension to perform image boundary extension ([PR #386](https://github.com/boostorg/gil/pull/386)).
- Project release notes maintained in Markdown file `RELEASES.md` ([PR #404](https://github.com/boostorg/gil/pull/404)).

### Changed
- Move all tests, core features and extensions, inside `test/` directory ([PR #302](https://github.com/boostorg/gil/pull/302)).

### Removed
- Replace Boost.MPL with Boost.MP11 ([PR #274](https://github.com/boostorg/gil/pull/274)).
- Removed use of Boost.TypeTraits ([PR #274](https://github.com/boostorg/gil/pull/274)).
- Dropped support for GCC <= 4.8 ([PR #296](https://github.com/boostorg/gil/pull/296)).
- Remove `include/boost/gil/version.hpp` file as unused ([PR #403](https://github.com/boostorg/gil/pull/403)).

### Fixed
- Undetermined value of default-initialized channel and pixel objects ([PR #273](https://github.com/boostorg/gil/pull/273)).
- Undefined behaviour due to `std::is_trivially_default_constructible` specializations ([PR #284](https://github.com/boostorg/gil/pull/284)).
- Crash when reading PNG files with an invalid header ([PR #385](https://github.com/boostorg/gil/pull/385)).
- Applied the [Rule of Three](https://en.wikipedia.org/wiki/Rule_of_three_(C%2B%2B_programming)) for numerous types.
- Removed uses of deprecated implicit definition of defaulted copy assignment operator or copy constructor.

## [1.68.0] - 2018-08-09

### Added
- The library now requires a C++11-compliant compiler.
- Added Toolbox extension following the [review and acceptance into Boost](https://lists.boost.org/boost-announce/2011/01/0281.php).

### Changed
- The I/O extensions have been entirely rewritten as I/O v2, [reviewed and accepted into Boost](https://lists.boost.org/boost-announce/2011/01/0281.php).
- Documentation has been reformatted and updated.

### Removed
- The existing I/O v1 extension has been replaced with I/O v2.

## [1.53.0] - 2013-02-04

### Fixed
- Fixed self-assignment warnings (Trac [#4919](https://svn.boost.org/trac10/ticket/4919)).

## [1.35.0] - 2008-03-29

### Added
- First Boost release of Generic Image Library developed by Lubomir Bourdev and Hailin Jin following the [review and acceptance into Boost](https://lists.boost.org/Archives/boost/2006/11/112896.php).

---------------------------------------------------------------------

## Pre-Boost History of the Generic Image Library (GIL) by Adobe

The log of changes prior the first release of GIL as part of Boost
was collected from https://stlab.adobe.com/gil/news.html site and
linked PDF documents with detailed changes.

---------------------------------------------------------------------

## [2.1.1] - 2007-09-15

### Changed
- Swapped template arguments for `color_element_type`, `color_element_reference_type` and `color_element_const_reference_type` to take `ColorBase` first for consistency with the other similar metafunctions.

### Fixed
- Minor bugs fixed.

## [2.1.0] - 2007-06-17

### Added
- Added support for accessing raw memory from image views by getting raw pointer to the beginning of the memory
  associated with a homogeneous image view using new functions `interleaved_view_get_raw_data` or `planar_view_get_raw_data`.
- Support for non-byte-aligned images (e.g. 6-bit RGB222, or 1-bit grayscale).
  To support bit distance, we are using the same classes that were providing byte distance (`byte_addressible_step_iterator`, `byte_addressible_2d_locator`, etc.)
  except that now they operate on memory units instead of bytes.
  A memory unit can currently be either a byte or a bit.
- New `byte_to_memunit` function required by the `MemoryBasedIteratorConcept`, which specifies the number of bits per memory unit (either 1 or 8).
- New classes for references and iterators over bit-aligned pixels: `bit_aligned_pixel_reference`, `bit_aligned_pixel_iterator`.
  The memory unit of bit aligned pixel iterators is a bit, i.e. `byte_to_unit<bit_aligned_pixel_iterator<T> >::value == 8`.
- The `value_type` of a bit-aligned image is a `packed_pixel` (new name, see below).
  A packed pixel is a pixel that is byte-aligned but whose channels may not be byte aligned.
  There is a strong analogy with the way interleaved and planar images are implemented, with `packed_pixel` corresponding
  to `pixel`, `bit_aligned_pixel_reference` corresponding to `planar_pixel_reference`
  and `bit_aligned_pixel_iterator` corresponding to `planar_pixel_iterator`.
- New metafunction `bit_aligned_image_type` for constructing bit-aligned images.
  A bit-aligned image is an image whose pixels may not be byte-aligned (such as an RGB222 image).
- New metafunction `pixel_value_type` for constructing homogenous pixel value from elements.
- New metafunction `packed_pixel_type` for constructing homogenous packed pixel from elements.
- New metafunction `packed_image_type` for constructing packed images with packed pixel as its `value_type`.

### Changed
- Renamed `heterogeneous_packed_pixel` to `packed_pixel`.
- Renamed `ByteAdvancableIteratorConcept` to `MemoryBasedIteratorConcept`.
- Renamed `byte_addressable_{step_iterator,2d_locator}` to `memory_based_{step_iterator,2d_locator}`.
- Renamed `byte_{advance,advanced,distance,step}` to `memunit_{advance,advanced,distance,step}`,
- Renamed `locator::row_bytes()` to `locator::row_size()` and `locator::pix_bytestep()` to `locator::pixel_size()`.
- Simplified `packed_channel_reference` and `packed_dynamic_channel_reference` by removing the `BitField` parameter (it is now computed automatically).
- Improved `channel_convert` - it is faster by switching to floating-point math only if necessary.

### Fixed
- Fixed a roundoff bug in the conversion (related to floating-point math switching).
- Fixed histogram regression tests.

## [2.0.x] - 2007-03-27

### Changed
- Minor bug fixes.
- Regression test improvements.

### Removed
- Removed any external dependencies from the regression tests.

## [2.0.0] - 2007-03-08

### Added
- Further Boost integration:
  - Directories follow the Boost convention.
  - Different models are usually now split in separate files.
  - Renamed some files for better consistency.
  - Renamed classes, functions and template arguments with longer but clearer and more consistent names.
- New `deprecated.hpp` - a file that maps many of the deprecated names to current ones.
  Including it will help porting your code to GIL 2.0. After porting to GIL 2.0, however,
  make sure that your code works when this file is not included.
- New `swap` function required for reference proxies, since the `std::swap` default does not do the right thing.
- Metafunctions `iterator_type_from_pixel` and `view_type_from_pixel` to allow creating standard iterators and views associated with a pixel type.
- New `scoped_channel_value`, a channel adaptor that changes the operational range of a channel. `bits32f` is defined as a `float` with range `0.0` to `1.0`.
- New `packed_channel_value`, `packed_channel_reference` and `packed_dynamic_channel_reference` which model channels operating on bit ranges.
- New `heterogeneous_packed_pixel`, a model of a pixel whose channels are bit ranges (e.g. 16-bit RGB pixel in the 565 format).
- Metafunctions to get the k-th element of a color base (or its reference): `kth_semantic_element_type`, `kth_semantic_element_reference_type`, `kth_semantic_element_const_reference_type`.
- Metafunctions to operate on pixel iterator: `const_iterator_type`, `iterator_is_mutable`, `is_iterator_adaptor`.
- New image view algorithms `uninitialized_fill_pixels`, `uninitialized_copy_pixels` and method `is_1d_traversable`.
- Added support for creating images with a new value to fill.

### Changed
- Updated the design guide and tutorial, updated syntax of concepts to the latest concepts proposal.
- In `image`, `image_view`, `any_image`, `any_image_view`:  
  There are no longer global functions `get_width()`, `get_height()`, `get_dimensions()`, `num_channels()`.
  Use  methods `width()`, `height()`, `dimensions()` instead.
- In models of pixel, pixel iterator, pixel locator, image view and image:  
  There used to be different ways of getting to a pixel, channel, color space, etc. of an image view,
  pixel, locator, iterator and image (e.g. traits, member typedefs).
  Now all pixel-based GIL constructs (pixels, pixel iterators, locators, image views and images) model
  `PixelBasedConcept`, which means they provide the following metafunctions: `color_space_type`, `channel_mapping_type`, `is_planar`, `num_channels`
  and for homogeneous constructs we also have: `channel_type`.
  To get the pixel type or pixel reference/const reference type of an image, image view, locator
  and pixel, use member typedefs `value_type`, `reference` and `const_reference`.
- In `locator`, `image`, `image_view`, `any_image` and `any_image_view`:  
  Removed `dynamic_x_step_t`, `dynamic_y_step_t`, `dynamic_xy_step_t` and `dynamic_xy_step_transposed_t`
  as member typedefs of locators and image views.
  Instead, there are separate concepts `HasDynamicXStepTypeConcept`, `HasDynamicYStepTypeConcept`,
  `HasTransposedTypeConcept` which all GIL-provided locators, views and images model.
  Those concepts require a metafunction to get the corresponding type.
  Analogously, all GIL pixel iterators model `HasDynamicXStepTypeConcept`.
- In channel, the min and max value is now part of the channel traits.
  For all built-in types the channel range equals the physical range (as determined by `std::numeric_traits<T>::max()`).
- Provide `channel_convert` support to convert between any of the GIL-provided channel types.
  The operation is also consistent - conversion is done as a linear mapping that maps the min/max to the min/max.
- In pixel, major redesign of pixel-level constructs.
  Renamed `color_base` to `homogeneous_color_base` and defined it once, not for each color space.
  The color base is a first-class concept and allows to model any bundle of color elements.
  Work needed to define a new color space has been simplified a lot.
  All former pixel-level algorithms and accessors now operate on color bases.
  The elements of a color base can be accessed by physical or semantic index or by name (channel names
  can no longer be accessed as members of the pixel e.g. `my_pixel.gray = 0`), use `get_color` instead).
- In color base, algorithms now can take heterogeneous pixels (i.e. pixels each channel of which may have a different type).
  The `color_convert` can operate on heterogeneous pixels with the exception of to/from RGBA.
- In image, the class `image` is no longer templated over the image view. It is now templated over pixel value.
- In dynamic image, instead of removed `cross_vector_image_types` and `cross_vector_image_view_types`, create MPL vector to enumerate types.
- Renamed algorithms `{copy,equal,fill,for_each,generate,max,min,transform}_channels` to `static_{copy,equal,fill,for_each,generate,max,min,transform}`.
- Rename metafunctions `channel` to `at_c`, `semantic_channel` to `semantic_at_c`, `get_nth_channel` to `dynamic_at_c`.
- Renamed `planar_{ptr,ref}` to `planar_pixel_{iterator,reference}`.
- Renamed `PixelConcept` to `HomogeneousPixelConcept`.
- Renamed `HeterogeneousPixelConcept` to `PixelConcept`.
- Renamed `pixel_image_iterator` to `iterator_from_2d`.
- Renamed `is_contiguous` to `is_1d_traversable`.
- Renamed `membased_2d_locator` to `byte_addressable_2d_locator`.
- Renamed `resize_clobber_image` to `image::recreate`.

### Fixed
- Now compiles with GCC 4.1.1.
- Fixed some bugs in defining reference proxies.

### Removed
- Flattened the `core` directory as part of Boost integration.
- Got rid of channel accessors from pixel types.
- Got rid of `pixel_traits`. Use nested typedefs `value_type`, `reference` and `const_reference` or metafunctions implementing `PixelBasedConcept`.
- Got rid of `pixel_iterator_traits`. Use `std::iterator_traits`, `PixelBasedConcept` metafunctions or the new metafunctions for pixel iterators.
- Got rid of the ability to directly access pixels through image, only through views. The image no longer models STL's random access container concept.
- No more LAB and HSB color space, because there is no color conversion support implemented for these. New color spaces can be added with just a few lines of code.

## [1.x] - 2007-01-03

### Added
- Restored back the ability to assign a channel to a grayscale pixel.

### Changed
- Fixed some minor issues with color converted views of dynamic images.

## [1.x] - 2006-11-07

GIL accepted to Boost.

GIL's Boost review was successful and GIL will be part of the Boost libraries.
It will most likely first appear in the 1.35 version of Boost.
In the future our web page will continue to provide you with the latest
improvements to GIL, as we have the flexibility to release more frequently than Boost.

## [1.x] - 2006-10-20

### Added
- New regression tests.
- Source code of usage examples is available to download from the website.

### Changed
- Minor changes to GIL core.

## [1.x] - 2006-10-02

### Added
- First version of the Numeric extension.
  The extension provides some basic image processing algorithms, such as convolution and resampling.
- Introduction of pixel traits.

### Changed
- Improved consistent use of MPL predicates and standardized template parameter names.

## [1.x] - 2006-09-20

### Added
- GIL now allows users to overload the default color conversion with one of their own.
- Section in the design guide describes how to overload default color conversion.

### Changed
- Color conversion improvements.

## [1.0] - 2006-08-29

### Added
- Pixel dereference adaptors are introduced.
- Locator concepts/models are made more generic.
- Example of creating the Mandelbrot set is described in the tutorial.

### Changed
- It is now easier to construct virtual image views.

## [] - 2006-06-27

### Added
- A GIL Flash presentation is posted (aka video lecture).

## [] - 2006-06-14

### Added
- GIL homepage goes live
