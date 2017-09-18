#ifndef VIENNACL_DEVICE_SPECIFIC_EXECUTION_HANDLER_HPP
#define VIENNACL_DEVICE_SPECIFIC_EXECUTION_HANDLER_HPP

/* =========================================================================
   Copyright (c) 2010-2016, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


/** @file viennacl/device_specific/execution_handler.hpp
    @brief Helper for handling fallbacks, lazy compilation, input-dependent kernels, etc
*/

#include <map>

#include "viennacl/tools/shared_ptr.hpp"

#include "viennacl/device_specific/lazy_program_compiler.hpp"
#include "viennacl/device_specific/templates/template_base.hpp"
#include "viennacl/device_specific/utils.hpp"

namespace viennacl
{
namespace device_specific
{

class execution_handler
{
public:
  typedef std::map< std::string, tools::shared_ptr<template_base> > container_type;

private:
  std::string append_prefix(std::string const & str)
  {
    return "_" + str;
  }

  std::string define_extension(std::string const & ext)
  {
    // Note: On devices without double precision support, 'ext' is an empty string.
    return (ext.length() > 1) ? std::string("#pragma OPENCL EXTENSION " + ext + " : enable\n") : std::string("\n");
  }

  void init_program_compiler(std::string const & name, bool force_recompilation)
  {
    lazy_programs_.push_back(lazy_program_compiler(&ctx_, name, force_recompilation));
    lazy_programs_.back().add(define_extension(device_.double_support_extension()));
  }

public:
  execution_handler(std::string const & program_name_base, viennacl::ocl::context & ctx, viennacl::ocl::device const & device, bool force_recompilation = false) : ctx_(ctx), device_(device), program_names_(2)
  {
    lazy_programs_.reserve(2);
    init_program_compiler(program_name_base + "_0", force_recompilation);
    init_program_compiler(program_name_base + "_1", force_recompilation);
  }

  void add(std::string const & key, template_base const & T, statements_container const & statements)
  {
    if (kernels_.insert(container_type::value_type(key, T.clone())).second)
    {
      std::vector<std::string> sources = at(kernels_, key)->generate(append_prefix(key), statements, device_);
      assert(sources.size()<=2);
      for (unsigned int i = 0; i < sources.size(); ++i)
        lazy_programs_[i].add(sources[i]);
    }
  }

  template_base * template_of(std::string const & key)
  {
    return at(kernels_, key).get();
  }

  void execute(container_type::key_type const & key, statements_container const & statements)
  {
    tools::shared_ptr<template_base> & template_pointer = at(kernels_, key);
    template_pointer->enqueue(append_prefix(key), lazy_programs_, statements);
  }

private:
  viennacl::ocl::context & ctx_;
  viennacl::ocl::device const & device_;
  container_type kernels_;
  std::vector<std::string> program_names_;
  std::vector<lazy_program_compiler> lazy_programs_;
};

}
}
#endif
