#ifndef VIENNACL_DEVICE_SPECIFIC_FORWARDS_H
#define VIENNACL_DEVICE_SPECIFIC_FORWARDS_H

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


/** @file viennacl/device_specific/forwards.h
    @brief Forwards declaration
*/

#include <list>
#include <map>
#include <set>
#include <stdexcept>

#include "viennacl/scheduler/io.hpp"

#include "viennacl/ocl/forwards.h"
#include "viennacl/tools/shared_ptr.hpp"
#include "viennacl/scheduler/forwards.h"

#include "viennacl/backend/mem_handle.hpp"

namespace viennacl
{
namespace device_specific
{

//Error codes
static const int TEMPLATE_VALID = 0;
static const int TEMPLATE_LOCAL_MEMORY_OVERFLOW = -1;
static const int TEMPLATE_WORK_GROUP_SIZE_OVERFLOW = -2;
static const int TEMPLATE_LOCAL_SIZE_0_OVERFLOW = -3;
static const int TEMPLATE_LOCAL_SIZE_1_OVERFLOW = -4;
static const int TEMPLATE_LOCAL_SIZE_2_OVERFLOW = -5;
static const int TEMPLATE_LOCAL_SIZE_NOT_WARP_MULTIPLE = -6;
static const int TEMPLATE_INVALID_SIMD_WIDTH = -7;
static const int TEMPLATE_INVALID_FETCHING_POLICY_TYPE= -9;

static const int TEMPLATE_GLOBAL_MEMORY_REQUIRES_ZERO_LOCAL_FETCH = -10;
static const int TEMPLATE_MS_NS_MUST_BE_SIMD_WIDTH_MULTIPLE = -11;
static const int TEMPLATE_KS_MUST_BE_SMALLER_THAN_KL = -12;
static const int TEMPLATE_SIMD_WIDTH_MUST_BE_ONE = -13;
static const int TEMPLATE_LOCAL_FETCH_PRODUCT_MUST_MATCH_LOCAL_SIZE_PRODUCT = -14;
static const int TEMPLATE_LOCAL_FETCH_0_MUST_BE_KL_MULTIPLE = -15;
static const int TEMPLATE_LOCAL_FETCH_0_MUST_BE_NL_MULTIPLE = -16;
static const int TEMPLATE_LOCAL_FETCH_1_MUST_BE_KL_MULTIPLE = -17;
static const int TEMPLATE_LOCAL_FETCH_1_MUST_BE_ML_MULTIPLE = -18;

struct index_tuple
{
  index_tuple(std::string const & _i, std::string const & _bound0) : i(_i), bound0(_bound0), j(""), bound1(""){ }
  index_tuple(std::string const & _i, std::string const & _bound0, std::string const & _j, std::string const & _bound1) : i(_i), bound0(_bound0), j(_j), bound1(_bound1){ }
  std::string i;
  std::string bound0;
  std::string j;
  std::string bound1;
};

inline bool is_scalar_reduction(scheduler::statement_node const & node)
{
  return node.op.type==scheduler::OPERATION_BINARY_INNER_PROD_TYPE || node.op.type_family==scheduler::OPERATION_VECTOR_REDUCTION_TYPE_FAMILY;
}

inline bool is_vector_reduction(scheduler::statement_node const & node)
{
  return node.op.type==scheduler::OPERATION_BINARY_MAT_VEC_PROD_TYPE
      || node.op.type_family==scheduler::OPERATION_ROWS_REDUCTION_TYPE_FAMILY
      || node.op.type_family==scheduler::OPERATION_COLUMNS_REDUCTION_TYPE_FAMILY;
}

inline scheduler::statement_node const & lhs_most(scheduler::statement::container_type const & array, vcl_size_t root)
{
  scheduler::statement_node const * current = &array[root];
  while (current->lhs.type_family==scheduler::COMPOSITE_OPERATION_FAMILY)
    current = &array[current->lhs.node_index];
  return *current;
}

enum expression_type
{
  SCALAR_AXPY_TYPE,
  VECTOR_AXPY_TYPE,
  MATRIX_AXPY_TYPE,
  REDUCTION_TYPE,
  ROW_WISE_REDUCTION_Nx_TYPE,
  ROW_WISE_REDUCTION_Tx_TYPE,
  MATRIX_PRODUCT_NN_TYPE,
  MATRIX_PRODUCT_TN_TYPE,
  MATRIX_PRODUCT_NT_TYPE,
  MATRIX_PRODUCT_TT_TYPE,
  INVALID_EXPRESSION_TYPE
};

inline const char * expression_type_to_string(expression_type type)
{
  switch (type)
  {
  case SCALAR_AXPY_TYPE : return "Scalar AXPY";
  case VECTOR_AXPY_TYPE : return "Vector AXPY";
  case MATRIX_AXPY_TYPE : return "Matrix AXPY";
  case REDUCTION_TYPE : return "Reduction";
  case ROW_WISE_REDUCTION_Nx_TYPE : return "Row-wise reduction: Ax";
  case ROW_WISE_REDUCTION_Tx_TYPE : return "Row-wise reduction : Tx";
  case MATRIX_PRODUCT_NN_TYPE : return "Matrix-Matrix Product : AA";
  case MATRIX_PRODUCT_TN_TYPE : return "Matrix-Matrix Product : TA";
  case MATRIX_PRODUCT_NT_TYPE : return "Matrix-Matrix Product : AT";
  case MATRIX_PRODUCT_TT_TYPE : return "Matrix-Matrix Product : TT";
  default : return "INVALID EXPRESSION";
  }
}

/** @brief generate the string for a pointer kernel argument */
static std::string generate_value_kernel_argument(std::string const & scalartype, std::string const & name)
{
  return scalartype + ' ' + name + ",";
}

/** @brief generate the string for a pointer kernel argument */
static std::string generate_pointer_kernel_argument(std::string const & address_space, std::string const & scalartype, std::string const & name)
{
  return address_space +  " " + scalartype + "* " + name + ",";
}

/** @brief Emulation of C++11's .at() member for std::map<>, const-version */
template<typename KeyT, typename ValueT>
ValueT const & at(std::map<KeyT, ValueT> const & map, KeyT const & key)
{
  typename std::map<KeyT, ValueT>::const_iterator it = map.find(key);
  if (it != map.end())
    return it->second;

  throw std::out_of_range("Generator: Key not found in map");
}

/** @brief Emulation of C++11's .at() member for std::map<>, non-const version */
template<typename KeyT, typename ValueT>
ValueT & at(std::map<KeyT, ValueT> & map, KeyT const & key)
{
  typename std::map<KeyT, ValueT>::iterator it = map.find(key);
  if (it != map.end())
    return it->second;

  throw std::out_of_range("Generator: Key not found in map");
}

/** @brief Exception for the case the generator is unable to deal with the operation */
class generator_not_supported_exception : public std::exception
{
public:
  generator_not_supported_exception() : message_() {}
  generator_not_supported_exception(std::string message) : message_("ViennaCL: Internal error: The generator cannot handle the statement provided: " + message) {}
  virtual const char* what() const throw() { return message_.c_str(); }
  virtual ~generator_not_supported_exception() throw() {}
private:
  std::string message_;
};

namespace utils
{
  class kernel_generation_stream;
}


enum leaf_t
{
  LHS_NODE_TYPE,
  PARENT_NODE_TYPE,
  RHS_NODE_TYPE
};

class mapped_object;
class template_base;

typedef std::pair<vcl_size_t, leaf_t> mapping_key;
typedef std::map<mapping_key, tools::shared_ptr<mapped_object> > mapping_type;


namespace tree_parsing
{

  template<class Fun>
  inline void traverse(scheduler::statement const & statement, vcl_size_t root_idx, Fun const & fun, bool inspect);

  inline void process(utils::kernel_generation_stream & stream, leaf_t leaf, std::string const & type_key, std::string const & to_process,
                      scheduler::statement const & statement, vcl_size_t root_idx, mapping_type const & mapping, std::set<std::string> & already_processed);
  inline std::string evaluate(leaf_t leaf, std::map<std::string, std::string> const & accessors, scheduler::statement const & statement, vcl_size_t root_idx,mapping_type const & mapping);
}

using scheduler::INT_TYPE;
using scheduler::UINT_TYPE;
using scheduler::ULONG_TYPE;
using scheduler::LONG_TYPE;
using scheduler::FLOAT_TYPE;
using scheduler::DOUBLE_TYPE;

typedef cl_uint vendor_id_type;
typedef cl_device_type device_type;
typedef std::string device_name_type;

class symbolic_binder
{
public:
  virtual ~symbolic_binder(){ }
  virtual bool bind(viennacl::backend::mem_handle const * ph) = 0;
  virtual unsigned int get(viennacl::backend::mem_handle const * ph) = 0;
};

class bind_to_handle : public symbolic_binder
{
public:
  bind_to_handle() : current_arg_(0){ }
  bool bind(viennacl::backend::mem_handle const * ph) {return (ph==NULL)?true:memory.insert(std::make_pair((void*)ph, current_arg_)).second; }
  unsigned int get(viennacl::backend::mem_handle const * ph){ return bind(ph) ? current_arg_++ : at(memory, (void*)ph); }
private:
  unsigned int current_arg_;
  std::map<void*,unsigned int> memory;
};

class bind_all_unique : public symbolic_binder
{
public:
  bind_all_unique() : current_arg_(0){ }
  bool bind(viennacl::backend::mem_handle const *) {return true; }
  unsigned int get(viennacl::backend::mem_handle const *){ return current_arg_++; }
private:
  unsigned int current_arg_;
  std::map<void*,unsigned int> memory;
};

enum binding_policy_t{
  BIND_ALL_UNIQUE,
  BIND_TO_HANDLE
};

inline tools::shared_ptr<symbolic_binder> make_binder(binding_policy_t policy)
{
  if (policy==BIND_TO_HANDLE)
    return tools::shared_ptr<symbolic_binder>(new bind_to_handle());
  else
    return tools::shared_ptr<symbolic_binder>(new bind_all_unique());
}

template<char C>
struct char_to_type{ };

class statements_container
{
public:
  typedef std::list<scheduler::statement> data_type;
  enum order_type { SEQUENTIAL, INDEPENDENT };

  statements_container(data_type const & data, order_type order) : data_(data), order_(order)
  { }

  statements_container(scheduler::statement const & s0) : order_(INDEPENDENT)
  {
    data_.push_back(s0);
  }

  statements_container(scheduler::statement const & s0, scheduler::statement const & s1, order_type order) : order_(order)
  {
    data_.push_back(s0);
    data_.push_back(s1);
  }

  std::list<scheduler::statement> const & data() const { return data_; }

  order_type order() const { return order_; }

private:
  std::list<scheduler::statement> data_;
  order_type order_;
};

}

}
#endif
