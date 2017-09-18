#ifndef VIENNACL_DEVICE_SPECIFIC_TEMPLATES_MATRIX_PRODUCT_HPP
#define VIENNACL_DEVICE_SPECIFIC_TEMPLATES_MATRIX_PRODUCT_HPP

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


/** @file viennacl/device_specific/templates/matrix_product_template.hpp
*
* Kernel template for the matrix product operation
*/

#include <vector>

#include "viennacl/scheduler/forwards.h"

#include "viennacl/detail/matrix_def.hpp"
#include "viennacl/matrix_proxy.hpp"

#include "viennacl/device_specific/templates/template_base.hpp"
#include "viennacl/device_specific/mapped_objects.hpp"
#include "viennacl/device_specific/utils.hpp"
#include "viennacl/device_specific/tree_parsing.hpp"
#include "viennacl/forwards.h"

#include "viennacl/tools/tools.hpp"

namespace viennacl
{
namespace device_specific
{

struct matrix_product_parameters : public template_base::parameters_type
{
  matrix_product_parameters(unsigned int simd_width
                            , unsigned int local_size_0, unsigned int KL, unsigned int local_size_1
                            , unsigned int ms, unsigned int ks, unsigned int ns
                            , fetching_policy_type A_fetching_policy_param, fetching_policy_type B_fetching_policy_param
                            , unsigned int local_fetch_0_param, unsigned int local_fetch_1_param): template_base::parameters_type(simd_width, local_size_0, local_size_1, 1),
    kL(KL), mS(ms), kS(ks), nS(ns), A_fetching_policy(A_fetching_policy_param), B_fetching_policy(B_fetching_policy_param),
    local_fetch_0(local_fetch_0_param), local_fetch_1(local_fetch_1_param),
    mL(ms*local_size_0), nL(ns*local_size_1){}

  unsigned int kL;

  unsigned int mS;
  unsigned int kS;
  unsigned int nS;

  fetching_policy_type A_fetching_policy;
  fetching_policy_type B_fetching_policy;

  unsigned int local_fetch_0;
  unsigned int local_fetch_1;

  unsigned int mL;
  unsigned int nL;
};

class matrix_product_template : public template_base_impl<matrix_product_template, matrix_product_parameters>
{

private:
  unsigned int n_lmem_elements() const
  {
    unsigned int N = 0;
    if (p_.A_fetching_policy==FETCH_FROM_LOCAL)
      N += p_.kL * (p_.mL+1);
    if (p_.B_fetching_policy==FETCH_FROM_LOCAL)
      N += p_.nL * (p_.kL+1);
    return N;
  }

  int check_invalid_impl(viennacl::ocl::device const & /*device*/) const
  {
    if (p_.A_fetching_policy!=FETCH_FROM_LOCAL && p_.B_fetching_policy!=FETCH_FROM_LOCAL&& (p_.local_fetch_0!=0 || p_.local_fetch_1!=0))
      return TEMPLATE_GLOBAL_MEMORY_REQUIRES_ZERO_LOCAL_FETCH;

    if ((p_.mS % p_.simd_width) > 0 || (p_.nS % p_.simd_width) > 0)
      return TEMPLATE_MS_NS_MUST_BE_SIMD_WIDTH_MULTIPLE;

    if (p_.kS > p_.kL)
      return TEMPLATE_KS_MUST_BE_SMALLER_THAN_KL;

    if (!(A_trans_=='N' && B_trans_=='T') && p_.simd_width>1)
      return TEMPLATE_SIMD_WIDTH_MUST_BE_ONE;

    if (p_.A_fetching_policy==FETCH_FROM_LOCAL || p_.B_fetching_policy==FETCH_FROM_LOCAL)
    {
      if ((p_.local_fetch_0*p_.local_fetch_1) !=(p_.local_size_0*p_.local_size_1))
        return TEMPLATE_LOCAL_FETCH_PRODUCT_MUST_MATCH_LOCAL_SIZE_PRODUCT;
    }

    if (p_.A_fetching_policy==FETCH_FROM_LOCAL)
    {
      unsigned int bound1 = (A_trans_=='N')?p_.kL:p_.mL;
      unsigned int bound0 = (A_trans_=='N')?p_.mL:p_.kL;

      if (p_.local_fetch_1>0 && (bound1 % p_.local_fetch_1)> 0)
        return A_trans_=='N'?TEMPLATE_LOCAL_FETCH_1_MUST_BE_KL_MULTIPLE:TEMPLATE_LOCAL_FETCH_1_MUST_BE_ML_MULTIPLE;

      if (p_.local_fetch_0>0 && (bound0 % (p_.local_fetch_0*p_.simd_width)) > 0)
        return A_trans_=='N'?TEMPLATE_LOCAL_FETCH_0_MUST_BE_NL_MULTIPLE:TEMPLATE_LOCAL_FETCH_0_MUST_BE_KL_MULTIPLE;

    }
    if (p_.B_fetching_policy==FETCH_FROM_LOCAL)
    {
      unsigned int bound1 = (B_trans_=='T')?p_.kL:p_.nL;
      unsigned int bound0 = (B_trans_=='T')?p_.nL:p_.kL;

      if (p_.local_fetch_1>0 && (bound1 % p_.local_fetch_1)> 0)
        return B_trans_=='T'?TEMPLATE_LOCAL_FETCH_1_MUST_BE_KL_MULTIPLE:TEMPLATE_LOCAL_FETCH_1_MUST_BE_ML_MULTIPLE;

      if (p_.local_fetch_0>0 && (bound0 % (p_.local_fetch_0*p_.simd_width)) > 0)
        return B_trans_=='T'?TEMPLATE_LOCAL_FETCH_1_MUST_BE_KL_MULTIPLE:TEMPLATE_LOCAL_FETCH_1_MUST_BE_ML_MULTIPLE;

    }

    return TEMPLATE_VALID;
  }

  static void parse(scheduler::statement const & s,
                    vcl_size_t & C_idx, leaf_t & C_leaf, vcl_size_t & alpha_idx, leaf_t & alpha_leaf,
                    vcl_size_t & A_idx, leaf_t & A_leaf, bool& A_trans, vcl_size_t & B_idx, leaf_t & B_leaf, bool& B_trans,
                    vcl_size_t & beta_idx, leaf_t & beta_leaf)
  {
    using namespace tree_parsing;
    using namespace scheduler;

    scheduler::statement::container_type const & array = s.array();
    vcl_size_t root_idx = s.root();

    C_idx = root_idx;
    C_leaf = LHS_NODE_TYPE;

    vcl_size_t node_add_idx = array[root_idx].rhs.node_index;

    vcl_size_t node_1_idx = array[node_add_idx].lhs.node_index;
    alpha_idx = node_1_idx;
    alpha_leaf = RHS_NODE_TYPE;

    vcl_size_t mat_prod_idx = array[node_1_idx].lhs.node_index;
    if (array[mat_prod_idx].lhs.type_family==MATRIX_TYPE_FAMILY)
    {
      A_trans = false;
      A_idx = mat_prod_idx;
    }
    else
    {
      A_trans = true;
      A_idx = array[mat_prod_idx].lhs.node_index;
    }
    A_leaf = LHS_NODE_TYPE;

    if (array[mat_prod_idx].rhs.type_family==MATRIX_TYPE_FAMILY)
    {
      B_trans = false;
      B_idx = mat_prod_idx;
      B_leaf = RHS_NODE_TYPE;
    }
    else
    {
      B_trans = true;
      B_idx = array[mat_prod_idx].rhs.node_index;
      B_leaf = LHS_NODE_TYPE;
    }

    vcl_size_t node_2_idx = array[node_add_idx].rhs.node_index;
    beta_idx = node_2_idx;
    beta_leaf = RHS_NODE_TYPE;
  }

  void VIENNACL_HANDLE_BOUNDS(bool fallback, utils::kernel_generation_stream & stream, std::string const & inbounds, std::string const & do_if, std::string do_else) const
  {
    if (fallback)
    {
      stream << "if (" << inbounds << ")" << std::endl;
      stream.inc_tab();
      stream << do_if << ";" << std::endl;
      stream.dec_tab();
      stream << "else" << std::endl;
      stream.inc_tab();
      stream << do_else << ";" << std::endl;
      stream.dec_tab();
    }
    else
      stream << do_if << ";" << std::endl;
  }


  std::string generate_impl(const std::string &kernel_prefix, const statements_container &statements, const std::vector<mapping_type> &mappings, bool fallback) const
  {
    using std::string;
    using tools::to_string;

    parameters_type pfallback(1, p_.local_size_0, p_.kL, p_.local_size_1, p_.mS, 1, p_.nS, p_.A_fetching_policy, p_.B_fetching_policy, p_.local_fetch_0, p_.local_fetch_1);
    parameters_type const & p = fallback?pfallback:p_;

#define VIENNACL_MUL_STRIDE1 string(fallback?"*#stride1":"")
#define VIENNACL_HANDLE_BOUNDS(in_bounds, to_load) (!fallback?string(to_load):string( string(in_bounds) + "?" + string(to_load) + ":0"))
#define VIENNACL_VSTORE(value, offset, ptr) vstore(p.simd_width, value, offset, ptr)

    string widthstr = tools::to_string(p.simd_width);

    //////////////////
    /// INIT
    /// //////////////
    utils::kernel_generation_stream stream;
    scheduler::statement const & st = statements.data().front();
    mapping_type const & mapping = mappings.front();

    bool A_trans = false, B_trans = false;
    vcl_size_t C_idx=0, alpha_idx=0, A_idx=0, B_idx=0, beta_idx=0;
    leaf_t C_leaf=LHS_NODE_TYPE, alpha_leaf=LHS_NODE_TYPE, A_leaf=LHS_NODE_TYPE, B_leaf=LHS_NODE_TYPE, beta_leaf=LHS_NODE_TYPE;
    parse(st, C_idx, C_leaf, alpha_idx, alpha_leaf, A_idx, A_leaf, A_trans, B_idx, B_leaf, B_trans, beta_idx, beta_leaf);

    mapped_matrix      * C     = (mapped_matrix*     )at(mapping, mapping_key(    C_idx,     C_leaf)).get();
    mapped_host_scalar * alpha = (mapped_host_scalar*)at(mapping, mapping_key(alpha_idx, alpha_leaf)).get();
    mapped_matrix      * A     = (mapped_matrix*     )at(mapping, mapping_key(    A_idx,     A_leaf)).get();
    mapped_matrix      * B     = (mapped_matrix*     )at(mapping, mapping_key(    B_idx,     B_leaf)).get();
    mapped_host_scalar * beta  = (mapped_host_scalar*)at(mapping, mapping_key( beta_idx,  beta_leaf)).get();

    //////////////////
    /// DECLARATIONS
    /// //////////////

    stream << " __attribute__((reqd_work_group_size(" << p.local_size_0 << "," << p.local_size_1 << ",1)))" << std::endl;
    std::map<std::string, unsigned int> widths;
    widths[A->name()] = p.simd_width;
    widths[B->name()] = p.simd_width;
    generate_prototype(stream, kernel_prefix, "unsigned int M, unsigned int N, unsigned int K, ", mappings, statements, widths);
    stream << "{" << std::endl;
    stream.inc_tab();
    if(!fallback)
    {
      stream << A->process("#start1 /= " + to_string(p.simd_width) + ";") << std::endl;
      stream << A->process("#ld /= " + to_string(p.simd_width) + ";") << std::endl;
      stream << B->process("#start1/= "  + to_string(p.simd_width) + ";") << std::endl;
      stream << B->process("#ld /= " + to_string(p.simd_width) + ";") << std::endl;
    }
    tree_parsing::process(stream, PARENT_NODE_TYPE, "matrix", "#pointer += $OFFSET{#start1, #start2};", statements, mappings);
    tree_parsing::process(stream, PARENT_NODE_TYPE, "matrix", "#ld *= #nldstride;", statements, mappings);

    ///Result Values
    stream << C->process("#scalartype rC[" + to_string(p.mS) + "][" + to_string(p.nS) + "] = {{(#scalartype)0}};") << std::endl;
    if (p.A_fetching_policy==FETCH_FROM_LOCAL)
      stream << A->process("#scalartype rA[" + to_string(p.kS) + "][" + to_string(p.mS) + "];") << std::endl;
    else
      stream << A->process(utils::append_width("#scalartype",p.simd_width) + " rA[" + to_string(p.kS) + "][" + to_string(p.mS/p.simd_width) + "];") << std::endl;
    if (p.B_fetching_policy==FETCH_FROM_LOCAL)
      stream << B->process("#scalartype rB[" + to_string(p.kS) + "][" + to_string(p.nS) + "];");
    else
      stream << B->process(utils::append_width("#scalartype",p.simd_width) + " rB[" + to_string(p.kS) + "][" + to_string(p.nS/p.simd_width) + "];") << std::endl;


    if (p.A_fetching_policy==FETCH_FROM_LOCAL)
      stream << A->process("__local #scalartype lA[" + to_string(p.kL*(p.mL+1)) + "];");
    if (p.B_fetching_policy==FETCH_FROM_LOCAL)
      stream << B->process("__local #scalartype lB[" + to_string(p.kL*(p.nL+1)) + "];");
    stream << std::endl;

    stream << "size_t gidx = get_group_id(0);" << std::endl;
    stream << "size_t gidy = get_group_id(1);" << std::endl;
    stream << "size_t idx = get_local_id(0);" << std::endl;
    stream << "size_t idy = get_local_id(1);" << std::endl;

    if (p.A_fetching_policy==FETCH_FROM_LOCAL || p.B_fetching_policy==FETCH_FROM_LOCAL)
    {
      stream << std::endl;
      stream << "size_t idt = " << p.local_size_0 << "*idy + idx;" << std::endl;
      stream << "size_t idxT = idt % " << p.local_fetch_0 << ";" << std::endl;
      stream << "size_t idyT = idt / " << p.local_fetch_0 << ";" << std::endl;
    }
    stream << std::endl;

    if (fallback)
    {
      //Bounds checking for M (in A, C)
      stream << "bool in_bounds_m[" << p.mS << "];" << std::endl;
      stream << "for(size_t m = 0; m < " << p.mS << "; m++)" << std::endl;
      stream.inc_tab();
      switch (p.A_fetching_policy)
      {
      case FETCH_FROM_GLOBAL_CONTIGUOUS:
        stream << "in_bounds_m[m] = gidx*" << p.mL << " + idx*" << p.mS << " + m < M;" << std::endl;
        break;
      default:
        stream << "in_bounds_m[m] = gidx*" << p.mL << " + idx + m*" << p.local_size_0 << " < M;" << std::endl;
        break;
      }
      stream.dec_tab();

      //Bounds checking for A if Local
      if (p.A_fetching_policy==FETCH_FROM_LOCAL)
      {
        unsigned int fetch_size = (A_trans_=='N'?p.local_fetch_0*p.simd_width:p.local_fetch_1);
        stream << "bool in_bounds_m_local[" << p.mL/fetch_size << "];" << std::endl;
        stream << "for(size_t m = 0; m < " << p.mL/fetch_size << "; m++)" << std::endl;
        stream.inc_tab();
        stream << "in_bounds_m_local[m] = gidx*" << p.mL << " + " << (A_trans_=='N'?"idxT":"idyT") << " + m*" << fetch_size << " < M;" << std::endl;
        stream.dec_tab();
      }

      //Bounds checking for N (in B, C)
      stream << "bool in_bounds_n[" << p.nS << "];" << std::endl;
      stream << "for(size_t n = 0; n < " << p.nS << "; n++)" << std::endl;
      stream.inc_tab();
      switch (p.B_fetching_policy)
      {
      case FETCH_FROM_GLOBAL_CONTIGUOUS:
        stream << "in_bounds_n[n] = gidy*" << p.nL << " + idy*" << p.nS << " + n < N;" << std::endl;
        break;
      default:
        stream << "in_bounds_n[n] = gidy*" << p.nL << " + idy + n*" << p.local_size_1 << " < N;" << std::endl;
        break;
      }
      stream.dec_tab();

      //Bounds checking for B if Local
      if (p.B_fetching_policy==FETCH_FROM_LOCAL)
      {
        unsigned int fetch_size = (B_trans_=='T'?p.local_fetch_0*p.simd_width:p.local_fetch_1);
        stream << "bool in_bounds_n_local[" << p.nL/fetch_size << "];" << std::endl;
        stream << "for(size_t n = 0; n < " <<  p.nL/fetch_size << "; n++)" << std::endl;
        stream.inc_tab();
        stream << "in_bounds_n_local[n] = gidy*" << p.nL << " + " << (B_trans_=='T'?"idxT":"idyT") << " + n*" << fetch_size << " < N;" << std::endl;
        stream.dec_tab();
      }
    }

    switch (p.A_fetching_policy)
    {
    case FETCH_FROM_LOCAL:
      if (A_trans_=='N')
        stream << A->process("#pointer += (gidx*" + to_string(p.mL/p.simd_width) + " + idxT)" + VIENNACL_MUL_STRIDE1 + " + idyT*#ld;") << std::endl;
      else
        stream << A->process("#pointer += idxT" + VIENNACL_MUL_STRIDE1 + " + gidx*" + to_string(p.mL/p.simd_width) + "*#ld + idyT*#ld;") << std::endl;
      break;

    case FETCH_FROM_GLOBAL_CONTIGUOUS:
      if (A_trans_=='N')
        stream << A->process("#pointer += (gidx*" + to_string(p.mL/p.simd_width) + "+ idx*" + to_string(p.mS/p.simd_width) + ")" + VIENNACL_MUL_STRIDE1 + ";") << std::endl;
      else
        stream << A->process("#pointer += (gidx*" + to_string(p.mL/p.simd_width) + "+ idx*" + to_string(p.mS/p.simd_width) + ")*#ld;") << std::endl;
      break;

    case FETCH_FROM_GLOBAL_STRIDED:
      if (A_trans_=='N')
        stream << A->process("#pointer += (gidx*" + to_string(p.mL/p.simd_width) + "+ idx" + ")" + VIENNACL_MUL_STRIDE1 + ";") << std::endl;
      else
        stream << A->process("#pointer += (gidx*" + to_string(p.mL/p.simd_width) + "+ idx)*#ld;") << std::endl;
      break;

    //default: break;
    }

    switch (p.B_fetching_policy)
    {
    case FETCH_FROM_LOCAL:
      if (B_trans_=='T')
        stream << B->process("#pointer += (gidy*" + to_string(p.nL/p.simd_width) + " + idxT" + ")" + VIENNACL_MUL_STRIDE1 + " + idyT*#ld;") << std::endl;
      else
        stream << B->process("#pointer += idxT" + VIENNACL_MUL_STRIDE1 + " + gidy*" + to_string(p.nL/p.simd_width) + "*#ld + idyT*#ld;") << std::endl;
      break;

    case FETCH_FROM_GLOBAL_CONTIGUOUS:
      if (B_trans_=='T')
        stream << B->process("#pointer += (gidy*" + to_string(p.nL/p.simd_width) + "+ idy*" + to_string(p.nS/p.simd_width) + ")" + VIENNACL_MUL_STRIDE1 + ";") << std::endl;
      else
        stream << B->process("#pointer += (gidy*" + to_string(p.nL/p.simd_width) + "+ idy*" + to_string(p.nS/p.simd_width) + ")*#ld;") << std::endl;
      break;

    case FETCH_FROM_GLOBAL_STRIDED:
      if (B_trans_=='T')
        stream << B->process("#pointer += (gidy*" + to_string(p.nL/p.simd_width) + "+ idy" + ")" + VIENNACL_MUL_STRIDE1 + ";") << std::endl;
      else
        stream << B->process("#pointer += (gidy*" + to_string(p.nL/p.simd_width) + "+ idy)*#ld;") << std::endl;
      break;

    //default: break;
    }

    stream << std::endl;
    stream << "size_t K_size_t = K;" << std::endl;
    stream << "for(size_t block_k=0; block_k < K_size_t; block_k+=" << p.kL << "){" << std::endl;
    stream.inc_tab();

    if (p.A_fetching_policy==FETCH_FROM_LOCAL)
    {
      if (A_trans_=='N')
        stream << A->process("__local #scalartype* plA = lA + idyT*" + to_string(p.mL + 1) + " + " + to_string(p.simd_width) + "*idxT;") << std::endl;
      else
        stream << A->process("__local #scalartype* plA = lA + idxT*" + to_string(p.mL + 1) + " + idyT;") << std::endl;
    }


    if (p.B_fetching_policy==FETCH_FROM_LOCAL)
    {
      if (B_trans_=='T')
        stream  << B->process("__local #scalartype* plB = lB + idyT*" + to_string(p.nL+1) + " + " + to_string(p.simd_width) + "*idxT;") << std::endl;
      else
        stream << B->process("__local #scalartype* plB = lB + idxT*" + to_string(p.nL+1) + "+ idyT;") <<std::endl;
    }


    if (p.A_fetching_policy==FETCH_FROM_LOCAL || p.B_fetching_policy==FETCH_FROM_LOCAL)
      stream << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;

    ///Fetch LHS to Local Memory
    if (p.A_fetching_policy==FETCH_FROM_LOCAL && A_trans_=='N')
      for (unsigned int k = 0; k < p.kL; k += p.local_fetch_1)
        for (unsigned int m = 0; m < p.mL; m += p.local_fetch_0*p.simd_width)
        {
          string in_bounds = "in_bounds_m_local[" + to_string(m/(p.local_fetch_0*p.simd_width)) + "]";
          string to_load = "#pointer[" + to_string(k) + "*#ld + " + to_string(m/p.simd_width) + VIENNACL_MUL_STRIDE1 + "]";
          stream << A->process(VIENNACL_VSTORE(VIENNACL_HANDLE_BOUNDS(in_bounds, to_load), "0", "plA + " + to_string(k*(p.mL+1)+m))) << ";" << std::endl;
        }
    else if (p.A_fetching_policy==FETCH_FROM_LOCAL && A_trans_=='T')
      for (unsigned int k = 0; k < p.mL; k += p.local_fetch_1)
        for (unsigned int m = 0; m < p.kL; m += p.local_fetch_0*p.simd_width)
        {
          string in_bounds = "in_bounds_m_local[" + to_string(k/p.local_fetch_1) + "]";
          string to_load = "#pointer[" + to_string(k) + "*#ld + " + to_string(m/p.simd_width) + VIENNACL_MUL_STRIDE1 + "]";
          stream << A->process(VIENNACL_VSTORE(VIENNACL_HANDLE_BOUNDS(in_bounds, to_load), "0", "plA + " + to_string(m*(p.mL+1)+k))) << ";" << std::endl;
        }

    if (p.B_fetching_policy==FETCH_FROM_LOCAL && B_trans_=='T')
      for (unsigned int k = 0; k < p.kL; k += p.local_fetch_1)
        for (unsigned int n = 0; n < p.nL; n += p.local_fetch_0*p.simd_width)
        {
          string in_bounds = "in_bounds_n_local[" + to_string(n/(p.local_fetch_0*p.simd_width)) + "]";
          string to_load = "#pointer[" + to_string(k) + "*#ld + " + to_string(n/p.simd_width) + VIENNACL_MUL_STRIDE1 + "]";
          stream << B->process(VIENNACL_VSTORE(VIENNACL_HANDLE_BOUNDS(in_bounds, to_load), "0", "plB + " + to_string(k*(p.nL+1)+n))) << ";" << std::endl;
        }
    else if (p.B_fetching_policy==FETCH_FROM_LOCAL && B_trans_=='N')
      for (unsigned int k = 0; k < p.nL; k += p.local_fetch_1)
        for (unsigned int n = 0; n < p.kL; n += p.local_fetch_0*p.simd_width)
        {
          string in_bounds = "in_bounds_n_local[" + to_string(k/p.local_fetch_1) + "]";
          string to_load = "#pointer[" + to_string(k) + "*#ld + " + to_string(n/p.simd_width) + VIENNACL_MUL_STRIDE1 + "]";
          stream << B->process(VIENNACL_VSTORE(VIENNACL_HANDLE_BOUNDS(in_bounds, to_load), "0", "plB + " + to_string(n*(p.nL+1)+k))) << ";" << std::endl;
        }

    if (p.A_fetching_policy==FETCH_FROM_LOCAL || p.B_fetching_policy == FETCH_FROM_LOCAL)
    {
      stream << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      stream << "size_t offA = " << p.simd_width << "*idx;" << std::endl;
      stream << "size_t offB = " << p.simd_width << "*idy;" << std::endl;
    }

    if (fallback)
      stream << "for(size_t k = 0; k < " << p.kL << " && (block_k + k < K_size_t); k+=" << p.kS << "){" << std::endl;
    else
      stream << "for(size_t k = 0; k < " << p.kL << "; k+=" << p.kS << "){" << std::endl;
    stream.inc_tab();

    ///Fetch LHS to registers
    stream << "#pragma unroll " << p.kS <<  std::endl;
    stream << "for(size_t kk = 0; kk < " << p.kS << "; kk++)" << std::endl;
    stream << "#pragma unroll " << p.mS/p.simd_width << std::endl;
    stream << "for(size_t mm = 0; mm < " << p.mS/p.simd_width << "; mm++)" << std::endl;
    stream << "{" << std::endl;
    stream.inc_tab();
    switch (p.A_fetching_policy)
    {
    case FETCH_FROM_LOCAL:
      for (unsigned int ss = 0; ss < p.simd_width; ++ss)
        stream << "rA[kk][mm*" << p.simd_width << "+" << ss << "] = lA[offA + mm*" << p.local_size_0*p.simd_width << "+" << ss << "+ kk*" << (p.mL+1) << "];" << std::endl;
      break;

    case FETCH_FROM_GLOBAL_CONTIGUOUS:
    {
      if (A_trans_=='N')
        stream << "rA[kk][mm] = " << A->process(VIENNACL_HANDLE_BOUNDS("in_bounds_m[mm]", "#pointer[kk*#ld + mm" + VIENNACL_MUL_STRIDE1 + "]")) << ";" << std::endl;
      else
        stream << "rA[kk][mm] = " << A->process(VIENNACL_HANDLE_BOUNDS("in_bounds_m[mm]", "#pointer[mm*#ld + kk" + VIENNACL_MUL_STRIDE1 + "]")) << ";" << std::endl;
      break;
    }

    case FETCH_FROM_GLOBAL_STRIDED:
    {
      if (A_trans_=='N')
        stream << "rA[kk][mm] = " << A->process(VIENNACL_HANDLE_BOUNDS("in_bounds_m[mm]", "#pointer[kk*#ld + mm*" + to_string(p.local_size_0) + VIENNACL_MUL_STRIDE1 + "]")) << ";" << std::endl;
      else
        stream << "rA[kk][mm] = " << A->process(VIENNACL_HANDLE_BOUNDS("in_bounds_m[mm]", "#pointer[mm*#ld*" + to_string(p.local_size_0) + " + kk" + VIENNACL_MUL_STRIDE1 + "]")) << ";" << std::endl;
      break;
    }

    //default: break;
    }
    stream.dec_tab();
    stream << "}" << std::endl;

    stream << "#pragma unroll " << p.kS << std::endl;
    stream << "for(size_t kk = 0; kk < " << p.kS << "; kk++)" << std::endl;
    stream << "#pragma unroll " << p.nS/p.simd_width << std::endl;
    stream << "for(size_t nn = 0; nn < " << p.nS/p.simd_width << "; nn++)" << std::endl;
    stream << "{" << std::endl;
    stream.inc_tab();
    switch (p.B_fetching_policy)
    {
    case FETCH_FROM_LOCAL:
      for (unsigned int ss = 0; ss < p.simd_width; ++ss)
        stream << "rB[kk][nn*" << p.simd_width << "+" << ss << "] = lB[offB + nn*" << p.local_size_1*p.simd_width << "+" << ss  << "+ kk*" << (p.nL+1) << "];" << std::endl;
      break;

    case FETCH_FROM_GLOBAL_CONTIGUOUS:
    {
      if (B_trans_=='T')
        stream << "rB[kk][nn] = " << B->process(VIENNACL_HANDLE_BOUNDS("in_bounds_n[nn]", "#pointer[kk*#ld + nn" + VIENNACL_MUL_STRIDE1 + "]")) << ";" << std::endl;
      else
        stream << "rB[kk][nn] = " << B->process(VIENNACL_HANDLE_BOUNDS("in_bounds_n[nn]", "#pointer[nn*#ld + kk" + VIENNACL_MUL_STRIDE1 + "]")) << ";" << std::endl;
      break;
    }

    case FETCH_FROM_GLOBAL_STRIDED:
    {
      if (B_trans_=='T')
        stream << "rB[kk][nn] = " << B->process(VIENNACL_HANDLE_BOUNDS("in_bounds_n[nn]", "#pointer[kk*#ld + nn*" + to_string(p.local_size_1) + VIENNACL_MUL_STRIDE1 + "]")) << ";" << std::endl;
      else
        stream << "rB[kk][nn] = " << B->process(VIENNACL_HANDLE_BOUNDS("in_bounds_n[nn]", "#pointer[nn*#ld*" + to_string(p.local_size_1) + " + kk" + VIENNACL_MUL_STRIDE1 + "]")) << ";" << std::endl;
      break;
    }

    //default: break;
    }
    stream.dec_tab();
    stream << "}" << std::endl;


    ///Increment pointers
    switch (p.A_fetching_policy)
    {
    case FETCH_FROM_LOCAL:
      stream << "offA += " << p.kS*(p.mL+1) << ";" << std::endl;
      break;

    default:
      if (A_trans_=='N')
        stream << A->process("#pointer += " + to_string(p.kS) + "*#ld;") << std::endl;
      else
        stream << A->process("#pointer += " + to_string(p.kS) + "" + VIENNACL_MUL_STRIDE1 + ";") << std::endl;
      break;
    }


    switch (p.B_fetching_policy)
    {
    case FETCH_FROM_LOCAL:
      stream << "offB += " << p.kS*(p.nL+1) << ";" << std::endl;
      break;

    default:
      if (B_trans_=='T')
        stream << B->process("#pointer += " + to_string(p.kS) + "*#ld;") << std::endl;
      else
        stream << B->process("#pointer += " + to_string(p.kS) + "" + VIENNACL_MUL_STRIDE1 + ";") << std::endl;
      break;
    }


    stream << "#pragma unroll " << p.kS << std::endl;
    stream << "for(size_t kk = 0; kk <" << p.kS << "; ++kk)" << std::endl;
    stream << "{" << std::endl;
    stream.inc_tab();
    for (unsigned int nn=0; nn < p.nS; ++nn)
      for (unsigned int mm=0; mm < p.mS; ++mm)
      {
        string res_str, lhs_str, rhs_str;
        res_str = "rC[" + tools::to_string(mm) + "][" + tools::to_string(nn) + "]";
        if (p.A_fetching_policy==FETCH_FROM_LOCAL || p.simd_width==1)
          lhs_str = "rA[kk][" + tools::to_string(mm) + "]";
        else
          lhs_str = "rA[kk][" + tools::to_string(mm/p.simd_width) + "].s" + tools::to_string(mm%p.simd_width);
        if (p.B_fetching_policy==FETCH_FROM_LOCAL || p.simd_width==1)
          rhs_str = "rB[kk]["+tools::to_string(nn)+"]";
        else
          rhs_str = "rB[kk]["+tools::to_string(nn/p.simd_width)+"].s"+tools::to_string(nn%p.simd_width);
        stream << res_str << "=" << "fma(" << lhs_str << "," << rhs_str << "," << res_str << ");" << std::endl;
      }
    stream.dec_tab();
    stream << "}" << std::endl;




    stream.dec_tab();
    stream << "}" << std::endl;

    //Increment global pointer if local memory is used
    //Else, it's incremented directly when fetching
    if (p.A_fetching_policy==FETCH_FROM_LOCAL)
    {
      if (A_trans_=='N')
        stream << A->process("#pointer += " + to_string(p.kL) + "*#ld;") << std::endl;
      else
        stream << A->process("#pointer += " + to_string(p.kL) + "" + VIENNACL_MUL_STRIDE1 + ";") << std::endl;
    }

    if (p.B_fetching_policy==FETCH_FROM_LOCAL)
    {
      if (B_trans_=='T')
        stream << B->process("#pointer += " + to_string(p.kL) + "*#ld;") << std::endl;
      else
        stream << B->process("#pointer += " + to_string(p.kL) + "" + VIENNACL_MUL_STRIDE1 + ";") << std::endl;
    }

    stream.dec_tab();
    stream << "}" << std::endl;


    if (C->row_major())
    {
      unsigned int ministartstride0 = p.A_fetching_policy==FETCH_FROM_GLOBAL_CONTIGUOUS?p.mS:p.simd_width;
      unsigned int ministartstride1 = p.B_fetching_policy==FETCH_FROM_GLOBAL_CONTIGUOUS?p.nS:p.simd_width;

      stream << C->process("#pointer += gidx*" + to_string(p.mL) + "*#ld;") << std::endl;
      stream << C->process("#pointer += idx*" + to_string(ministartstride0) + "*#ld;") << std::endl;
      stream << C->process("#pointer += gidy*" + to_string(p.nL) + "*#stride2;") << std::endl;
      stream << C->process("#pointer += idy*" + to_string(ministartstride1) + "*#stride2;") << std::endl;

      for (unsigned int n=0; n < p.nS; ++n)
      {
        for (unsigned int m=0; m < p.mS; ++m)
        {
          unsigned int ministride1 = p.A_fetching_policy==FETCH_FROM_GLOBAL_CONTIGUOUS?1:p.local_size_0;
          string Cj = to_string((m/p.simd_width)*(ministride1*p.simd_width) + m%p.simd_width);
          if (fallback)
          {
            stream << "if (in_bounds_m[" + to_string(m) + "] && in_bounds_n[" + to_string(n) + "])" << std::endl;
            stream.inc_tab();
          }
          stream << C->process("#pointer[" + Cj + "*#ld] = rC[" + to_string(m) + "][" + to_string(n) + "]*" + alpha->name() + "+ #pointer[" + Cj + "*#ld]*" + beta->name() + ";") << std::endl;
          if (fallback)
            stream.dec_tab();
        }
        if ((n+1)%p.simd_width>0 || p.B_fetching_policy==FETCH_FROM_GLOBAL_CONTIGUOUS)
          stream << C->process("#pointer += #stride2;") << std::endl;
        else
          stream << C->process("#pointer += " + to_string((p.local_size_1*p.simd_width) - (p.simd_width-1)) + "*#stride2;") << std::endl;
      }

    }
    else
    {
      unsigned int ministartstride0 = p.A_fetching_policy==FETCH_FROM_GLOBAL_CONTIGUOUS?p.mS:p.simd_width;
      unsigned int ministartstride1 = p.B_fetching_policy==FETCH_FROM_GLOBAL_CONTIGUOUS?p.nS:p.simd_width;

      stream << C->process("#pointer += gidx*" + to_string(p.mL) + "*#stride1;") << std::endl;
      stream << C->process("#pointer += idx*" + to_string(ministartstride0) + "*#stride1;") << std::endl;
      stream << C->process("#pointer += gidy*" + to_string(p.nL) + "*#ld;") << std::endl;
      stream << C->process("#pointer += idy*" + to_string(ministartstride1) + "*#ld;") << std::endl;

      for (unsigned int m=0; m < p.mS; ++m)
      {
        for (unsigned int n=0; n < p.nS; ++n)
        {
          unsigned int ministride1 = p.B_fetching_policy==FETCH_FROM_GLOBAL_CONTIGUOUS?1:p.local_size_1;
          string Cj = to_string((n/p.simd_width)*(ministride1*p.simd_width) + n%p.simd_width);
          if (fallback)
          {
            stream << "if (in_bounds_m[" + to_string(m) + "] && in_bounds_n[" + to_string(n) + "])" << std::endl;
            stream.inc_tab();
          }
          stream << C->process("#pointer[" + Cj + "*#ld] = rC[" + to_string(m) + "][" + to_string(n) + "]*" + alpha->name() + " + #pointer[" + Cj + "*#ld]*" + beta->name() + ";") << std::endl;
          if (fallback)
            stream.dec_tab();
        }

        if ((m+1)%p.simd_width>0 || p.A_fetching_policy==FETCH_FROM_GLOBAL_CONTIGUOUS)
          stream << C->process("#pointer += #stride1;") << std::endl;
        else
          stream << C->process("#pointer += " + to_string((p.local_size_0*p.simd_width) - (p.simd_width-1)) + "*#stride1;") << std::endl;
      }
    }

    stream.dec_tab();
    stream << "}" << std::endl;

    return stream.str();

#undef VIENNACL_MUL_STRIDE1
#undef VIENNACL_HANDLE_BOUNDS
#undef VIENNACL_VSTORE
  }

  std::vector<std::string> generate_impl(std::string const & kernel_prefix, statements_container const & statements, std::vector<mapping_type> const & mappings) const
  {
    std::vector<std::string> res;
    res.push_back(generate_impl(kernel_prefix, statements, mappings, false));
    res.push_back(generate_impl(kernel_prefix, statements, mappings, true));
    return res;
  }

  template<class NumericT>
  void enqueue_block(scheduler::statement & statement,
                     scheduler::lhs_rhs_element& eA, scheduler::lhs_rhs_element& eB, scheduler::lhs_rhs_element& eC, scheduler::lhs_rhs_element& ebeta,
                     matrix_base<NumericT> const & A, matrix_base<NumericT> const & B, matrix_base<NumericT> const & C, NumericT beta,
                     std::vector<lazy_program_compiler> & programs, std::string const & kernel_prefix, vcl_size_t id)
  {
    if (A.size1()==0 || A.size2()==0 || B.size1()==0 || B.size2()==0 || C.size1()==0 || C.size2()==0)
      return;

    viennacl::ocl::kernel& kernel = programs[id].program().get_kernel(kernel_prefix);

    kernel.local_work_size(0, p_.local_size_0);
    kernel.local_work_size(1, p_.local_size_1);

    scheduler::statement::assign_element(eA, A);
    scheduler::statement::assign_element(eB, B);
    scheduler::statement::assign_element(eC, C);
    scheduler::statement::assign_element(ebeta, beta);

    if (id==1)
    {
      kernel.global_work_size(0, tools::align_to_multiple(tools::align_to_multiple((unsigned int)C.size1(),p_.mS)/p_.mS, p_.local_size_0));
      kernel.global_work_size(1, tools::align_to_multiple(tools::align_to_multiple((unsigned int)C.size2(),p_.nS)/p_.nS, p_.local_size_1));
    }
    else
    {
      kernel.global_work_size(0, C.size1()/p_.mS);
      kernel.global_work_size(1, C.size2()/p_.nS);
    }
    unsigned int current_arg = 0;
    kernel.arg(current_arg++, cl_uint(C.size1()));
    kernel.arg(current_arg++, cl_uint(C.size2()));
    if (A.row_major())
      kernel.arg(current_arg++, cl_uint(A_trans_=='T'?A.size2():A.size1()));
    else
      kernel.arg(current_arg++, cl_uint(A_trans_=='N'?A.size2():A.size1()));
    set_arguments(statement, kernel, current_arg);
    viennacl::ocl::enqueue(kernel);

  }

  template<class NumericT>
  matrix_slice< viennacl::matrix_base<NumericT> >  create_slice(viennacl::matrix_base<NumericT>* scheduler::lhs_rhs_element::*ptr, scheduler::lhs_rhs_element const & element,
                                                                          vcl_size_t s0_0, vcl_size_t s0_1, vcl_size_t s1_0, vcl_size_t s1_1, bool swap)
  {
    matrix_base<NumericT> & M = *(element.*ptr);
    slice s0(s0_0, 1, s0_1 - s0_0);
    slice s1(s1_0, 1, s1_1 - s1_0);
    if (swap)
      std::swap(s0, s1);
    return matrix_slice<viennacl::matrix_base<NumericT> >(M, s0, s1);
  }

  template<class NumericT>
  void enqueue_impl(viennacl::matrix_base<NumericT>* scheduler::lhs_rhs_element::*ptr_matrix,
                    scheduler::statement & statement, scheduler::lhs_rhs_element & A, scheduler::lhs_rhs_element & B, scheduler::lhs_rhs_element & C, scheduler::lhs_rhs_element & beta,
                    NumericT beta_value, std::vector<lazy_program_compiler> & programs, std::string const & kernel_prefix)
  {
    using namespace device_specific::utils;
    vcl_size_t ldstrideA = call_on_matrix(A, leading_stride());
    vcl_size_t ldstrideB = call_on_matrix(B, leading_stride());
    vcl_size_t ldstrideC = call_on_matrix(C, leading_stride());
    vcl_size_t ldstartA = call_on_matrix(A, leading_start());
    vcl_size_t ldstartB = call_on_matrix(B, leading_start());
    bool swap_A = ((A_trans_=='T') ^ utils::call_on_matrix(A, row_major_fun()));
    bool swap_B = ((B_trans_=='T') ^ utils::call_on_matrix(B, row_major_fun()));

    vcl_size_t M = call_on_matrix(C, size1_fun());
    vcl_size_t N = call_on_matrix(C, size2_fun());
    vcl_size_t K;
    if (utils::call_on_matrix(A, row_major_fun()))
      K = A_trans_=='T'?call_on_matrix(A, size2_fun()):call_on_matrix(A, size1_fun());
    else
      K = A_trans_=='N'?call_on_matrix(A, size2_fun()):call_on_matrix(A, size1_fun());

    if (M < p_.mL || N < p_.nL || K < p_.kL || ldstrideA> 1 || ldstrideB > 1 || ldstrideC > 1 ||
        (p_.simd_width>1 && (ldstartA % p_.simd_width > 0 || ldstartB % p_.simd_width > 0)))
    {
      enqueue_block(statement, A, B, C, beta, create_slice(ptr_matrix, A, 0, M, 0, K, swap_A),
                    create_slice(ptr_matrix, B, 0, K, 0, N,  swap_B),
                    create_slice(ptr_matrix, C, 0, M, 0, N, false), beta_value, programs, kernel_prefix, 1);
      return;
    }


    scheduler::lhs_rhs_element Acopy = A;
    scheduler::lhs_rhs_element Bcopy = B;
    scheduler::lhs_rhs_element Ccopy = C;

    vcl_size_t lM = M / p_.mL * p_.mL;
    vcl_size_t lN = N / p_.nL * p_.nL;
    vcl_size_t lK = K / p_.kL * p_.kL;


    enqueue_block(statement, A, B, C, beta, create_slice<NumericT>(ptr_matrix, Acopy, 0, lM, 0, lK, swap_A), create_slice<NumericT>(ptr_matrix, Bcopy, 0, lK, 0, lN, swap_B), create_slice<NumericT>(ptr_matrix, Ccopy, 0, lM, 0, lN, false), beta_value, programs, kernel_prefix, 0);
    enqueue_block(statement, A, B, C, beta, create_slice<NumericT>(ptr_matrix, Acopy, 0, lM, lK, K, swap_A), create_slice<NumericT>(ptr_matrix, Bcopy, lK, K, 0, lN, swap_B), create_slice<NumericT>(ptr_matrix, Ccopy, 0, lM, 0, lN, false), (NumericT)1, programs, kernel_prefix, 1);

    enqueue_block(statement, A, B, C, beta, create_slice<NumericT>(ptr_matrix, Acopy, 0, lM, 0, lK, swap_A), create_slice<NumericT>(ptr_matrix, Bcopy, 0, lK, lN, N, swap_B), create_slice<NumericT>(ptr_matrix, Ccopy, 0, lM, lN, N, false), beta_value, programs, kernel_prefix, 1);
    enqueue_block(statement, A, B, C, beta, create_slice<NumericT>(ptr_matrix, Acopy, 0, lM, lK, K, swap_A), create_slice<NumericT>(ptr_matrix, Bcopy, lK, K, lN, N, swap_B), create_slice<NumericT>(ptr_matrix, Ccopy, 0, lM, lN, N, false), (NumericT)1, programs, kernel_prefix, 1);

    enqueue_block(statement, A, B, C, beta, create_slice<NumericT>(ptr_matrix, Acopy, lM, M, 0, lK, swap_A), create_slice<NumericT>(ptr_matrix, Bcopy, 0, lK, 0, lN, swap_B), create_slice<NumericT>(ptr_matrix, Ccopy, lM, M, 0, lN, false), beta_value, programs, kernel_prefix, 1);
    enqueue_block(statement, A, B, C, beta, create_slice<NumericT>(ptr_matrix, Acopy, lM, M, lK, K, swap_A), create_slice<NumericT>(ptr_matrix, Bcopy, lK, K, 0, lN, swap_B), create_slice<NumericT>(ptr_matrix, Ccopy, lM, M, 0, lN, false), (NumericT)1, programs, kernel_prefix, 1);

    enqueue_block(statement, A, B, C, beta, create_slice<NumericT>(ptr_matrix, Acopy, lM, M, 0, lK, swap_A), create_slice<NumericT>(ptr_matrix, Bcopy, 0, lK, lN, N, swap_B), create_slice<NumericT>(ptr_matrix, Ccopy, lM, M, lN, N, false), beta_value, programs, kernel_prefix, 1);
    enqueue_block(statement, A, B, C, beta, create_slice<NumericT>(ptr_matrix, Acopy, lM, M, lK, K, swap_A), create_slice<NumericT>(ptr_matrix, Bcopy, lK, K, lN, N, swap_B), create_slice<NumericT>(ptr_matrix, Ccopy, lM, M, lN, N, false), (NumericT)1, programs, kernel_prefix, 1);
  }

public:
  matrix_product_template(matrix_product_template::parameters_type const & parameters, char A_trans, char B_trans) : template_base_impl<matrix_product_template, matrix_product_parameters>(parameters, BIND_ALL_UNIQUE), A_trans_(A_trans), B_trans_(B_trans){ }

  virtual void enqueue(std::string const & kernel_prefix, std::vector<lazy_program_compiler> & programs, statements_container const & statements)
  {
    using namespace device_specific::utils;
    using namespace tree_parsing;

    scheduler::statement const & st = statements.data().front();
    bool A_trans, B_trans;
    vcl_size_t C_idx=0, A_idx=0, B_idx=0, alpha_idx=0, beta_idx = 0;
    leaf_t C_leaf=LHS_NODE_TYPE, A_leaf=LHS_NODE_TYPE, B_leaf=LHS_NODE_TYPE, alpha_leaf=LHS_NODE_TYPE, beta_leaf=LHS_NODE_TYPE;
    parse(st, C_idx, C_leaf, alpha_idx, alpha_leaf, A_idx, A_leaf, A_trans, B_idx, B_leaf, B_trans, beta_idx, beta_leaf);

    scheduler::statement stcopy = st;
    scheduler::lhs_rhs_element& A = utils::lhs_rhs_element(stcopy, A_idx, A_leaf);
    scheduler::lhs_rhs_element& B = utils::lhs_rhs_element(stcopy, B_idx, B_leaf);
    scheduler::lhs_rhs_element& C = utils::lhs_rhs_element(stcopy, C_idx, C_leaf);
    scheduler::lhs_rhs_element& beta = utils::lhs_rhs_element(stcopy, beta_idx, beta_leaf);






    if (C.numeric_type==scheduler::FLOAT_TYPE)
      enqueue_impl<float>(&scheduler::lhs_rhs_element::matrix_float, stcopy, A, B, C, beta, beta.host_float, programs, kernel_prefix);
    else if (C.numeric_type==scheduler::DOUBLE_TYPE)
      enqueue_impl<double>(&scheduler::lhs_rhs_element::matrix_double, stcopy, A, B, C, beta, beta.host_double, programs, kernel_prefix);
    else
      throw generator_not_supported_exception("GEMM only supported for float/double");

  }

private:
  const char A_trans_;
  const char B_trans_;
};

}

}

#endif
