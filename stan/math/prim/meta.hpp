#ifndef STAN_MATH_PRIM_META_HPP
#define STAN_MATH_PRIM_META_HPP

#include <iostream>

#include <stan/math/prim/meta/get.hpp>
#include <stan/math/prim/meta/index_type.hpp>

#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/length.hpp>

#include <stan/math/prim/meta/append_return_type.hpp>
#include <stan/math/prim/meta/as_array_or_scalar.hpp>
#include <stan/math/prim/meta/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/meta/as_scalar.hpp>
#include <stan/math/prim/mat/meta/broadcast_array.hpp>
//#include <stan/math/prim/meta/is_constant.hpp>  // come back to it
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <stan/math/prim/meta/length_mvt.hpp>
//#include <stan/math/prim/meta/operands_and_partials.hpp>  // come back to it
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/seq_view.hpp>
//#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/vector_seq_view.hpp>

#include <stan/math/prim/arr/meta/contains_std_vector.hpp>
//#include <stan/math/prim/arr/meta/is_constant.hpp>
//#include <stan/math/prim/arr/meta/scalar_type.hpp>
//#include <stan/math/prim/arr/meta/value_type.hpp>
#include <stan/math/prim/meta/VectorBuilderHelper.hpp>

#include <stan/math/prim/meta/ad_promotable.hpp>
//#include <stan/math/prim/meta/child_type.hpp>
#include <stan/math/prim/meta/contains_fvar.hpp>
//#include <stan/math/prim/meta/contains_std_vector.hpp>
#include <stan/math/prim/meta/contains_vector.hpp>
#include <stan/math/prim/meta/error_index.hpp>
//#include <stan/math/prim/scal/meta/include_summand.hpp>  // causing bug
//#include <stan/math/prim/scal/meta/is_constant.hpp>
#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>
//#include <stan/math/prim/scal/meta/max_size_mvt.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/partials_type.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
//#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/scalar_type_pre.hpp>
#include <stan/math/prim/scal/meta/size_of.hpp>
#include <stan/math/prim/scal/meta/value_type.hpp>
#include <stan/math/prim/scal/meta/StdVectorBuilder.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>

#endif
