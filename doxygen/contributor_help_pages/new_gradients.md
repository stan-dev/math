## Adding A New Function With Known Gradients {#new_grad}

-----------------------------------------------------

### NOTE:

These docs are a bit outdated. For the most up to date guide on contributing to Stan Math see the [Getting Started Guide](@ref getting_started), [Common Pitfalls](@ref common_pitfalls), and [Adding New Distributions Guide](@ref new_distribution).


------------------------------------------------------

If you have a function f and you know the gradients for it, it is straightforward to add the function to Stan's math library.


#### Templated function

If your function is templated separately on all of its arguments and bottoms out only in functions built into the C++ standard library or into Stan, there's nothing to do other than add it to the Stan namespace.  For example, if we want to code in a simple z-score calculation, we can do it very easily as follows, because the `mean` and `sd` function are built into Stan's math library:


```cpp
#include <stan/math/rev/core.hpp>

namespace stan {
  namespace math {
    vector<T> z(const vector<T>& x) {
      T mean = mean(x);
      T sd = sd(x);
      vector<T> result;
      for (size_t i = 0; i < x.size(); ++i)
        result[i] = (x[i] - mean) / sd;
      return result;
    }
  }
}
```

This will even work for higher-order autodiff if you include the top-level header `<stan/math.hpp>`.

#### Simple univariate example with known derivatives

Suppose have a code to calculate a univariate function and its derivative:

```cpp
namespace bar {
  double foo(double x);
  double d_foo(double x);
}
```

We can implement foo for Stan as follows, making sure to put it in the `stan::math` namespace so that argument-dependnet lookup (ADL) can work.  We'll assume the two relevant functions are defined in `my_foo.hpp`

```cpp
#include "my_foo.hpp"
#include <stan/math/rev/core.hpp>

namespace stan {
  namespace math {
    double foo(double x) {
      return bar::foo(x);
    }

    var foo(const var& x) {
      double a = x.val();
      double fa = bar::foo(x_d);
      double dfa_da = bar::d_foo(a);
      return precomp_v_vari(fa, x.vi_, dfa_da);
    }
  }
}
```

There are similar functions `precomp_vv_vari` for two-argument functions and so on.



#### Functions of more than one argument

For most efficiency, each argument should independently allow `var` or `double` arguments (`int` arguments don't need gradients---those can just get carried along).  So if you're writing a two-argument function of scalars, you want to do this:

```cpp
namespace bar {
  double foo(double x, double y);
  double dfoo_dx(double x, double y);
  double dfoo_dy(double x, double y);
}
```

This can be coded as follows

```cpp
#include "my_foo.hpp"
#include <stan/math/rev/core.hpp>

namespace stan {
  namespace math {
    double foo(double x, double y) {
      return bar::foo(x, y);
    }

    var foo(const var& x, double y) {
      double a = x.val();
      double fay = bar::foo(a, y);
      double dfay_da = bar::foo_dx(a, y);
      return return precomp_v_vari(fay, x.vi_, dfay_da);
    }

    var foo(double x, const var& y) {
      double b = y.val();
      double fxb = bar::foo(x, b);
      double dfxb_db = bar::foo_dy(x, b);
      return return precomp_v_vari(fxb, y.vi_, dfxb_db);
    }

    var foo(const var& x, const var& y) {
      double a = x.val();
      double b = y.val();
      double fab = bar::foo(a, b);
      double dfab_da = bar::foo_dx(a, b);
      double dfab_db = bar::foo_dy(a, b);
      return precomp_vv_vari(fab, x.vi_, y.vi_, dfab_da, dfab_db);
    }
  }
}
```

Same thing works for three scalars.


#### Vector function with one output

Now suppose the signatures you have give you a function and a gradient using standard vectors (you can also use Eigen vectors or matrices in a similar way):

```cpp
namespace bar {
  double foo(const vector<double>& x);
  vector<double> grad_foo(const vector<double>& x);
}
```

Of course, it's often more efficient to compute these at the same time---if you can do that, replace the lines for `fa` and `grad_fa` witha more efficient call.  Otherwise, the code follows the previous version:

```cpp
#include "my_foo.hpp"
#include <stan/math/rev/core.hpp>

namespace stan {
  namespace math {
    vector<double> foo(const vector<double>& x) {
      return bar::foo(x);
    }

    vector<var> foo(const vector<var>& x) {
      vector<double> a = value_of(x);
      double fa = bar::foo(a);
      vector<double> grad_fa = bar::grad_foo(a);
      return precomputed_gradients(fa, x, grad_fa);
    }
  }
}
```

The `precomputed_gradients` class can be used to deal with any form of inputs and outputs.  All you need to do is pack all the `var` arguments into one vector and their matching gradients into another.

#### Functions with more than one output

If you have a function with more than one output, you'll need the full Jacobian matrix.  You create each output `var` the same way you create a univariate output, using `precomputed_gradients`.  Then you just put them all together.  This time, I'll use Eigen matrices assuming that's how you get your Jacobian:

```cpp
namespace bar {
  Eigen::VectorXd foo(const Eigen::VectorXd& x);
  Eigen::MatrixXd foo_jacobian(const Eigen::VectorXd& x);
}
```

You can code this up the same way as before (with the same includes):

```cpp
namespace stan {
  namespace math {
    Eigen::VectorXd foo(const Eigen::VectorXd& x) {
      return bar::foo(x);
    }
    Eigen::Matrix<var, -1, 1> foo(const Eigen::Matrix<var,-1,1>& x,
                                  std::ostream* pstream__) {
      // extract parameter values from x
      const var * x_data = x.data();
      int x_size = x.size();
      vector<var> x_std_var(x_data, x_data + x_size );
      VectorXd a = value_of( x );

      // evaluate f
      VectorXd f_val = bar::foo(a);
      // f_val_jacobian[i][j] is the partial derivative of f_i w.r.t. parameter j
      vector<vector<double>> f_val_jacobian = bar::foo_jacobian(a);

      int f_size = f_val_jacobian.size();
      Eigen::Matrix<var, -1, 1> f_var_jacobian(f_size);
      for (int i=0; i < f_size; i++) {
        f_var_jacobian(i) = precomputed_gradients(f_val(i), x_std_var, f_val_jacobian[i]);
      }
      return f_var_jacobian;
    }
  }
}
```

We could obviously use a form of `precomputed_gradients` that takes in an entire Jacobian to remove the error-prone and non-memory-local loops (by default, matrices are stored column-major in Eigen).
