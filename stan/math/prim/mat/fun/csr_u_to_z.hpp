#ifndef STAN_MATH_PRIM_MAT_FUN_CSR_U_TO_Z
#define STAN_MATH_PRIM_MAT_FUN_CSR_U_TO_Z

#include <stdexcept>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Convert the u data structure to the z data structure for CSR
     * storage of matrices.
     *
     * @param[in] u List of row starts.
     * @return List of non-zero entries in rows.
     * @throws If list of indexes is empty.
     */
    vector<int> csr_u_to_z(const vector<int>& u) {
      if (u.size() == 0)
        throw std::domain_error("csr_u_to_z: u.size() == 0,"
                                " but must be positive");
      vector<int> z(u.size() - 1);
      for (size_t i = 0; i < z.size(); ++i)
        z[i] = u[i + 1] - u[i];
      return z;
    }

    void print_vec(const vector<int>& x, const std::string& name) {
      std::cout << std::endl
                << name << ".size()=" << x.size()
                << std::endl;
      for (size_t i = 0; i < x.size(); ++i)
        std::cout << name << "[" << i << "]=" << x[i] << std::endl;
      std::cout << std::endl;
    }

  }
}
#endif
