#ifndef GONS_QUASSINEWTON_H
#define GONS_QUASSINEWTON_H

#include "function.h"
#include "matrix.h"
#include "utilites.h"
#include "vector.h"

#include "linearsearch.h"

namespace gons {
namespace qussinewton {

using namespace gons::utilites::LOG_MSG;
template <typename Function, typename X> class QuasiNewtonMethod {

public:
  struct QuasiNewtonParameters {
    double tolerance = 1.0e-6;
    GONS_UINT max_iterations = 1000;
    bool verbose = true;
  };
  enum class QuasiNewtonStatus { Success, Failure, MaxIterationReached };
  QuasiNewtonMethod(Function &f, X &x) : m_f(f), m_x(x), m_wolfe_search(f, x) {}
  ~QuasiNewtonMethod() = default;
  template <typename T, GONS_UINT Size>
  Matrix<T, Size, Size> SROneUpdate(const Matrix<T, Size, Size> &H,
                                    const Vector<T, Size> &s,
                                    const Vector<T, Size> &y) {
    Matrix<T, Size, Size> result;
    // Compute the outer product of s and y
    Vector<T, Size> yk_min_BkSk = y - H * s;
    Vector<T, Size> yk_min_BkSk_T = yk_min_BkSk.transpose();
    // Update the Hessian approximation
    return result;
  }

  QuasiNewtonStatus Optimize() {

    GONS_UINT iter = 0;

    Matrix Bk = m_f.hessian(m_x).Identity(); // Initial Hessian approximation
    while (true) {
      // calc direnction
      X pk = -Bk * m_f.gradient(m_x);
      // line search
      double alpha = m_wolfe_search.Search(m_f, m_x, pk);
      // update x
      X x_new = m_x + alpha * pk;
      // update B
      X sk = x_new - m_x;
      X yk = m_f.gradient(x_new) - m_f.gradient(m_x);
      Bk = SROneUpdate(Bk, sk, yk);

      m_x = x_new;
      ++iter;

      if (m_f.gradient(m_x).Norm2() < m_params.tolerance) {

        LOG("QuasiNewton: Optimization finished.");
        LOG("QuasiNewton: Function value: " << m_f.value(m_x));
        LOG("QuasiNewton: Iterations: " << iter);

        return QuasiNewtonStatus::Success;
      }
      if (iter++ >= m_params.max_iterations) {
        LOG_ERROR("QuasiNewton: Maximum number of iterations reached.");
        return QuasiNewtonStatus::MaxIterationReached;
      }
      if (m_params.verbose) {
        LOG("Iteration: " << iter);
        LOG("Function value: " << m_f.value(m_x));
        LOG("Gradient norm: " << m_f.gradient(m_x).Norm2());
        LOG("x: " << m_x);
      }
    }
  }

private:
  Function &m_f;
  X &m_x;
  QuasiNewtonParameters m_params;
  // Line search method
  linearsearch::WolfeSearch<Function, X> m_wolfe_search;
  //
};
} // namespace qussinewton
} // namespace gons

#endif // GONS_QUASSINEWTON_H