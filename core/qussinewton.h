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
template <typename Function, typename X> class SRMethod {

public:
  enum class SRMethodStatus { Success, Failure, MaxIterationReached };
  enum class SRMethodType { SR1, SR2 };
  struct SRMethodParameters {
    double tolerance = 1.0e-6;
    GONS_UINT max_iterations = 1000;
    GONS_BOOL verbose = false;
    SRMethodType MethodType = SRMethodType::SR1;
  };

  SRMethod(Function &f, X &x) : m_f(f), m_x(x) {}
  ~SRMethod() = default;

  void set_parameters(const SRMethodParameters &params) { m_params = params; }

  X get_x() const { return m_x; }
  double get_function_value() const { return m_f(m_x); }
  template <typename T, GONS_UINT Size>
  Matrix<T, Size, Size> SROneUpdate(const Matrix<T, Size, Size> &B,
                                    const Vector<T, Size> &s,
                                    const Vector<T, Size> &y) {
    Matrix<T, Size, Size> result;
    // Compute the outer product of s and y
    Vector<T, Size> yk_min_BkSk = y - B * s.transpose();
    Vector<T, Size> yk_min_BkSk_T = yk_min_BkSk.transpose();
    T denominator = yk_min_BkSk_T.dot(s);
    // CHECK(FLT_EQUAL(denominator, 0.0), "Denominator is zero");
    Matrix<T, Size, Size> nominator = yk_min_BkSk.outerProduct(yk_min_BkSk);

    result = B + nominator * (1.0 / denominator);
    // Update the Hessian approximation
    return result;
  }
  template <typename T, GONS_UINT Size>
  Matrix<T, Size, Size> SRTwoMethod(const Matrix<T, Size, Size> &H,
                                    const Vector<T, Size> &s,
                                    const Vector<T, Size> &y) {
    Matrix<T, Size, Size> result;
    // Compute the outer product of s and y
    Vector<T, Size> sk_min_Bkyk = s - H * y.transpose();
    Vector<T, Size> sk_min_Bkyk_T = sk_min_Bkyk.transpose();
    T denominator = sk_min_Bkyk_T.dot(y);
    // CHECK(FLT_EQUAL(denominator, 0.0), "Denominator is zero");
    Matrix<T, Size, Size> nominator = sk_min_Bkyk.outerProduct(sk_min_Bkyk);

    result = H + nominator * (1.0 / denominator);
    // Update the Hessian approximation
    return result;
  }

  SRMethodStatus Optimize() {

    GONS_UINT iter = 0;

    Matrix QuassiH =
        m_f.hessian(m_x).Identity(); // Initial Hessian approximation

    while (true) {
      // calc direnction

      X pk = m_params.MethodType == SRMethodType::SR1
                 ? -1 * QuassiH.Inverse() * m_f.gradient(m_x).transpose()
                 : -1 * QuassiH * m_f.gradient(m_x).transpose();
      // line search to find step length Todo: change to wolfe search
      // update x
      X x_new = m_x + pk;
      // update B
      X sk = x_new - m_x;

      X yk = m_f.gradient(x_new) - m_f.gradient(m_x);

      QuassiH = m_params.MethodType == SRMethodType::SR1
                    ? SROneUpdate(QuassiH, sk, yk)
                    : SRTwoMethod(QuassiH, sk, yk);

      m_x = x_new;
      ++iter;

      if (m_f.gradient(m_x).Norm2() < m_params.tolerance) {

        LOG("SRMethod: Optimization using "
            << (m_params.MethodType == SRMethodType::SR1 ? "SR1" : "SR2")
            << " finished.");
        LOG("SRMethod: After Iterations: " << iter);
        LOG("SRMethod: Final X : " << m_x);
        LOG("SRMethod: Function Function value: " << m_f(m_x) << '\n');

        return SRMethodStatus::Success;
      }
      if (iter++ >= m_params.max_iterations) {
        LOG_ERROR("SRMethod: Maximum number of iterations reached.");

        return SRMethodStatus::MaxIterationReached;
      }
      if (m_params.verbose) {
        LOG("Iteration: " << iter);
        LOG("x: " << m_x);
      }
    }
  }

private:
  Function &m_f;
  X &m_x;
  SRMethodParameters m_params;
  // Line search method Todo: implement different line search methods
  // linearsearch::WolfeSearch<Function, X> m_wolfe_search;
};
} // namespace qussinewton
} // namespace gons

#endif // GONS_QUASSINEWTON_H