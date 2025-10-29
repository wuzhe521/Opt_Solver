#ifndef GONS_NEWTON_H
#define GONS_NEWTON_H

#include "function.h"
#include "matrix.h"
#include "vector.h"
#include "utilites.h"
#include "algebra.h"

namespace gons {
using namespace gons::utilites::LOG_MSG;

template <typename Function, typename X> class NewtonMethod {
public:
  struct NewtonParameters {
    double tolerance = 1.0e-6;
    GONS_UINT max_iterations = 3;
    bool verbose = true;
  };
  enum class NewtonStatus { Success, Failure, MaxIterationReached };

public:
  NewtonMethod(const Function &f, const X &x) : f_(f), x_(x) {}

  // compute search direction
  X Search(const X &x, const Function &f) {
    X gradient = f.gradient(x);
   
    gons::Matrix hessian = f.hessian(x);
    
    // gons::Matrix hessian_inv = hessian.Inverse();
    X neg_gradient = gradient * -1;
    // Solve for step: hessian * step = -gradient
    auto step = SolveLinearSystem(hessian, neg_gradient);
    return step;
  }

  void Optimize() {
    for (int i = 0; i < param_.max_iterations; ++i) {
      X gradient = Search(x_, f_);
      x_ += gradient;
      if (gradient.Norm2() < param_.tolerance) {
        LOG("After " << i << " iterations, ");
        LOG("Newton method converged.");
        LOG("Final solution: " << x_);
        LOG("Final gradient: " << gradient);
        LOG("Final function value: " << f_(x_));
        return;
      }
      if (param_.verbose) {
        LOG("Iteration: " << i);
        LOG("x = " << x_);
        LOG("f(x) = " << f_(x_));
        LOG("Gradient: " << gradient);
      }
    }
    LOG("After " << param_.max_iterations << " iterations, ");
    LOG("Newton method did not converge.");
  }

private:
  Function f_;
  X x_;
  NewtonParameters param_;
};

} // namespace gons

#endif // GONS_NEWTON_H