#ifndef GONS_QUAD_PENALTY_H
#define GONS_QUAD_PENALTY_H

#include "algebra.h"
#include "config.h"
#include "constant.h"
#include "function.h"
#include "gradientsearch.h"
#include "linearsearch.h"
#include "matrix.h"
#include "newton.h"
#include "qussinewton.h"
#include "utilites.h"
#include "vector.h"

namespace gons {
namespace penalty_function {
using namespace gons::utilites::LOG_MSG;
using namespace gons::gradientsearch;

template <typename DataType, GONS_UINT Var_Size, GONS_UINT Con_Size>
class function {
private:
  GONS_FLOAT rho_ = 10.0;
  Matrix<DataType, Var_Size, Var_Size> P_;
  Vector<DataType, Var_Size> q_;
  Matrix<DataType, Con_Size, Var_Size> A_;
  Vector<DataType, Con_Size> b_;

public:
  function(Matrix<DataType, Var_Size, Var_Size> &P,
           Vector<DataType, Var_Size> &q, 
           Matrix<DataType, Con_Size, Var_Size> &A,
           Vector<DataType, Con_Size> &b)
      : P_(P), q_(q), A_(A), b_(b) {} // Fixed order here
  DataType operator()(Vector<DataType, Var_Size> &x) const {
    return 0.5 * CostFunction(x) + 0.5 * rho_ * penaltyFunction(x);
  }
  Vector<DataType, Var_Size> gradient(Vector<DataType, Var_Size> &x) {
    return CostFunction_grad(x) + rho_ * grad_penaltyFunction(x);
  }
  void update_rho(GONS_FLOAT rho) { rho_ = rho; }

public:
  double CostFunction(Vector<DataType, Var_Size> &x) const {
    return x.dot(P_ *x.transpose()) + q_.dot(x);
  }
  Vector<DataType, Var_Size>
  CostFunction_grad(Vector<DataType, Var_Size> &x) const {
    return x * P_ + q_;
  }
  double penaltyFunction(Vector<DataType, Var_Size> &x) const {
    DataType sum = 0;
    Vector<DataType, Con_Size> Conflict = A_ * x.transpose() - b_;
    for (size_t i = 0; i < Con_Size; i++) {
      sum += Conflict(i) * Conflict(i);
    }
    return sum;
  }
  Vector<DataType, Var_Size> // To be Done  ... !!!
  grad_penaltyFunction(Vector<DataType, Var_Size> &x) {
    Vector<DataType, Con_Size> Ax = A_ * x.transpose();
    Vector<DataType, Con_Size> Conflict = Ax - b_;
    Vector<DataType, Var_Size> grad = Conflict * A_;
    return grad;
  }
};
// quadratic const function
// linear constraint
// DataType :  float or double
// Var_Size :  number of variables
// Con_Size :  number of constraints
// description:
//             f(x) = 0.5 * x' * P * x + q' * x
//             s.t.
//             Ax = b

template <typename DataType, GONS_UINT Var_Size, GONS_UINT Con_Size>
class QuadPenaltyFunction {
  using X = Vector<DataType, Var_Size>;
  using F = function<DataType, Var_Size, Con_Size>;

public:
  struct quad_penalty_parmeters {
    DataType rho = 1;
    DataType alpha = 2; // penalty growth rate
    GONS_UINT max_iter = 1000;
    GONS_BOOL verbose = false;
    GONS_BOOL war_start = true;
    DataType tol = 1e-5;
  };
  enum class PenaltyOptStatus { Failure, MaxIterationReached, Success };

public:
  // for quadratic function:
  // f(x) = 0.5 * x' * P * x + q' * x
  // for linear inequality constraint:
  // l <= g(x) = A * x + b <= u

  QuadPenaltyFunction(Matrix<DataType, Var_Size, Var_Size> &P,
                      Vector<DataType, Var_Size> &q,
                      Matrix<DataType, Con_Size, Var_Size> &A,
                      Vector<DataType, Con_Size> &b)
      : P_(P), q_(q), A_(A), b_(b), f_(P, q, A, b) {}
  Vector<DataType, Var_Size> get_result() const {
    return x_;
  }
  PenaltyOptStatus Optimize() {

    GONS_UINT iter = 0;
    do {
      qussinewton::BFGSMethod<F, X> BFGS(f_, x_);
      BFGS.Optimize();
      x_ = BFGS.get_x();
      // 判断是否结束循环
      if (f_.grad_penaltyFunction(x_).Norm2() < parameters_.tol) {
        LOG("满足结束条件")
        LOG("迭代次数为：" << iter)
        LOG("当前解为：" << x_.transpose())
        return PenaltyOptStatus::Success;
      }
      // 不满足结束条件， 继续增大rho
      parameters_.rho *= parameters_.alpha;
      f_.update_rho(parameters_.rho);

    } while (iter < parameters_.max_iter);
    LOG_ERROR("迭代次数超过最大限制")
    return PenaltyOptStatus::MaxIterationReached;
  }

private:
  Matrix<DataType, Var_Size, Var_Size> P_;
  Vector<DataType, Var_Size> q_;
  Matrix<DataType, Con_Size, Var_Size> A_;
  Vector<DataType, Con_Size> b_;
  // optimal variables
  Vector<DataType, Var_Size> x_;

  quad_penalty_parmeters parameters_;
  // Create a function object
  function<DataType, Var_Size, Con_Size> f_;
  // unconstrained search method
};

} // namespace penalty_function
} // namespace gons

#endif