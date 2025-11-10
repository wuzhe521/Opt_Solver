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
// quadratic const function
// linear constraint
// DataType :  float or double
// Var_Size :  number of variables
// Con_Size :  number of constraints
// description:
//             f(x) = 0.5 * x' * P * x + q' * x
//             s.t.
//             l <= A * x <= u


template <typename DataType, GONS_UINT Var_Size, GONS_UINT Con_Size>
class QuadPenaltyFunction {
public:
class function {
  private:
    GONS_FLOAT rho = 10.0;
    Matrix<DataType, Var_Size, Var_Size> P_;
    Vector<DataType, Var_Size> q_;
    Vector<DataType, Con_Size> l_;
    Vector<DataType, Con_Size> u_;
    Matrix<DataType, Con_Size, Var_Size> A_;

  public:
    function(Matrix<DataType, Var_Size, Var_Size> &P,
             Vector<DataType, Var_Size> &q, Vector<DataType, Con_Size> &l,
             Vector<DataType, Con_Size> &u,
             Matrix<DataType, Con_Size, Var_Size> &A)
        : P_(P), q_(q), l_(l), u_(u), A_(A) {} // Fixed order here
    DataType operator()(Vector<DataType, Var_Size> &x) const {
      return CostFunction(x) + 0.5 * rho * penaltyFunction(x);
    }
    Vector<DataType, Var_Size> gradient(Vector<DataType, Var_Size> &x) {
      return CostFunction_grad(x) + rho * grad_penaltyFunction(x);
    }
    double CostFunction(Vector<DataType, Var_Size> &x) const {
      return x.dot(P_ * x) + q_.dot(x);
    }
    Vector<DataType, Var_Size>
    CostFunction_grad(Vector<DataType, Var_Size> &x) const {
      return P_ * x + q_;
    }
    double penaltyFunction(Vector<DataType, Var_Size> &x) const {
      DataType sum = 0;
      Vector<DataType, Con_Size> Ax_ = A_ * x.transpose();
      for (size_t i = 0; i < Con_Size; i++) {
        DataType temp = l_(i) - Ax_(i); // g(x)_lower : l - Ax
        if (temp > 0)
          sum += temp * temp;
        temp = Ax_(i) - u_(i); // g(x)_upper : Ax - u
        if (temp > 0)
          sum += temp * temp; // max(0, g(x))^2
      }
      return sum;
    }
    Vector<DataType, Var_Size> // To be Done  ... !!!
    grad_penaltyFunction(Vector<DataType, Var_Size> &x) {
      Vector<DataType, Var_Size> grad;
      Vector<DataType, Con_Size> Ax = A_ * x.transpose();
      for (size_t i = 0; i < Con_Size; i++) {
        DataType temp = l_(i) - Ax(i);
        if (temp > 0)
          return -1* (A_ * x.transpose() + l_) * A_;
        temp = Ax(i) - u_(i);
        if (temp < 0)
          return  (A_ * x.transpose() - u_) * A_;
      }
      return grad;
    }
  };
public:
  struct quad_penalty_parmeters {
    DataType rho = 0.1;
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
                      Vector<DataType, Con_Size> &l,
                      Vector<DataType, Con_Size> &u,
                      Matrix<DataType, Con_Size, Var_Size> &A)
      : P_(P), q_(q), l_(l), u_(u), A_(A) {}

  PenaltyOptStatus Optimize() {

    GONS_UINT iter = 0;
    do {
      // Create a function object
      function f(P_, q_, l_, u_, A_, parameters_.rho);
      BarzilaiBorwein<function, Vector<GONS_FLOAT, Var_Size>>
          barzilai_borwein_search(f, x_);
      barzilai_borwein_search.Optimize();
      x_ = barzilai_borwein_search.get_x();
      // 判断是否结束循环
      if (f.grad_penaltyFunction(x_).Norm2() < parameters_.tol) {
        LOG("满足结束条件")
        LOG("迭代次数为：" << iter)
        LOG("当前解为：" << x_.transpose())

        return PenaltyOptStatus::Success;
      }
      // 不满足结束条件， 继续增大rho
      parameters_.rho *= parameters_.alpha;

    } while (iter < parameters_.max_iter);
    LOG_ERROR("迭代次数超过最大限制")
    return PenaltyOptStatus::MaxIterationReached;
  }

private:
  Matrix<DataType, Var_Size, Var_Size> P_;
  Vector<DataType, Var_Size> q_;
  Vector<DataType, Con_Size> l_;
  Vector<DataType, Con_Size> u_;
  Matrix<DataType, Con_Size, Var_Size> A_;
  // optimal variables
  Vector<DataType, Var_Size> x_;

  quad_penalty_parmeters parameters_;
  // unconstrained search method
};

} // namespace penalty_function
} // namespace gons

#endif