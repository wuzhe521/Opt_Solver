#include "../core/quad_penalty.h"

using namespace gons;
using namespace ::utilites::LOG_MSG;

using X = ::Vector<double, 2>;

int main() {

  Matrix<double, 2, 2> P = {{4.0, 1.0}, {1.0, 2.0}};
  Vector<double, 2> q = {1.0, 1.0};
  Matrix<double, 3, 2> A = {{1.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}};

  Vector<double, 3> l = {1.0, 0.0, 0.0};
  Vector<double, 3> u = {1.0, 0.7, 0.7};

  Vector<double, 2> x = {0.0, 0.0};

  /* penalty_function::QuadPenaltyFunction<double, 2, 3> PenaltyFunc(P, q, l, u,
                                                                  A);*/
  penalty_function::QuadPenaltyFunction<double, 2, 3>::function f(P, q, l, u,
                                                                  A);

  f(x);
  f.gradient(x);
  // PenaltyFunc.Optimize();
  return 0;
}