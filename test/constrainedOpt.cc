#include "../core/quad_penalty.h"

using namespace gons;
using namespace ::utilites::LOG_MSG;

using X = ::Vector<double, 2>;

int main() {

  Matrix<double, 2, 2> P = {{4.0, 1.0}, {1.0, 2.0}};
  Vector<double, 2> q = {1.0, 1.0};
  Matrix<double, 1, 2> A = {{1.0, 1.0}};
  Vector<double, 1> b = {1.0};
  penalty_function::QuadPenaltyFunction_ECQP<double, 2, 1> PenaltyFunc(P, q, A,
                                                                       b);
  PenaltyFunc.Optimize();
  X rlt = PenaltyFunc.get_result();
  std::cout << "result: " << rlt << std::endl;

  Matrix<double, 2, 2> P2 = {{4.0, 1.0}, {1.0, 2.0}};
  Vector<double, 2> q2 = {1.0, 1.0};
  Matrix<double, 3, 2> A2 = {{1.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}};
  Vector<double, 3> l2 = {1.0, 0.0, 0.0};
  Vector<double, 3> u2 = {1.0, 0.7, 0.7};

  penalty_function::QuadPenaltyFunction_IECQP<double, 2, 3> PenaltyFunc2(P2, q2, A2, l2, u2);
  PenaltyFunc2.Optimize();
  std::cout << "result: " << PenaltyFunc2.get_result() << std::endl;
  return 0;
}