#include "../core/quad_penalty.h"

using namespace gons;
using namespace ::utilites::LOG_MSG;

using X = ::Vector<double, 2>;

int main() {

  Matrix<double, 2, 2> P = {{4.0, 1.0}, {1.0, 2.0}};
  Vector<double, 2> q = {1.0, 1.0};
  Matrix<double, 1, 2> A = {{1.0, 1.0}};
  Vector<double, 1> b = {1.0};
  penalty_function::QuadPenaltyFunction<double, 2, 1> PenaltyFunc(P, q, A, b);
  PenaltyFunc.Optimize();
  X rlt = PenaltyFunc.get_result();
  std::cout << "result: " << rlt << std::endl;
  
  return 0;
}