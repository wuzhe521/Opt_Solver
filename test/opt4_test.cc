#include "../core/qussinewton.h"

using namespace gons::utilites::LOG_MSG;

using X = gons::Vector<double, 2>;
// 定义测试函数类，用于计算目标函数值和梯度
class TestFunction {
public:
  inline double operator()(const X &x) const {
    return x(0) * x(0) + x(1) * x(1);
  }
  inline X gradient(const X &x) { return {2 * x(0), 2 * x(1)}; }
  inline gons::Matrix<double, 2, 2> hessian(const X &x) {
    return gons::Matrix<double, 2, 2>({{2, 0}, {0, 2}});
  }
};
int main() {
  TestFunction f;
  X x0 = {10, 10};
  gons::qussinewton::QuasiNewtonMethod<TestFunction, X> qussi_newton(f, x0);
  qussi_newton.Optimize();
  LOG(gons::utilites::GetVersion());
  LOG("QuasiNewtonMethod");
  return 0;
}