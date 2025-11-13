#include "../core/newton.h"

// Test Function
// f(x) = x^2 + y^2
using X = gons::Vector<double, 2>;
using M = gons::Matrix<double, 2, 2>;
class TestFunction {
public:
  TestFunction() = default;
  ~TestFunction() = default;

  double operator()(const X &x) const { return x(0) * x(0) + 10 * x(1) * x(1); }
  X gradient(const X &x) const { return {2 * x(0), 20 * x(1)}; }
  M hessian(const X &x) const { return {{2, 0}, {0, 20}}; }
};

int main(int argc, char *argv[]) {
  UNUSED(argc);
  UNUSED(argv);

  TestFunction f;
  X x = {10, 10};
  gons::newton::NewtonMethod<TestFunction, X> nm(f, x);
  nm.Optimize();

  return 0;
}
