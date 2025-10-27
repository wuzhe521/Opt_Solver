#include "../core/linearsearch.h"

using namespace gons;

using namespace gons::utilites::LOG_MSG;

using X = gons::Vector<GONS_FLOAT, 2>;

// 定义测试函数类，用于计算目标函数值和梯度
// 目标函数为x1^2 + x2^2

class TestFunction {
public:
  inline GONS_FLOAT operator()(const X &x) const {
    return x(0) * x(0) + x(1) * x(1);
  }
  inline X gradient(const X &x) { return {2 * x(0), 2 * x(1)}; }
};
int main() {
  TestFunction f;
  X x = {10, 10};
  // 使用Armijo搜索优化器
  gons::ArmijoSearch<TestFunction, X> armijo_search(f, x);
  armijo_search.Optimize();
  LOG("Final x: " << armijo_search.get_x());
  LOG("Final f(x): " << armijo_search.get_function_value());
  // 使用Goldstein搜索优化器
  gons::GoldsteinSearch<TestFunction, X> goldstein_search(f, x);
  gons::GoldsteinSearch<TestFunction, X>::GoldsteinParameters goldstein_params;
  goldstein_search.Optimize();
  LOG("Final x: " << goldstein_search.get_x());
  LOG("Final f(x): " << goldstein_search.get_function_value());

  // 使用Wolfe搜索优化器
  gons::WolfeSearch<TestFunction, X> wolfe_search(f, x);
  gons::WolfeSearch<TestFunction, X>::WolfeParameters wolfe_params;
  wolfe_search.Optimize();
  LOG("Final x: " << wolfe_search.get_x());
  LOG("Final f(x): " << wolfe_search.get_function_value());
  return 0;
}
