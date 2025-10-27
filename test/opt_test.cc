#include "../core/linearsearch.h"


using namespace gons;

using namespace gons::utilites::LOG_MSG;


using X = gons::Vector<GONS_FLOAT, 2 >;

// 定义测试函数类，用于计算目标函数值和梯度
// 目标函数为x1^2 + x2^2

class TestFunction {
public:
    inline GONS_FLOAT operator()(const X &x) const
    {
        return x(0) * x(0) + x(1) * x(1);
    }
    inline X gradient(const X &x)
    {
        return {2 * x(0), 2 * x(1)};
    }
};
int main()
{ 
    TestFunction f;
    X x = {10, 10};
    gons::ArmijoSearch<TestFunction, X> armijo_search(f, x);

    gons::ArmijoSearch<TestFunction, X>::ArmijoParameters params;

    params.alpha = 0.01;
    params.beta = 0.2;
    params.gamma = 0.333;
    params.enable_max_iter = true;
    params.epsilon = 1e-6;
    params.print_info = false;
    armijo_search.set_params(params);
    armijo_search.Optimize();
    LOG("Final x: " << armijo_search.get_x());
    LOG("Final f(x): " << armijo_search.get_function_value());
    return 0;
}
