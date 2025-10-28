#include "../core/gradientsearch.h"

// Test Function
// f(x) = x^2 + y^2
using X = gons::Vector<double, 2>;
class TestFunction {
public:
    TestFunction() = default;
    ~TestFunction() = default;
    
    double operator()(const X& x) const {
        return x(0) * x(0) + x(1) * x(1);
    }
    X gradient(const X& x) const {
        
        return {2 * x(0), 2 * x(1)};
    }
};

int main( int argc, char* argv[] )
{ 
    UNUSED(argc);
    UNUSED(argv);

    TestFunction f;
    X x = {1.0, 1.0};
    gons::GradientDescentSearch<TestFunction, X> gd(f, x);

    gd.Optimize();

    return 0;
}
