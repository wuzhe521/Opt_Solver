#ifndef GONS_CORE_FUNCTION_H_
#define GONS_CORE_FUNCTION_H_

#include "config.h"
#include "matrix.h"
#include "utilites.h"
#include "vector.h"

namespace gons {
    template <typename T>
    class Rosenbrock {
    public:
        Rosenbrock();
        ~Rosenbrock();
        
        // Function evaluation
        T operator()const Vector<T, 2> & x) const{
            return (1 - x(0)) * (1 - x(0)) + 100 * (x(1) - x(0) * x(0)) * (x(1) - x(0) * x(0));
        }
        
        // Gradient evaluation
        Vector<T, 2> gradient(const Vector<T, 2>& x) const {
            Vector<T, 2> grad;
            grad(0) = -2 * (1 - x(0)) - 200 * x(0) * (x(1) - x(0) * x(0));
            grad(1) = 200 * (x(1) - x(0) * x(0));
            return grad;
        }
    };

} // namespace gons

#endif // GONS_CORE_FUNCTION_H_