# Gradient Based Search

梯度类算法，其本质是仅仅使用函数的一阶导数信息选取下降方向 $d^k$. 其中最基本的梯度类算法为梯度下降法，直接选取负梯度作为下降方向．
梯度下降法的方向选取非常直观，实际应用范围非常广，因此它在优化算法中的地位可相当于高斯消元法在线性方程组算法中的地位．


对于函数 $f(x, y) = x^2 + 10*y^2$， 其一阶导数为 $\frac{\partial f}{\partial x} = 2x$, $\frac{\partial f}{\partial y} = 20y$, 因此其负梯度为 $-2x, -20y$, 因此其下降方向为 $-2x, -20y$, 其下降方向的选择是基于一阶导数的。
梯度下降法的演示代码如下：
```python
import numpy as np
import matplotlib.pyplot as plt


plt.figure()
x = np.linspace(-12, 12, 100)
y = np.linspace(-12, 12, 100)
X, Y = np.meshgrid(x, y)
Z = X**2 + Y**2
plt.contour(X, Y, Z, 20)

def f(x):
    return x[0]**2 + 10 * x[1]**2
    
def df(x):
    return np.array([2*x[0], 20*x[1]])
x0 = np.array([-10.0, -10.0])
x = []
x.append(x0)
max_iter = 100
step_size = 0.085
for i in range(max_iter):
    x_new = x[-1] - step_size * df(x[-1])
    x.append(x_new)
    plt.plot([x[-2][0], x_new[0]], [x[-2][1], x_new[1]], 'r-')
    plt.plot(x_new[0], x_new[1], '*', )
    plt.text(x_new[0], x_new[1], f'{i}')

    print(f'Iteration {i}: x = {x_new}, f(x) = {f(x_new)}')
    if np.linalg.norm(df(x[-1])) < 1e-6:
        break
    if i == max_iter - 1:
        print('Maximum number of iterations reached')
    
plt.plot(x[-1][0], x[-1][1], '*')
plt.show()
```
在等高线视图中，依次画出梯度下降路径，如下：

![](./pics/gradient_py.png)


在固定算法的基础上，将固定步长改为自适应步长，步长的求解方式采用线搜索方法的Armijo方法，就实现了梯度类算法，
总体逻辑：
1. 输入初始点 $x_0$，初始步长 $\alpha_0$，最大迭代次数 $max\_iter$，步长缩减因子 $\beta$，步长增加因子 $\gamma$，精度 $\epsilon$
2. 迭代 $i = 0, 1, 2, ..., max\_iter$ 
3. 计算当前点 $x_i$ 的一阶导数 $g_i$, 函数值 $f_i$,
4. 根据当前点 $x_i$ 的函数值 $f_i$ 和导数 $g_i$，用线搜算法计算当前点 $x_i$ 的下降步长 $\beta_i$
5. 更新 $x_{i+1} = x_i - \beta_i g_i$
6. 判断当前点 $x_{i+1}$ 的函数值 $f_{i+1}$ 和上一次点 $x_i$ 的函数值 $f_i$差是否小于精度 $\epsilon$，小于则迭代结束，否则继续迭代
梯度类算法的实现代码如下：
```cpp
template <typename Function, typename X> class GradientDescentSearch {
  struct gradient_descent_parameters {
    double tolerance = 1.0e-6;
    GONS_UINT max_iterations = 1000;
    bool verbose = false;
  };
  enum class gradient_descent_state { SUCCESS, FAILURE, MAX_ITERATION_REACHED };

public:
  // Corrected constructor with matching initialization order
  GradientDescentSearch(const Function &f, const X &x) 
      : f_(f), x_(x), linear_search_(f, x){}

  ~GradientDescentSearch() = default;
  void set_params(const gradient_descent_parameters &parameters) {
    parameters_ = parameters;
  }
  double SearchStep(const X &x, const Function &f) {
    return linear_search_.Search(x, f);
  }
  gradient_descent_state Optimize() {
    GONS_UINT iter = 0;
    GONS_UINT max_iter = parameters_.max_iterations;
    GONS_FLOAT tolerance = parameters_.tolerance;

    while (iter < max_iter) {
      // 计算梯度
      X gradient = f_.gradient(x_);

      // 计算步长
      double step = SearchStep(x_, f_);

      // 更新x
      X x_new = x_ - step * gradient;
      if (parameters_.verbose) // 打印信息
      {
        LOG("Iteration: " << iter);
        LOG("x = " << x_);
        LOG("f(x) = " << f_(x_));
        LOG("Step = " << step);
      }
      if (std::abs(f_(x_new) - f_(x_)) < tolerance) {
        LOG_WARNING("After " << iter << " iterations");
        LOG_WARNING("Gradient descent converged.");
        LOG_WARNING("x = " << x_);
        LOG_WARNING("f(x) = " << f_(x_));
        return gradient_descent_state::SUCCESS;
      }

      x_ = x_new;
      iter++;
    }
    LOG_WARNING("After " << iter << " iterations");
    LOG_WARNING("Gradient descent did not converge.");
    LOG_WARNING("x = " << x_);
    LOG_WARNING("f(x) = " << f_(x_));
    return gradient_descent_state::MAX_ITERATION_REACHED;
  }

private:
  Function f_;
  X x_;
  ArmijoSearch<Function, X> linear_search_;
  gradient_descent_parameters parameters_;

};
```
测试代码如下：
```cpp
TestFunction f;
X x = {1.0, 1.0};
gons::GradientDescentSearch<TestFunction, X> gd(f, x);
gd.Optimize();
```
运行结果如下：
```bash
After 14 iterations
Gradient descent converged.
x = 0.000783642 0.000783642 
f(x) = 1.22819e-06
```

## 条件数对求解的影响
我们依次画出函数 $f(x) = x_1^2 + y^2$,   $f(x) = x_1^2 + 2* y^2$,   $f(x) = x_1^2 + 5 *y^2$ 和$f(x) = x_1^2 + 10* y^2$ 的梯度下降路径:
<div align="center">
<img src="./pics/gradient_py_1.png" width="300" height="300" alt="条件数1"> 
<img src="./pics/gradient_py_2.png" width="300" height="300" alt="条件数2"> 
<img src="./pics/gradient_py_5.png" width="300" height="300" alt="条件数5"> 
<img src="./pics/gradient_py_10.png" width="300" height="300" alt="条件数10"> 
</div>
可以看到，当条件数越大，梯度下降法的收敛速度越慢，当条件数越小时，梯度下降法的收敛速度越快。因为条件数越大，意味着在某个维度上，函数的变化越剧烈，在梯度下降法中，梯度越剧烈，下降的步长越大，因此梯度下降法的收敛速度越慢。

## Barzilai-Borwein 方法
当问题的条件数很大，也即问题比较病态时，梯度下降法的收敛性质会受到很大影响。
Barzilai-Borwein 方法是一种基于梯度下降法的方法，其基本思路是使用梯度下降法的方向和步长来更新参数，而不是使用一阶导数来更新参数。