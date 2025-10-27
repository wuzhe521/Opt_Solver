# Linear Search

我们将求解 f ( x ) 的最小值点的过程比喻成下山的过程．假设一个人处于某点 x 处， f ( x ) 表示此地的高度，为了寻找最低
点，在点 x 处需要确定如下两件事情：第一，下一步该向哪一方向行走；第
二，沿着该方向行走多远后停下以便选取下一个下山方向．以上这两个因素
确定后，便可以一直重复，直至到达 f ( x ) 的最小值点．  
线搜索类算法的数学表述为：给定当前迭代点 $x^k$ ，首先通过某种算法选取向量 $d^k$，之后确定正数 $α^k$ ，则下一步的迭代点可写作
$x^{k+1} = x^k + α^k d^k$．

$d^k$ 是一个向量，表示下一步的搜索方向, 在求最小值问题上，要求$d^k$是一个下降方向，即$(d^k)^T*\bigtriangledown f(x^k) < 0$。$α^k$ 是一个正数，表示沿着 $d^k$ 方向行走的距离．

线搜索类算法的关键是如何选取一个好的方向 $d^k$ 和一个好的步长 $α^k$ ．

方向 $d^k$ 的选取可以通过梯度下降法来实现，即 $d^k = -\bigtriangledown f(x^k)$ ．
步长 $α^k$ 的选取有多种方法，首先，构造辅助函数 $g(α^k) = f(x^k + α^k d^k)$ , 然后选取 $α^k$ 使得 $g(α^k)$ 尽可能的小，如果是求最小值问题，则 $g(α^k)$ 的最小值点即为 $x^k$ 的最小值点．这种方法称为精确搜索．但是这一方法的缺点是，计算复杂度较高，且无法保证找到全局最小值点．
我们在实际应用中常使用近似搜索方法，步长$α^k$ 无需 精确求解，只需满足一些条件即可．这种方法也常成为近似搜索/非精确搜索．
常见的搜索方法有：
1. 固定步长 $α^k$ ．
2. Armijo 条件
3. Goldstein 条件
4. wolfe 条件
这些方法在计算复杂度上较低，且能有效找到局部最小值点。

这类算法的实现原理为：  
1. 选取一个初始点 $x^0$ 作为起点，并确定一个初始的搜索方向 $d^0$ ．
2. 沿着 $d^0$ 方向 按照某个步长 行走，直到遇到一个点 $x^1$．
3. 判断点 $x^1$, 是否满足<font size = 3, color = "red"> 某个停止条件 </font>，否则将 $x^0$ 更新为 $x^1$ ，并重新选取一个新的搜索方向 $d^1$ ．
4. 重复步骤 2 和 3 ，直到满足某个停止条件．


## Armijo 条件
线搜索类算法中，Armijo 条件是常用的搜索条件．即步长满足如下条件：
$f(x^k + α^k d^k) < f(x^k) + c_1 α^k \bigtriangledown f(x^k)^T d^k$
其中 $c_1$ 是一个正数，常取 $c_1 = 10^{-4}$ 。
在gons中，Armijo 搜索条件类定义如下, 其中选取下降方向 $d^k$ 为 $d^k = -\bigtriangledown f(x^k)$ ．：
```cpp
X ArmijoRule() {
    const double min_beta =
        1e-10; // Minimum alpha to avoid numerical instability
    GONS_UINT iteration = 0;
    while (iteration < params_.max_inner_iter) {
      double beta = params_.beta;
      // Check if alpha is too small
      if (params_.beta < min_beta) {
        break;
      }
      // Cache gradient to avoid repeated computation
      X gradient = f_.gradient(x_);
      double f_x = f_(x_);
      double f_x_new = f_(x_ - beta * gradient);
      if (f_x_new < f_x - params_.alpha * beta * gradient.Norm2()) {
        // Armijo condition is satisfied
        params_.beta = beta;
        return x_ - beta * gradient;
      }
      // Update alpha for next iteration
      beta *= params_.gamma;
      if (params_.enable_max_iter)
        iteration++;
    }
    return x_ - params_.beta * f_.gradient(x_);
  }
  ```
得到步长以后，我们就需要更新迭代点了，在gons中，更新迭代点的方法如下：
```cpp
ArmijoStatus Optimize() {
    LOG(SOLVER_HEADER);
    GONS_UINT i = 0;
    while (true) {
      ++i;
      if (params_.enable_max_iter && i >= params_.max_iter) {
        LOG_ERROR("Max iterations reached: " << params_.max_iter);
        break;
      }
      X x_new = ArmijoRule();
      double f_old = f_(x_);
      double f_new = f_(x_new);
      double grad_norm = f_.gradient(x_).Norm2();
      if (params_.print_info) {
        LOG("Iteration: " << i);
        LOG_WARNING("Function value last: " << f_old);
        LOG_WARNING("Function value new: " << f_new);
        LOG_WARNING("Gradient Norm2: " << grad_norm);
        LOG_ERROR("x: " << x_);
      }
      if (std::abs(f_new - f_old) < params_.epsilon) {
        LOG_ERROR("After " << i << " iterations, ");
        LOG_ERROR("Armijo search converged to a local minimum");
        return ArmijoStatus::Success;
      }
      x_ = x_new;
    }
    return ArmijoStatus::Failure;
  }
  ```
以函数 $f(x_1, x_2) = x^2_1 + x^2_2$ 为例：
```cpp
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
```
设初始点为 $x^0 = (10, 10)$ , 搜索条件为 $c_1 = 10^{-4}$ , 搜索参数为 $max\_iter = 1000$ , $beta = 0.4$ , $gamma = 0.333$ ：
```cpp
int main()
{ 
    TestFunction f;
    X x = {10, 10};
    gons::ArmijoSearch<TestFunction, X> armijo_search(f, x);

    gons::ArmijoSearch<TestFunction, X>::ArmijoParameters params;

    params.alpha = 0.0001;
    params.beta = 0.4;
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
```
运行结果为：
> After 7 iterations,   
> Armijo search converged to a local minimum  
> Final x: 0.00064, 0.00064  
> Final f(x): 8.192e-07  

## Goldstein 搜索条件
线搜索类算法中，Goldstein 搜索条件是常用的搜索条件．即步长满足如下条件：
$f(x^k + α^k d^k) < f(x^k) + c_1 α^k \bigtriangledown f(x^k)^T d^k$
$f(x^k + α^k d^k) > f(x^k) + c_2 α^k \bigtriangledown f(x^k)^T d^k$
其中 $c_1$ 和 $c_2$ 都是正数，常取 $c_1 = 10^{-4}$ ， $c_2 = 0.9$ 。
在gons中，Goldstein 搜索条件类定义如下, 其中选取下降方向 $d^k$ 为 $d^k = -\bigtriangledown f(x^k)$ ．：
```cpp
X Search() {
    // Implementation of Goldstein search
    const double min_beta =
        1e-10; // Minimum beta to avoid numerical instability
    GONS_UINT iteration = 0;
    double beta = params_.beta;
    while (true) {
      ++iteration;
      if (iteration > params_.max_inner_iter) {
        break;
      }
      if (beta < min_beta) { // Check if beta is too small and break if it is
        break;
      }
      X gradient = f_.gradient(x_);
      double f_x = f_(x_);
      double f_x_new = f_(x_ - params_.beta * gradient);
      if (f_x_new < f_x - params_.beta *  params_.alpha * gradient.Norm2() &&
          f_x_new > f_x - (1 -  params_.alpha) * params_.beta * gradient.Norm2()) {
        params_.beta = beta;
        return x_ - params_.alpha * gradient;
      }
      if (f_x_new <=
          f_x - (1 - params_.alpha) * params_.beta * gradient.Norm2()) {
        //
        beta *= params_.gamma * 1.5;
        continue;
      }
      beta *= params_.gamma;
    }
    return x_ - params_.alpha * f_.gradient(x_);
  }
  ```
  得到步长以后，我们就需要更新迭代点了，在gons中，更新迭代点方法如下：
  ```cpp
GoldsteinStatus Optimize() {
    LOG(SOLVER_HEADER);
    GONS_UINT i = 0;
    while (true) {
      ++i;
      if (params_.enable_max_iter)
        if (i >= params_.max_iter)
          break;
      X x_new = Search();
      if (params_.print_info) {
        LOG("Iteration: " << i);
        LOG_WARNING("Function value last: " << f_(x_));
        LOG_WARNING("Function value new: " << f_(x_new));
        LOG_WARNING("Gradient Norm2: " << f_.gradient(x_).Norm2());
        LOG_ERROR("x: " << x_);
      }
      if (std::abs(f_(x_new) - f_(x_)) < params_.epsilon) {
        LOG_ERROR("After " << i << " iterations, ");
        LOG_ERROR("Goldstein search converged to a local minimum");
        return GoldsteinStatus::Success;
      }
      x_ = x_new;
    }
    return GoldsteinStatus::Failure;
  }
  ```
  以函数 $f(x_1, x_2) = x^2_1 + x^2_2$ 为例, 我们可以使用Goldstein搜索算法来找到其最小值。
  ```cpp
  gons::GoldsteinSearch<TestFunction, X> goldstein_search(f, x);
  gons::GoldsteinSearch<TestFunction, X>::GoldsteinParameters goldstein_params;
  goldstein_search.Optimize();
  LOG("Final x: " << goldstein_search.get_x());
  LOG("Final f(x): " << goldstein_search.get_function_value());
  ```
  运行结果为：
  > After 15 iterations,  
  > Goldstein search converged to a local minimum   
  > Final x: 0.00061, 0.000641   
  > Final f(x): 7.45058e-07

  ## wolfe 搜索条件
  线搜索类算法中，wolfe 搜索条件是常用的搜索条件．即步长满足如下条件：  
  1. $f(x^k + α^k d^k) < f(x^k) + c_1 α^k \bigtriangledown f(x^k)^T d^k$  
  2. $\bigtriangledown f(x^k + α^k d^k)^T d^k > c_2 \bigtriangledown f(x^k)^T d^k$   
  
  其中 $c_1$ 和 $c_2$ 都是正数，常取 $c_1 = 10^{-4}$ ， $c_2 = 0.9$ 。可以看出，wolfe 搜索条件是Armijo 搜索条件 的一个扩展,
  wolfe 条件的第一条就是 Armijo 条件，第二条是 curvature 条件, 要求在新的点 $x^k + α^k d^k$ 处，切线的斜率大于 $c_2 \bigtriangledown f(x^k)^T d^k$ 。
  在gons中，wolfe 搜索条件类定义如下, 其中选取下降方向 $d^k$ 为 $d^k = -\bigtriangledown f(x^k)$ ．：
  ```cpp
  X Search() {
    const double min_beta = 1e-10;
    GONS_UINT iteration = 0;
    double beta = params_.beta;
    while (true) {
      ++iteration;
      if (iteration > params_.max_inner_iter) {
        break;
      }
      if (beta < min_beta) { // Check if beta is too small and break if it is
        break;
      }
      X gradient = f_.gradient(x_);
      double f_x = f_(x_);
      // Wolfe conditions
      double f_x_new = f_(x_ - beta * gradient);
      if (f_x_new <= f_x - params_.alpha_1 * beta * gradient.Norm2()) {
        // Armijo condition satisfied
        X gradient_new = f_.gradient(x_ - beta * gradient);
        if (-gradient_new.Norm2() > - params_.alpha_2 * gradient.Norm2()) {
          // Wolfe condition satisfied
          params_.beta = beta;
          return x_ - beta * gradient;
        } else {
          // Only Armijo condition satisfied, decrease beta
          beta *= params_.gamma;
        }
      } else {
        // Neither condition satisfied, increase beta
        beta *= params_.beta;
      }
    }
    params_.beta = beta;
    return x_ - beta * f_.gradient(x_);
  }
  ```
  得步长以后，我们就需要更新迭代点了，在gons中，更新迭代点方法如下：
  ```cpp
  WolfeStatus Optimize() {
    LOG(SOLVER_HEADER);
    GONS_UINT i = 0;
    while (true) {
      ++i;
      if (params_.enable_max_iter)
        if (i >= params_.max_iter)
          break;
      X x_new = Search();
      if (params_.print_info) { 
        LOG("Iteration: " << i);
        LOG_WARNING("Function value last: " << f_(x_));
        LOG_WARNING("Function value new: " << f_(x_new));
        LOG_WARNING("Gradient Norm2: " << f_.gradient(x_).Norm2());
        LOG_ERROR("x: " << x_);
      }
      if (std::abs(f_(x_new) - f_(x_)) < params_.epsilon) {
        LOG_ERROR("After " << i << " iterations, ");
        LOG_ERROR("Wolfe search converged to a local minimum");
        return WolfeStatus::Success;
      }
      x_ = x_new;
    }
    return WolfeStatus::Failure;
  }
  ```
  以函数 $f(x_1, x_2) = x^2_1 + x^2_2$ 为例, 我们可以使用wolfe搜索算法来找到其最小值。
  ```cpp
  gons::WolfeSearch<TestFunction, X> wolfe_search(f, x);
  gons::WolfeSearch<TestFunction, X>::WolfeParameters wolfe_params;
  wolfe_search.Optimize();
  LOG("Final x: " << wolfe_search.get_x());
  LOG("Final f(x): " << wolfe_search.get_function_value());
  return 0;
  ```
  运行结果为：
 > After 20 iterations,   
 > Wolfe search converged to a local minimum  
 > Final x: 0.00060936 0.00060936   
 > Final f(x): 7.42639e-07  