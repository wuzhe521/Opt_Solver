# Newton's Method for Optimization

## 牛顿法的基本思想
梯度法仅仅依赖函数值和梯度的信息（即一阶信息），如果函数 f ( x ) 充分光滑，则可以利用二阶导数信息构造下降方向 dk ．牛顿类算法就是利用二阶导数信息来构造迭代格式的算法。  
对于连续二次可微的函数 $f(x)$， 在迭代点 $x^k$ 处， 泰勒展示如下：
$$f(x^k + d^k) = f(x^k) + g^k(d^k)+ \frac{1}{2}(d^k)^T H^k(d^k) + O(||d^k||^3)$$
其中 $g^k(d^k) = \bigtriangledown f(x^k)^T d^k$ ， $H^k$ 为 $f(x^k)$ 的二阶导数矩阵．

牛顿类算法目的是根据这个二阶近似来选取合适的下降方向 $d^k$ ．忽略高阶项 $O(||d^k||^3)$ ，则牛顿类算法迭代公式为：
$$f(x^k + d^k) \approx f(x^k) + \bigtriangledown f(x^k)(d^k)+ \frac{1}{2}(d^k)^T \bigtriangledown^2 f(x^k)(d^k)$$
将右侧视为 $d^k$ 的函数，求最适合的 $d^k$ ，即是求合适的 $d$ ，使得 $d$ 满足：
$$ \arg \min\limits_{d} f(x^k) + \bigtriangledown f(x^k)(d) + \frac{1}{2}d^T \bigtriangledown^2 f(x^k(d)$$
可得： $d^k = -\bigtriangledown^2 f(x^k)^{-1} \bigtriangledown f(x^k)$
因此经典牛顿法的更新格式为：
$$x^{k+1} = x^k - \bigtriangledown^2 f(x^k)^{-1} \bigtriangledown f(x^k)$$
步长 $α^k$ 恒为 1，即可以不额外考虑步长的选取.

## 经典牛顿法

经典牛顿法，其迭代公式为：
$$x^{k+1} = x^k - \bigtriangledown^2 f^k(x)^{-1} \bigtriangledown f^k(x)$$
迭代公式中， $x^k$ 为当前迭代点， $\bigtriangledown^2 f^k(x)^{-1}$ 为当前迭代点 $x^k$ 的二阶导数矩阵的逆， $\bigtriangledown f^k(x)$ 为当前迭代点 $x^k$ 的梯度向量．

给出 代码实现：
```cpp
template <typename Function, typename X> class NewtonMethod {
public:
  struct NewtonParameters {
    double tolerance = 1.0e-6;
    GONS_UINT max_iterations = 3;
    bool verbose = true;
  };
  enum class NewtonStatus { Success, Failure, MaxIterationReached };

public:
  NewtonMethod(const Function &f, const X &x) : f_(f), x_(x) {}

  // compute search direction
  X Search(const X &x, const Function &f) {
    X gradient = f.gradient(x);
   
    gons::Matrix hessian = f.hessian(x);
    
    // gons::Matrix hessian_inv = hessian.Inverse();
    X neg_gradient = gradient * -1;
    // Solve for step: hessian * step = -gradient
    auto step = SolveLinearSystem(hessian, neg_gradient);
    return step;
  }

  void Optimize() {
    for (int i = 0; i < param_.max_iterations; ++i) {
      X gradient = Search(x_, f_);
      x_ += gradient;
      if (gradient.Norm2() < param_.tolerance) {
        LOG("After " << i << " iterations, ");
        LOG("Newton method converged.");
        LOG("Final solution: " << x_);
        LOG("Final gradient: " << gradient);
        LOG("Final function value: " << f_(x_));
        return;
      }
      if (param_.verbose) {
        LOG("Iteration: " << i);
        LOG("x = " << x_);
        LOG("f(x) = " << f_(x_));
        LOG("Gradient: " << gradient);
      }
    }
    LOG("After " << param_.max_iterations << " iterations, ");
    LOG("Newton method did not converge.");
  }

private:
  Function f_;
  X x_;
  NewtonParameters param_;
};
```
测试函数仍采用 $f(x) = x^2 + 10y^2 $,   测试代码：
```cpp
TestFunction f;
X x = {10, 10};
gons::NewtonMethod<TestFunction, X> nm(f, x);
nm.Optimize();
```
输出结果：
```bash
After 1 iterations, 
Newton method converged.
Final solution: 0 0 
Final gradient: 0 0 
Final function value: 0
```
测试结果表明，Newton方法在给定初始点x = [10, 10]时，在1次迭代后收敛，最终的解为[0, 0]，并且函数值最小为0。


## 修正牛顿法

### 修正牛顿法的基本思想
[经典牛顿法存在缺陷](http://faculty.bicmr.pku.edu.cn/~wenzw/optbook/opt1-short.pdf)：
1. 每一步迭代需要求解一个 n 维线性方程组，这导致在高维问题中计算量很大．海瑟矩阵 $H^k$ 既不容易计算又不容易储存．
2. 牛顿法的迭代方向 $d^k$ 是基于 $H^k$ 的逆矩阵，如果 $H^k$ 的逆矩阵不存在，则迭代方向不存在，迭代过程无法继续．
3. 当迭代点距最优值较远时，直接选取步长 αk = 1 会使得迭代极其不稳定，在有些情况下迭代点列会发散. 

带线搜索的修正牛顿法，其基本思想是对牛顿方程中的海瑟矩阵 $H^k$ 进行修正，使其变成**正定矩阵**；同时引入线搜索以改善算法稳定性。
### 修正矩阵的选取
矩阵 $E^k$ 确定，一般有两种方法：
1. 单位矩阵放大；选取 $E^k$ 为 $\lambda I$, 其中 $\lambda$ 为正定常数，$I$ 为单位矩阵，$\lambda$ 的选取需要根据具体问题确定，不能只是一味的增大。

2. 修正Cholesky分解; 即计算 $ \bigtriangledown^2 f^k(x^k) + E $，其中 $E$ 为修正矩阵。

### 带线搜索的修正牛顿法  
引入修正矩阵使得： $B = \bigtriangledown^2 f^k(x^k) + E$为正定矩阵，且条件数最小。  
一般的框架算法：

1. 给定初始点 $x_0$
2. 循环迭代：
4. 求解修正的方程 $B^k d^k = - \bigtriangledown^k f^k(x^k)$, 得到 $d^k$
5. 线搜索确定步长 $\alpha^k$
6. 更新 $x^{k+1} = x^k + \alpha^k d^k$
7. 如果 $||d^k|| < \epsilon$, 则停止迭代，否则返回 2

### 代码实现
<font color="red">To be continued. </font>


## 非精确牛顿法
在经典牛顿法中，计算牛顿方向 dk 依赖于求解线性方程组， 当海瑟矩阵 Hk 不易计算或不易储存时， 非精确牛顿法可以有效地解决这个问题．

常用的非精确牛顿法是牛顿共轭梯度法．由于共轭梯度法在求解线性方程组方面有不错的表现，因此在多数问题上牛顿共轭梯度法都有较好的数值表现。

<font color="red">To be continued. </font>