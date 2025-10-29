# Newton's Method for Optimization
梯度法仅仅依赖函数值和梯度的信息（即一阶信息），如果函数 f ( x ) 充分光滑，则可以利用二阶导数信息构造下降方向 dk ．牛顿类算法就是利用二阶导数信息来构造迭代格式的算法。  
对于连续二次可微的函数 $f(x)$， 在迭代点 $x^k$ 处， 泰勒展示如下：
$$f(x^k + d^k) = f(x^k) + g^k(d^k)+ \frac{1}{2}(d^k)^T H^k(d^k) + O(||d^k||^3)$$
其中 $g^k(d^k) = \bigtriangledown f(x^k)^T d^k$ ， $H^k$ 为 $f(x^k)$ 的二阶导数矩阵．

牛顿类算法目的是根据这个二阶近似来选取合适的下降方向 $d^k$ ．忽略高阶项 $O(||d^k||^3)$ ，则牛顿类算法迭代公式为：
$$f(x^k + d^k) \approx f(x^k) + \bigtriangledown f^k(x)(d^k)+ \frac{1}{2}(d^k)^T \bigtriangledown^2 f^k(x)(d^k)$$
将右侧视为 $d^k$ 的函数，求最适合的 $d^k$ ，即是求合适的 $d$ ，使得 $d$ 满足：
$$ \arg \min\limits_{d} f(x^k) + \bigtriangledown f^k(x)(d) + \frac{1}{2}d^T \bigtriangledown^2 f^k(x)(d)$$
可得： $d^k = -\bigtriangledown^2 f^k(x)^{-1} \bigtriangledown f^k(x)$
因此经典牛顿法的更新格式为：
$$x^{k+1} = x^k - \bigtriangledown^2 f^k(x)^{-1} \bigtriangledown f^k(x)$$
步长 $α^k$ 恒为 1，即可以不额外考虑步长的选取.

## 经典牛顿法

经典牛顿法，其迭代公式为：
$$x^{k+1} = x^k - \bigtriangledown^2 f^k(x)^{-1} \bigtriangledown f^k(x)$$
迭代公式中， $x^k$ 为当前迭代点， $\bigtriangledown^2 f^k(x)^{-1}$ 为当前迭代点 $x^k$ 的二阶导数矩阵的逆， $\bigtriangledown f^k(x)$ 为当前迭代点 $x^k$ 的梯度向量．

给出 代码实现：
```cpp
