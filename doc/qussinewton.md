# Quasi-Newton Methods
牛顿法在理论上和实践中均取得很好的效果．然而对于大规模问题，函数的海瑟矩阵计算代价特别大或者难以得到，即便得到海瑟矩阵我们还需要求解一个大规模线性方程组．那么能否使用海瑟矩阵或其逆矩阵的近似来进行牛顿迭代呢？拟牛顿法便是这样的算法，它能够在每一步以较小的计算代价生成近似矩阵，并且使用近似矩阵代替海瑟矩阵而产生的迭代序列仍具有超线性收敛的性质．  
拟牛顿方法不计算海瑟矩阵 $\nabla^2 f(x)$，而是构造其近似矩阵 $B^k$ 或其逆的近似矩阵 $H^k$ ．我们希望 $B^k$或$H^k$仍然保留海瑟矩阵的部分性质，例如使得 $d^k 仍然为下降方向．

在牛顿法中已知，在第k 步迭代中， $x^k$ 处， 梯度函数 $\nabla f(x)$ 满足：
$$f(x^k + d^k) \approx f(x^k) + \bigtriangledown f(x^k)(d^k)+ \frac{1}{2}(d^k)^T \bigtriangledown^2 f(x^k)(d^k)$$
步长 $d^k$ 最佳的选取是：$d^k = -\bigtriangledown^2 f(x^k)^{-1} \bigtriangledown f(x^k)$

拟牛顿法的理论基础就是找到一个矩阵 $B^k$，使得 $B^k d^k = - \bigtriangledown f(x^k)$

# 割线法



连续二阶可微函数 $f(x)$ 的 梯度函数 $\nabla f(x)$， 在 $x^{k+1}$ 处， 泰勒展式为：
$$
\nabla f(x) = \nabla f(x^{k+1})^T  + \nabla^2 f(x^{k+1}) (x - x^{k+1}) + O(||x - x^{k+1}||^2)
$$
转换一下：
$$
\nabla^2 f(x^{k+1}) (x - x^{k+1}) + O(||x - x^{k+1}||^2) = \nabla f(x) - \nabla f(x^{k+1})^T (x - x^k)
$$
令 $x = x^k$ ， $s^k = x^{k+1} - x^k$，$ y^k = \nabla f(x^{k+1}) - \nabla f(x^k)^T$，则可以得到：
$$
\nabla^2 f(x^{k+1}) s^k + O(||s^k||^2) = y^kx
$$

忽略高阶项：
$$
\nabla^2 f(x^{k+1}) s^k \approx y^kx
$$

我们构造的$f(x)$的 Hessian矩阵的近似矩阵$B$, 和近似逆矩阵$H$，仍满足：  
$$
B^k s^k = y^kx \\

s^k = H^k y^k  \tag{1}
$$

方程（1） 也被成为割线方程。


另外一个角度看割线方程：  
牛顿法的本质是对 $f(x)$在 点$x^k$处二阶泰勒展开， 然后根据展开后的公式，固定 $f(x)$ 和 $\nabla f(x)$， 求解 $d^{k}$ ， 求使得 $f(x+ d^{k})$ 最小化的 $d^{k}$  
当我们用割线法，用 $B^k$近似代替 $ \nabla^2 f(x^{k+1})$ 后， 将整个步长函数$M(d)$ 在点$k+1$ 的二阶泰勒近似：
$$
M_{k+1}(d) = f(x^{k+1}) + \nabla f(x^{k+1})^T d^k + \frac{1}{2}(d^k)^T B^k d^k
$$
$$
\nabla M_{k+1}(d) = \nabla f(x^{k+1}) + B^k d^k
$$
我们要求 在 $M_{k+1}(d)$ 在 $d = -s^k$和 $d = 0$ 处的梯度，和 $f(x)$ 在 $x=x^k$和$x = x^{k+1}$处的梯度保持一致。
$$
\nabla M_{k+1}(-s^k) = M_{k+1}(-s^k) =  \nabla f(x^{k+1}) + B^k s^k  = \nabla f(x^k) \\
\nabla M_{k+1}(0) = \nabla f(x^k) 
$$
可以看出，$\nabla M_{k+1}(0)$ 必然满足，主要 $\nabla f(x^{k+1}) + B^k s^k  = \nabla f(x^k) $成立，即可满足要求。
可得：
$$
B^k s^k  = \nabla f(x^k) - \nabla f(x^{k+1})  = y^k
$$

## 算法框架
1. 给定初始点 $x_0$， $B^0$ 或者 $H^0$
2. 循环迭代：
3. 计算方向 $ d^k = -(B^k)^{-1} \nabla f(x^k)$ 或者 $d^k = -H^k \nabla f(x^k)$
4. 线搜索确定步长 $\alpha^k$
5. 更新 $x^{k+1} = x^k + \alpha^k d^k$
6. 更新 $B^{k+1}$ 或者 $H^{k+1}$
7. 判断是否收敛，如果 $||d^k|| < \epsilon$, 则停止迭代，否则返回 2


在实际应用中基于$H^{k}$的拟牛顿法因其无需求解线性方程组，更加实用。基于 $B^{k}$ 的方法稳定性更好，但计算量更大。

那么，剩下的问题就成了如何选取 $B^{k+1}$ 或者 $H^{k+1}$ ？

## 拟牛顿矩阵更新形式

### [秩一更新](https://zhuanlan.zhihu.com/p/144736223)

这个方法继承的思路是：我希望在满足割线方程的条件下，每一步对于$B^k$的更新能够尽可能的小。这个“尽可能小”得含义，其实就是只修改矩阵的一个特征值：
$$
\begin{aligned}
&B^{k+1} = \min ||B- B^k||^2 \\
&s.t.\\
&B^{k+1} s^k = y^k
\end{aligned}
$$
有一种矩阵的变换方式叫做秩一更新，$ A + auu^T $, 在正定矩阵$A的基础上，加上一个秩为1的矩阵。
$$
B^{k+1} = B^k + auu^T \\
$$
其中, $auu^T$ 是一个正定，秩为1的矩阵， $a$ 是一个正数， $u$ 是一个单位向量。$a$和$u$是待定矩阵。
这种方法的好处是，可以套用 [Sherman-Morrison公式](https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula)来求出 $B^{k+1}$ 的逆矩阵。将求逆 O(n^3) 的复杂度，降低为 O(n^2)。但前提是，海瑟矩阵 $B^k$ 是正定的 ，且满足割线条件。

根据割线方程：
$$
B^{k+1}s^k = B^ks^k + \alpha u  u^T s^k = -y^k \tag{1}
$$
$$
\alpha u^T  u  s^k = y^k - B^ks^k \tag{2}
$$
令 $u = y^k - B^ks^k$ ，则 (2) 可以写成：
$$
\alpha (y^k - B^ks^k)^Ts^k(y^k - B^ks^k) = y^k - B^ks^k \tag{3}
$$
得到， $\alpha$的更新格式:
$$
\alpha = \frac{(y^k - B^ks^k)(y^k - B^ks^k)^T}{(y^k - B^ks^k)^Ts^k}
$$
同理，可以得到基于$H^k$的更新格式：
$$
H^{k+1} = H^k + \frac{(s^k - H^ky^k)(s^k - H^ky^k)^T}{(s^k - H^ky^k)^T y^k}  \tag{4}
$$

SR1 公式虽然结构简单，但是有一个重大缺陷：它不能保证矩阵在迭代过程中保持正定