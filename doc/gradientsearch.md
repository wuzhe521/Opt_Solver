# Gradient Based Search

梯度类算法，其本质是仅仅使用函数的一阶导数信息选取下降方向 $d^k$. 其中最基本的梯度类算法为梯度下降法，直接选取负梯度作为下降方向．
梯度下降法的方向选取非常直观，实际应用范围非常广，因此它在优化算法中的地位可相当于高斯消元法在线性方程组算法中的地位．

梯度下降法的演示代码如下：
```python
import numpy as np
import matplotlib.pyplot as plt
from opt_solver import gradient_descent
x = np.linspace(-5, 5, 100)

```