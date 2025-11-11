
## 向量二范数的求导


$$ \begin{aligned} 
f(x) &= \frac{1}{2} || Ax - b ||^2  \\
     &= \frac{1}{2} (Ax - b)^T (Ax - b) \\
     &= \frac{1}{2} \{x^T A^T A x - x^T A^T b - b^T A x + b^T b\} \\
\end{aligned}$$
利用以下的性质：
$$ \frac{\partial x^T D x}{\partial x} = D x + x D^T$$
$$ \frac{\partial x^T D}{\partial x} = D y$$
$$ \frac{\partial D^T x }{\partial x} = D$$

可以得到：
$$ \begin{aligned} 
\frac{\partial f}{\partial x} 
    &= \frac{1}{2} (A^T A + (A^T A)^T)x - \frac{1}{2} A^T b - \frac{1}{2}(b^T A)^T\\ 
    & = \frac{1}{2} (2A^T A )x - \frac{1}{2} A^T b - \frac{1}{2}A^T b\\
    &= A^TAx - A^Tb \\
    &= A^T (Ax - b)\\
\end{aligned}$$