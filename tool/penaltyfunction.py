# 导入包
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# 支持中文
matplotlib.rcParams["font.sans-serif"] = ["SimSun"]
matplotlib.rcParams["axes.unicode_minus"] = False
# 渲染公式
plt.rcParams["text.usetex"] = True


# 定义原函数
def f(X):
    x = X[0]
    y = X[1]
    return x**2 + 5 * y**2 + x * y


# 定义约束函数
def h(X):
    x = X[0]
    y = X[1]
    return np.array([x**2 + y**2 - 1])


# 定义梯度函数
def grad_f(X):
    x = X[0]
    y = X[1]
    return np.array([2 * x + y, 10 * y + x])


# 定义模长函数
def norm(vector):
    return np.sqrt(sum([i**2 for i in vector]))


# 回溯线性搜索实现梯度下降
alpha = 0.4
beta = 0.5
epsilon = 1e-6
k = 0

X = np.array([2, 2])
t = 1

k_list = [k]
X_all = X
# 定义罚函数
def g(X):
    return max(0, h(X)) ** 2

# 定义罚函数的梯度
def grad_g(X):
    if h(X) > 0:
        return np.array([2 * X[0], 2 * X[1]])
    else:
        return np.array([0, 0])


# 定义加入罚函数之后的目标函数
def penalty(X, u):
    return f(X) + u * g(X)


# 定义加入罚函数之后的目标函数的梯度
def grad_penalty(X, u):
    return grad_f(X) + u * grad_g(X)

# 初始化参数
u = 1
k = 0

X = np.array([2, 2])
t = 1

k_list = [k]
X_all = X
# 当梯度的模仍然高于精度要求，或者约束函数不被满足时，需要一直迭代
# 通过无约束优化算法求解加入罚函数之后的目标函数的极小值点
while norm(grad_f(X)) > epsilon or h(X) > epsilon:
    while norm(grad_penalty(X, u)) > epsilon:
        # 回溯线性搜索，检查目前的步长是否满足要求，若步长太大，则需要缩小步长
        while penalty(np.subtract(X, t * grad_penalty(X, u)), u) > penalty(
            X, u
        ) - alpha * t * grad_penalty(X, u).T @ grad_penalty(X, u):
            t = beta * t
        # 更新自变量值
        X = X - t * grad_penalty(X, u)
        # 记录此次的自变量值，它是加入罚函数之后的目标函数的极小值点
        k = k + 1
        k_list.append(k)
        X_all = np.vstack((X_all, X))
    # 增大惩罚参数
    u = 10 * u
# 退出循环，说明已经得到了满足精度要求的梯度，并且也满足约束条件。此时的自变量值就是极小值点
print("极小值点为：x* = {:.2f}, y* = {:.2f}".format(X[0], X[1]))


