import numpy as np
import matplotlib.pyplot as plt


plt.figure()
x = np.linspace(-12, 12, 100)
y = np.linspace(-12, 12, 100)
X, Y = np.meshgrid(x, y)
Z = X**2 + 10*Y**2
plt.contour(X, Y, Z, 20)

def f(x):
    return x[0]**2 + 10*x[1]**2

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
plt.title('f(x, y) = x^2 + 10 * y^2 Gradient Descent Path')
plt.show()