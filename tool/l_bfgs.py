import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-12, 12, 100)
y = np.linspace(-12, 12, 100)
X, Y = np.meshgrid(x, y)

Z = X**2 + 10*Y**2
plt.contour(X, Y, Z, 20)
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Contour plot of the x^2 + 10y^2")
plt.show()

