import numpy as np
import matplotlib.pyplot as plt


def barrier_function(x, t):
    return (-t)* np.log10(x)
def barrier_gradient(x):
    return -1 / (1 + np.exp(x))

if __name__ == '__main__':
    plt.xlim(-1, 12)
    plt.xlabel('x')
    plt.ylabel('Barrier function : -t*log10(-1.0 * x)')
    # plt.ylim(-5, 5)
    x = np.linspace(0.01 ,10, 100)
    for i in range(10):
        y = barrier_function(x, i)
        plt.plot(x, y, label=f'factor t={i}')
    plt.legend()
    plt.show()
    

