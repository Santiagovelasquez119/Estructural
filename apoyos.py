import math

import matplotlib.pyplot as plt
import numpy as np


def apoyo_rigido(x1:float, y1:float, B:float):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(B*np.array([0, 1])+x1, B*np.array([0, 0])+y1, color='k')
    plt.plot(B*np.array([0, 0.5])+x1, B*np.array([0, 1])+y1, color='k')
    plt.plot(B*np.array([0.5, 1])+x1, B*np.array([1, 0])+y1, color='k')
    plt.plot(B*np.array([-0.1, 1.1])+x1, B*np.array([0, 0])+y1, color='k')
    for i in range(11):
        plt.plot(B*np.array([0.0 + 0.1 * i, -0.1 + 0.1 * i]), B*np.array([0, -0.1]), color='k')
    plt.show()

def apoyo_simple(B:float):
    R= 1
    plt.figure()
    theta = np.linspace(0, 2*np.pi,100)
    x=R*np.cos(theta)
    y=R*np.sin(theta)-1
    plt.plot(B*x,B*y, color='k')
    plt.plot(B*np.array([-1.1, 1.1]),B*np.array([-2,-2]), color='k')
    for i in range(22):
        plt.plot(B*np.array([-1.0+0.1*i, -1.1+0.1*i]), B*np.array([-2, -2.1]), color='k')
    plt.show()

apoyo_rigido(0.5)

def apoyo_empotrado(x1,y1,B):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(B * np.array([-1, 1.1]) + x1, B * np.array([-2, -2]) + y1, color='k')
    for i in range(22):
        plt.plot(B*np.array([-1.0+0.1*i, -1.1+0.1*i])+x1, B*np.array([-2, -2.1])+y1, color='k')
    plt.show()
apoyo_empotrado(0,0,1)