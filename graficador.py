import numpy as np
import matplotlib.pylab as plt


datos = np.loadtxt("ejercicio28.dat")

xDrag = datos[:,0]
yDrag = datos[:,1]

plt.plot(xDrag,yDrag, color = "green", label = "Con Resistencia")
plt.grid()
plt.legend()
plt.savefig("ejercicio28.png")