import numpy as np
import matplotlib.pylab as plt

datos = np.loadtxt("ejercicio28.dat")

xDrag = datos[:,0]
yDrag = datos[:,1]
xSDrag = datos[:,2]
ySDrag = datos[:,3]

plt.plot(xDrag,yDrag, color = "green", label = "Con Resistencia")
plt.plot(xSDrag,ySDrag, color = "red", label = "Sin Resistencia")
plt.title("Proyectil")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid()
plt.legend()
plt.savefig("friccion.png")