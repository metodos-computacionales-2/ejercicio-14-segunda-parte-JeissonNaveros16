import numpy as np
import matplotlib.pylab as plt

d1 = np.loadtxt("EulerSolution.dat")
d2 = np.loadtxt("RungeKuttaSolution.dat")
d3 = np.loadtxt("EDOFrictionSolution.dat")
d4 = np.loadtxt("EDOLambda3Solution.dat")
d5 = np.loadtxt("EDOLambda5Solution.dat")

plt.figure(figsize=(15,8))
plt.subplot(1,2,1)
#Punto 4
plt.title("Solución EDO - Relación posición y tiempo")
plt.plot(d1[:,0], d1[:,1], label="Con método de Euler")
plt.plot(d2[:,0], d2[:,1], label="Con método de Runge-Kutta")
plt.xlabel("Tiempo")
plt.ylabel("Posición")
plt.legend()

plt.subplot(1,2,2)
#Punto 6
plt.title("Solución EDO - Relación velocidad y posición")
plt.plot(d2[:,1], d2[:,2], label="Con método de Runge-Kutta")
plt.xlabel("Posición")
plt.ylabel("Velocidad")
plt.legend()

plt.savefig("SolucionEDO.png")

plt.figure()
#Punto 7
plt.title("Solución EDO con fricción - Relación posición y tiempo")
plt.plot(d3[:,0], d3[:,1], label="Con gamma=1")
plt.xlabel("Tiempo")
plt.ylabel("Posición")
plt.legend()

plt.savefig("SolucionEDOConFriccion.png")

#Punto 8
plt.figure(figsize=(15,8))
#Punto 6
plt.title("Solución a ecuación diferencial para distintos lambdas")
plt.plot(d2[:,0], d2[:,1], label="lambda=1")
plt.plot(d4[:,0], d4[:,1], label="lambda=3")
plt.plot(d5[:,0], d5[:,1], label="lambda=5")
plt.ylabel("Posición")
plt.xlabel("Tiempo")
plt.legend()

plt.savefig("SolucionEDDistintosLambdas.png")