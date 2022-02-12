import numpy as np
import matplotlib.pyplot as plt
import PDEsolver as pde
import matplotlib.animation as animation

x0 = 0
xn = 1
dx = 1/16
Nx = int((xn-x0)/dx)
x = np.linspace(x0,xn,Nx+1)

t0 = 0
tn = 1
dt = 1/800
Nt = int((tn-t0)/dt)
t = np.linspace(t0,tn,Nt+1)

r = dt/(dx)**2

def f(x):
    return 1 + 2*x + np.sin(2*np.pi*x)

def g1(t):
    return 1

def g2(t):
    return 2

time,x,heat = pde.FTCS(x0,xn,dx,t0,tn,dt,f,g1,g2,"DN")

fig1,ax1 = plt.subplots()

def animate(i):
    ax1.clear()
    
    ax1.plot(x,heat[i,:],"r--")
    ax1.plot(x,heat[i,:],"ko")
    ax1.set_title("time = %1.2f sec"%(i))
    ax1.grid()
    ax1.set_xlabel("x")
    ax1.set_ylabel("U(x)")

ani = animation.FuncAnimation(fig1,animate,50,interval=10, blit=False)
#plt.show()
ani.save('Explicit-heat.gif', writer='imagemagick', fps=5)
