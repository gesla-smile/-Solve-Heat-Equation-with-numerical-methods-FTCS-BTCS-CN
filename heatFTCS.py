import numpy as np
import math as math

# Set initial condition
L = 1  # length of domain in x direction
tmax = 0.5  # end time
nx = 20  # number of nodes in x direction
nt = 10  # number of time steps
alpha = 0.1 #diffusivity, set it defaultly as 0.1

#compute mesh spacing and time step
dx = L / (nx - 1)
dt = tmax / (nt - 1)
r = alpha * dt / dx ** 2
r2 = 1 - 2 * r

#create arrays to save data to export
x = np.linspace(0, L, nx)
t = np.linspace(0, tmax, nt)
U = np.zeros((nx, nt))
ue = np.linspace(0, L, nx)


#Set initial condion and boundry condition
U[:,0] = (np.sin(np.dot(math.pi*1/L,x)))

#Loop over time steps
for m in range(2, nt):
    for i in range(2, nx-1):
        U[i,m] = r * U[i-1, m-1] + r2*U[i, m-1] + r*U[i+1, m-1]

#Compare with exact solution at end of simulation
#ue = math.sin(math.pi*x/L)*math.exp(-t*nt*alpha*(math.pi/L)**2
ue = np.dot(np.exp(np.dot(-t[nt-1],nt*alpha*(math.pi/L)**2)), np.sin(np.dot(math.pi/L,x)))
err = np.linalg.norm(U[:,nt-1] - ue)
print(err)