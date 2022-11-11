import numpy as np
import math 

nx=10
nt=100
h=1/nx
k=1/nt
r=k/(h*h)
alpha = 1
dt=100
t=np.arange(0,(dt+.5)*k,k)
x=np.arange(0,1.0001,h)
w=np.zeros((nx+1,dt+1))
b=np.zeros(nx-1)

# Initial Condition
for i in range (1,nx):
    w[i,0]=x[i]*x[i]
    
# Boundary Condition
for j in range (0,dt+1):
    w[0,j]=t[j]
    w[nx,j]=2-math.exp(-t[j])

A=np.zeros((nx-1,nx-1))
B=np.zeros((nx-1,nx-1))
for i in range (0,nx-1):
    A[i,i]=2+2*r
    B[i,i]=2-2*r

for i in range (0,nx-2):           
    A[i+1,i]=-r
    A[i,i+1]=-r
    B[i+1,i]=r
    B[i,i+1]=r
    
Ainv=np.linalg.inv(A)   

for j in range (1,dt+1):
    b[0]=r*w[0,j-1]+r*w[0,j]
    b[nx-2]=r*w[nx,j-1]+r*w[nx,j]
    v=np.dot(B,w[1:nx,j-1])
    w[1:nx,j]=np.dot(Ainv,v+b)

ue = np.dot(np.exp(np.dot(-t[dt],dt*alpha*(math.pi/1.0001)**2)), np.sin(np.dot(math.pi/1.0001,x)))
#err = np.linalg.norm(w[:,dt] - ue)
print(w)
