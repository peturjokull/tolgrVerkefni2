import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

c=299792.458; #Ljoshradi km/s
r1=[15600,7540,20140,0.07074]; #Gervitungl 1   
r2=[18760,2750,18610,0.07220] #Gervitungl 2
r3=[17610,14630,13480,0.07690] #Gervitungl 3
r4=[19170,610,18390,0.07242] #Gervitungl 4
r0=[0,0,6370,0.0001] #hnit mottakara
rho = 26570 #km
phi = [0, math.pi/8,math.pi/4, math.pi/4]
theta = [0, math.pi/2,math.pi,3*math.pi/2]
R=np.zeros(4)
t_m=np.zeros(4)
A=[r1[0],r2[0],r3[0],r4[0]]
B=[r1[0],r2[1],r3[2],r4[3]]
C=[r1[0],r2[1],r3[2],r4[3]]
t_m = [r1[3],r2[3],r3[3],r4[3]]

for i in range(4):
    A[i] = rho*math.cos(phi[i])*math.cos(theta[i])
    B[i] = rho*math.cos(phi[i])*math.sin(theta[i])
    C[i] = rho*math.cos(phi[i])
    R[i] = math.sqrt(A[i]**2 + B[i]**2 + (C[i]-r0[3])**2)
    t_m[i]=r0[3] + R[i]/c
print(A,'\n',B,'\n',C,'\n',R)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(r0[0],r0[1],r0[2],c='b',s=1000)
ax.scatter(A[0],A[1],A[2],c='r')
ax.scatter(B[0],B[1],B[2],c='r')
ax.scatter(C[0],C[1],C[2],c='r')
ax.scatter(R[0],R[1],R[2],c='r')
# ax.scatter(19170,610,18390,c='r')
ax.set_xlabel('x-ás (km)')
ax.set_ylabel('y-ás (km)')
ax.set_zlabel('z-ás (km)')

ax.invert_xaxis()

plt.show()

