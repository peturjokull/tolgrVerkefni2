import numpy.linalg as lin
import numpy as np
import itertools
import matplotlib.pyplot as plt
c = 299792.458  # Ljoshradi km/s
r0 = np.asarray([0, 0, 6370, 0.0001])  # hnit mottakara
rho = 26570  # km

# Aðferð til að fá hnit dreifð jafnt yfir jörðina
def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec.T

R = sample_spherical(8)

dist=np.zeros(len(R))
# print(R)
for i in range(len(R)):
    dist[i] = np.sqrt(((R[i,0]**2 + R[i,1]**2 + (R[i,2]-6370)**2)))/c
# print(dist.shape)
A = R[:,0].T # x-hnit gervihnatta
B = R[:,1].T # y-hnit gervihnatta
C = R[:,2].T # z-hnit gervihnatta
t_m = dist.T # fjarlægðar fylki

#### TEIKNA HER #### 
teikna = False
if teikna is True:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(r0[0], r0[1], r0[2], c='b', s=1000) #jörðin
    ax.scatter(A, B, C, c='r') # gervihnettir
    ax.set_xlabel('x-ás (km)')
    ax.set_ylabel('y-ás (km)')
    ax.set_zlabel('z-ás (km)')

    ax.invert_xaxis()

    plt.show()

####            ####

#errorcoef eru 
# allar mögulegu uppraðanir á gervihnöttunum
errorcoef = []
for i in itertools.product([0,1],repeat=8):
    errorcoef.append(i)
errorcoef = (np.asarray(errorcoef))

dt = 10**-8  # nákvæmni klukku í gervihnöttum
errormagcoef = 0
maxFE = 0

for i in range(0,len(errorcoef)):
    t = [t_m[0]+dt*errorcoef[i,0],
    t_m[1]+dt*errorcoef[i,1],
    t_m[2]+dt*errorcoef[i,2],
    t_m[3]+dt*errorcoef[i,3],
    t_m[4]+dt*errorcoef[i,4],
    t_m[5]+dt*errorcoef[i,5],
    t_m[6]+dt*errorcoef[i,6],
    t_m[7]+dt*errorcoef[i,7]
    ]

    t_diff = [
        t[0] - t_m[0],
        t[1] - t_m[1],
        t[2] - t_m[2],
        t[3] - t_m[3],
        t[4] - t_m[4],
        t[5] - t_m[5],
        t[6] - t_m[6],
        t[7] - t_m[7]
    ]

    for j in range(0,10):
        DF = np.array([
            [(r0[0] - A[0])*2, (r0[1] - B[0])*2,
                (r0[2] - C[0])*2, (2*c**2)*(t[0] - r0[3])],
            [(r0[0] - A[1])*2, (r0[1] - B[1])*2,
                (r0[2] - C[1])*2, (2*c**2)*(t[1] - r0[3])],
            [(r0[0] - A[2])*2, (r0[1] - B[2])*2,
                (r0[2] - C[2])*2, (2*c**2)*(t[2] - r0[3])],
            [(r0[0] - A[3])*2, (r0[1] - B[3])*2,
                (r0[2] - C[3])*2, (2*c**2)*(t[3] - r0[3])],
            [(r0[0] - A[4])*2, (r0[1] - B[4])*2,
                (r0[2] - C[4])*2, (2*c**2)*(t[4] - r0[3])],
            [(r0[0] - A[5])*2, (r0[1] - B[5])*2,
                (r0[2] - C[5])*2, (2*c**2)*(t[5] - r0[3])],
            [(r0[0] - A[6])*2, (r0[1] - B[6])*2,
                (r0[2] - C[6])*2, (2*c**2)*(t[6] - r0[3])],
            [(r0[0] - A[7])*2, (r0[1] - B[7])*2,
                (r0[2] - C[7])*2, (2*c**2)*(t[7] - r0[3])]    

        ])

        F = np.array([
        [(r0[0] - A[0])**2 + (r0[1] - B[0])**2 +
            (r0[2] - C[0])**2 - (c**2)*(t[0] - r0[3])**2],
        [(r0[0] - A[1])**2 + (r0[1] - B[1])**2 +
            (r0[2] - C[1])**2 - (c**2)*(t[1] - r0[3])**2],
        [(r0[0] - A[2])**2 + (r0[1] - B[2])**2 +
            (r0[2] - C[2])**2 - (c**2)*(t[2] - r0[3])**2],
        [(r0[0] - A[3])**2 + (r0[1] - B[3])**2 +
            (r0[2] - C[3])**2 - (c**2)*(t[3] - r0[3])**2],
        [(r0[0] - A[4])**2 + (r0[1] - B[4])**2 +
            (r0[2] - C[4])**2 - (c**2)*(t[4] - r0[3])**2],
        [(r0[0] - A[5])**2 + (r0[1] - B[5])**2 +
            (r0[2] - C[5])**2 - (c**2)*(t[5] - r0[3])**2],
        [(r0[0] - A[6])**2 + (r0[1] - B[6])**2 +
            (r0[2] - C[6])**2 - (c**2)*(t[6] - r0[3])**2],
        [(r0[0] - A[7])**2 + (r0[1] - B[7])**2 +
            (r0[2] - C[7])**2 - (c**2)*(t[7] - r0[3])**2]
        ])
        # print('DF.shape',DF.shape)
        # print('F.shape',F.shape)
        v = lin.solve(np.dot(DF.T,DF) ,np.dot(-DF.T,F))
        # print('v',v.T.shape)
        # print('r',r0.shape)
        for i in range(len(r0)):
            r0[i] = r0[i] + v[i,0]
    
    forwarderror = lin.norm(np.asarray([r0[0], r0[1], r0[2]-6370]), ord=2)
    backwarderror = c*lin.norm([t[0]-t_m[0], t[1] -
                                t_m[1], t[2]-t_m[2], t[3]-t_m[3]], ord=2)

    if maxFE <= np.absolute(forwarderror):
        maxFE = np.absolute(forwarderror)

if errormagcoef <= np.absolute(forwarderror/backwarderror):
    errormagcoef = np.absolute(forwarderror/backwarderror)

print("error magniﬁcation factor")
print(errormagcoef)
print("max distance error")
print(maxFE)