import math
import numpy as np
import numpy.linalg as lin

def daemi4(r0,r1,r2,r3,r4,c,rho,theta,phi):

	A= {'a' : r1['x'], 'b' : r2['x'], 'c' : r3['x'], 'd' : r4['x']}
	B= {'a' : r1['y'], 'b' : r2['y'], 'c' : r3['y'], 'd' : r4['y']}
	C= {'a' : r1['z'], 'b' : r2['z'], 'c' : r3['z'], 'd' : r4['z']}
	t= {'a' : r1['d'], 'b' : r2['d'], 'c' : r3['d'], 'd' : r4['d']}
	R=np.zeros((4,), dtype=float)



	A['a']=rho*math.cos(phi['a'])*math.cos(theta['a'])
	B['a']=rho*math.cos(phi['a'])*math.sin(theta['a'])
	C['a']=rho*math.sin(phi['a'])
	R[0]=math.sqrt(A['a']**2 + B['a']**2 + (C['a']-r0['z'])**2)
	t['a']=r0['z'] + R[0]/c

	A['b']=rho*math.cos(phi['b'])*math.cos(theta['b'])
	B['b']=rho*math.cos(phi['b'])*math.sin(theta['b'])
	C['b']=rho*math.sin(phi['b'])
	R[1]=math.sqrt(A['b']**2 + B['b']**2 + (C['b']-r0['z'])**2)
	t['b']=r0['z'] + R[1]/c

	A['c']=rho*math.cos(phi['c'])*math.cos(theta['c'])
	B['c']=rho*math.cos(phi['c'])*math.sin(theta['c'])
	C['c']=rho*math.sin(phi['c'])
	R[2]=math.sqrt(A['c']**2 + B['c']**2 + (C['c']-r0['z'])**2)
	t['c']=r0['z'] + R[2]/c

	A['d']=rho*math.cos(phi['d'])*math.cos(theta['d'])
	B['d']=rho*math.cos(phi['d'])*math.sin(theta['d'])
	C['d']=rho*math.sin(phi['d'])
	R[3]=math.sqrt(A['d']**2 + B['d']**2 + (C['d']-r0['z'])**2)
	t['d']=r0['z'] + R[3]/c

	dt=10**-8 #nákvæmni klukku í gervihnöttum
	errorcoef=[[0,0,0,0],[0,0,0,1],[0,0,1,0],[0,0,1,1],[0,1,0,0],[0,1,0,1],[0,1,1,0],[0,1,1,1],[1,0,0,0],[1,0,0,1],[1,0,1,0],[1,0,1,1],[1,1,0,0],[1,1,0,1],[1,1,1,0],[1,1,1,1]]
	errorcoef=np.array(errorcoef)
	errormagcoef = 0
	maxFE = 0

	for i in range(0,15):
		#t_new = t + np.multiply(dt,errorcoef[i])
		t_new = t + dt*errorcoef[i]
		for j in range(0,9):
			DF = np.array([
			[(r0[0] - A[0])*2 , (r0[1] - B[0])*2 , (r0[2] - C[0])*2 , (2*c**2)*(t_new[0]- r0[3])],
			[(r0[0] - A[1])*2 , (r0[1] - B[1])*2 , (r0[2] - C[1])*2 , (2*c**2)*(t_new[1] - r0[3])],
			[(r0[0] - A[2])*2 , (r0[1] - B[2])*2 , (r0[2] - C[2])*2 , (2*c**2)*(t_new[2]- r0[3])],
			[(r0[0] - A[3])*2 , (r0[1] - B[3])*2 , (r0[2] - C[3])*2 , (2*c**2)*(t_new[3] - r0[3])]
			])

			F = np.array([
			[(r0[0] - A[0])**2 + (r0[1] - B[0])**2 + (r0[2] -  C[0])**2 - (c*(t_new[0] - r0[3]))**2],
			[(r0[0] - A[1])**2 + (r0[1] - B[1])**2 + (r0[2] -  C[1])**2 - (c*(t_new[0] - r0[3]))**2],
			[(r0[0] - A[2])**2 + (r0[1] - B[2])**2 + (r0[2] -  C[2])**2 - (c*(t_new[0] - r0[3]))**2],
			[(r0[0] - A[3])**2 + (r0[1] - B[3])**2 + (r0[2] -  C[3])**2 - (c*(t_new[0] - r0[3]))**2]
			])

			res=lin.solve(DF,F)
			res = lin.solve(DF,F)
		r0['x'] = r0['x'] - res[0][0]
		r0['y'] = r0['y'] - res[1][0]
		r0['z'] = r0['z'] - res[2][0]
		r0['d'] = r0['d'] - res[3][0]

		forwarderror=lin.norm([r0['x'],r0['y'],r0['z']-6370],np.inf)
		backwarderror=c*lin.norm(t_new-t,np.inf)

		if maxFE <= abs(forwarderror):
			maxFE = abs(forwarderror)

		if errormagcoef <= abs(forwarderror/backwarderror):
			errorcoeff = abs(forwarderror/backwarderror)
	#return 1,2
	return errorcoeff,maxFE	

r1 = {'x':15600, 'y':7540 ,'z':20140,'d':0.07074 } #Gervitungl 1
r2 = {'x':18760, 'y':2750 ,'z':18610,'d':0.07220 } #Gervitungl 2
r3 = {'x':17610, 'y':14630,'z':13480,'d':0.07690 } #Gervitungl 3
r4 = {'x':19170, 'y':610  ,'z':18390,'d':0.07242 } #Gervitungl 4
r0 = {'x':0,'y':0,'z':6370,'d':0} #upphafsvigur

err, maxEf = daemi4(r0,r1,r2,r3,r4,c,rho,0.24,0.24)
