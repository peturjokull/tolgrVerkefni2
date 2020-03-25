import numpy as np
import numpy.linalg as lin
import math




def daemi1(r0,r1,r2,r3,r4,c):
	for i in range(10):
		F = np.array([
		[(r0['x'] - r1['x'])**2 + (r0['y'] - r1['y'])**2 + (r0['z'] - r1['z'])**2 - (c*(r1['d'] - r0['d']))**2],
		[(r0['x'] - r2['x'])**2 + (r0['y'] - r2['y'])**2 + (r0['z'] - r2['z'])**2 - (c*(r2['d'] - r0['d']))**2],
		[(r0['x'] - r3['x'])**2 + (r0['y'] - r3['y'])**2 + (r0['z'] - r3['z'])**2 - (c*(r3['d'] - r0['d']))**2],
		[(r0['x'] - r4['x'])**2 + (r0['y'] - r4['y'])**2 + (r0['z'] - r4['z'])**2 - (c*(r4['d'] - r0['d']))**2]
		])
		
		DF = np.array([
		[(r0['x'] - r1['x'])*2 , (r0['y'] - r1['y'])*2 , (r0['z'] - r1['z'])*2 , (2*c**2)*(r1['d'] - r0['d'])],
		[(r0['x'] - r2['x'])*2 , (r0['y'] - r2['y'])*2 , (r0['z'] - r2['z'])*2 , (2*c**2)*(r2['d'] - r0['d'])],
		[(r0['x'] - r3['x'])*2 , (r0['y'] - r3['y'])*2 , (r0['z'] - r3['z'])*2 , (2*c**2)*(r3['d'] - r0['d'])],
		[(r0['x'] - r4['x'])*2 , (r0['y'] - r4['y'])*2 , (r0['z'] - r4['z'])*2 , (2*c**2)*(r4['d'] - r0['d'])]
		])
		print(DF)
		print("nytt")
	
		
		res = lin.solve(DF,F)

		r0['x'] = r0['x'] - res[0][0]
		r0['y'] = r0['y'] - res[1][0]
		r0['z'] = r0['z'] - res[2][0]
		r0['d'] = r0['d'] - res[3][0]
	
	return r0


def daemi2(r0,r1,r2,r3,r4,c):
	ux = [2*(r2['x']-r1['x']),2*(r3['x']-r1['x']),2*(r4['x']-r1['x'])]
	uy = [2*(r2['y']-r1['y']),2*(r3['y']-r1['y']),2*(r4['y']-r1['y'])]
	uz = [2*(r2['z']-r1['z']),2*(r3['z']-r1['z']),2*(r4['z']-r1['z'])]
	ud = [2*(c**2)*(r1['d']-r2['d']),2*(c**2)*(r1['d']-r3['d']),2*(c**2)*(r1['d']-r4['d'])]

	w = [
	(r1['x']**2 - r2['x']**2)+(r1['y']**2 - r2['y']**2)+(r1['z']**2 - r2['z']**2) - (c**2)*(r1['d']**2 - r2['d']**2),
	(r1['x']**2 - r3['x']**2)+(r1['y']**2 - r3['y']**2)+(r1['z']**2 - r3['z']**2) - (c**2)*(r1['d']**2 - r3['d']**2),
	(r1['x']**2 - r4['x']**2)+(r1['y']**2 - r4['y']**2)+(r1['z']**2 - r4['z']**2) - (c**2)*(r1['d']**2 - r4['d']**2)
	]
	
	s1 = - (lin.det([uy,uz,ud])/lin.det([uy,uz,ux]))
	s2 = (lin.det([uy,uz,w])/lin.det([uy,uz,ux]))
	s3 = - (lin.det([ux,uz,ud])/lin.det([ux,uz,uy]))
	s4 = (lin.det([ux,uz,w])/lin.det([ux,uz,uy]))
	s5 = - (lin.det([ux,uy,ud])/lin.det([ux,uy,uz]))
	s6 = (lin.det([ux,uy,w])/lin.det([ux,uy,uz]))

	p1 = (s1**2 + s3**2 + s5**2 - c**2)
	p2 =  2*((c**2)*r1['d'] - s1*(s2+r1['x']) - s3*(s4+r1['y']) - s5*(s6+r1['z']))
	p3 = (s2 + r1['x'])**2 + (s4 + r1['y'])**2 + (s6 + r1['z'])**2 - (c*r1['d'])**2

	d = np.roots([p1,p2,p3])

	x = d*s1 - s2
	y = d*s3 - s4
	z = d*s5 - s6

	svar1 = {'x':x[0],'y':y[0],'z':z[0],'d':d[0]}
	svar2 = {'x':x[1],'y':y[1],'z':z[1],'d':d[1]}

	return svar1,svar2

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
	t['a']=r0['d'] + R[0]/c

	A['b']=rho*math.cos(phi['b'])*math.cos(theta['b'])
	B['b']=rho*math.cos(phi['b'])*math.sin(theta['b'])
	C['b']=rho*math.sin(phi['b'])
	R[1]=math.sqrt(A['b']**2 + B['b']**2 + (C['b']-r0['z'])**2)
	t['b']=r0['d'] + R[1]/c

	A['c']=rho*math.cos(phi['c'])*math.cos(theta['c'])
	B['c']=rho*math.cos(phi['c'])*math.sin(theta['c'])
	C['c']=rho*math.sin(phi['c'])
	R[2]=math.sqrt(A['c']**2 + B['c']**2 + (C['c']-r0['z'])**2)
	t['c']=r0['d'] + R[2]/c
	
	A['d']=rho*math.cos(phi['d'])*math.cos(theta['d'])
	B['d']=rho*math.cos(phi['d'])*math.sin(theta['d'])
	C['d']=rho*math.sin(phi['d'])
	R[3]=math.sqrt(A['d']**2 + B['d']**2 + (C['d']-r0['z'])**2)
	t['d']=r0['d'] + R[3]/c

	dt=10**-8 #nákvæmni klukku í gervihnöttum
	errorcoef=[[0,0,0,0],[0,0,0,1],[0,0,1,0],[0,0,1,1],[0,1,0,0],[0,1,0,1],[0,1,1,0],[0,1,1,1],[1,0,0,0],[1,0,0,1],[1,0,1,0],[1,0,1,1],[1,1,0,0],[1,1,0,1],[1,1,1,0],[1,1,1,1]]
	errorcoef=np.array(errorcoef)
	errormagcoef = 0
	maxFE = 0
	

	for i in range(0,16):
	
		t_new =[
			t['a'] + dt*errorcoef[i,0], 
			t['b'] + dt*errorcoef[i,1],
			t['c'] + dt*errorcoef[i,2],
			t['d'] + dt*errorcoef[i,3]
			]

		for j in range(0,10):
			DF = np.array([
			[(r0['x'] - A['a'])*2 , (r0['y'] - B['a'])*2 , (r0['z'] - C['a'])*2 , (2*c**2)*(t_new[0]- r0['d'])],
			[(r0['x'] - A['b'])*2 , (r0['y'] - B['b'])*2 , (r0['z'] - C['b'])*2 , (2*c**2)*(t_new[1] - r0['d'])],
			[(r0['x'] - A['c'])*2 , (r0['y'] - B['c'])*2 , (r0['z'] - C['c'])*2 , (2*c**2)*(t_new[2]- r0['d'])],
			[(r0['x'] - A['d'])*2 , (r0['y'] - B['d'])*2 , (r0['z'] - C['d'])*2 , (2*c**2)*(t_new[3] - r0['d'])]
			])
			
			F = np.array([
			[(r0['x'] - A['a'])**2 + (r0['y'] - B['a'])**2 + (r0['z'] -  C['a'])**2 - (c**2)*(t_new[0] - r0['d'])**2],
			[(r0['x'] - A['b'])**2 + (r0['y'] - B['b'])**2 + (r0['z'] -  C['b'])**2 - (c**2)*(t_new[1] - r0['d'])**2],
			[(r0['x'] - A['c'])**2 + (r0['y'] - B['c'])**2 + (r0['z'] -  C['c'])**2 - (c**2)*(t_new[2] - r0['d'])**2],
			[(r0['x'] - A['d'])**2 + (r0['y'] - B['d'])**2 + (r0['z'] -  C['d'])**2 - (c**2)*(t_new[3] - r0['d'])**2]
			])
			
			
		res=lin.solve(DF,F)
			
		r0['x'] = r0['x'] - res[0][0]
		r0['y'] = r0['y'] - res[1][0]
		r0['z'] = r0['z'] - res[2][0]
		r0['d'] = r0['d'] - res[3][0]
			
	forwarderror=lin.norm([r0['x'],r0['y'],r0['z']-6370],np.inf)
	backwarderror=c*lin.norm([t_new[0]-t['a'],t_new[1]-t['b'],t_new[2]-t['c'],t_new[3]-t['d']],np.inf)

	if maxFE <= abs(forwarderror):
		maxFE = abs(forwarderror)

	if errormagcoef <= abs(forwarderror/backwarderror):
		errormagcoef = abs(forwarderror/backwarderror)
	
	return errormagcoef,maxFE	

r1 = {'x':15600, 'y':7540 ,'z':20140,'d':0.07074 } #Gervitungl 1
r2 = {'x':18760, 'y':2750 ,'z':18610,'d':0.07220 } #Gervitungl 2
r3 = {'x':17610, 'y':14630,'z':13480,'d':0.07690 } #Gervitungl 3
r4 = {'x':19170, 'y':610  ,'z':18390,'d':0.07242 } #Gervitungl 4
r0 = {'x':0,'y':0,'z':6370,'d':0} #upphafsvigur
r02 = {'x':0,'y':0,'z':6370,'d':0.0001} #upphafsvigur2
phi = {'a':0, 'b': math.pi/8,'c': math.pi/4, 'd': math.pi/4}
theta = {'a':0, 'b': math.pi/2,'c': math.pi, 'd': 3*math.pi/2}
rho = 26570 #km
c=299792.458 #ljóshraði km/s
#svar1,svar2 = daemi1(r0,r1,r2,r3,r4,c)
svar1,svar2 = daemi4(r02,r1,r2,r3,r4,c,rho,theta,phi)

print(svar1)
print(svar2)


#daemi1(r0,r1,r2,r3,r4,c)