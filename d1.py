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
	
		
		res = lin.solve(DF,F)
		print(res)
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



r1 = {'x':15600, 'y':7540 ,'z':20140,'d':0.07074 } #Gervitungl 1
r2 = {'x':18760, 'y':2750 ,'z':18610,'d':0.07220 } #Gervitungl 2
r3 = {'x':17610, 'y':14630,'z':13480,'d':0.07690 } #Gervitungl 3
r4 = {'x':19170, 'y':610  ,'z':18390,'d':0.07242 } #Gervitungl 4
r0 = {'x':0,'y':0,'z':6370,'d':0} #upphafsvigur
c=299792.458 #ljóshraði km/s
svar1,svar2 = daemi(r0,r1,r2,r3,r4,c)

print(svar1)
print(svar2)


#daemi1(r0,r1,r2,r3,r4,c)