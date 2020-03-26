from sympy import symbols , Eq , solve
import numpy as np

c = 299792.458  # Ljoshradi km/s
r1 = {'x':15600, 'y':7540 ,'z':20140,'d':0.07074 } #Gervitungl 1
r2 = {'x':18760, 'y':2750 ,'z':18610,'d':0.07220 } #Gervitungl 2
r3 = {'x':17610, 'y':14630,'z':13480,'d':0.07690 } #Gervitungl 3
r4 = {'x':19170, 'y':610  ,'z':18390,'d':0.07242 } #Gervitungl 4
r0 = {'x':0,'y':0,'z':6370,'d':0} #upphafsvigur

x,y,z,d = symbols('x y z d')

eq1 = Eq((x-r1['x'])**2 + (y-r1['y'])**2 + (z-r1['z'])**2 - (c*(r1['d'] - d))**2, 0)
eq2 = Eq((x-r2['x'])**2 + (y-r2['y'])**2 + (z-r2['z'])**2 - (c*(r2['d'] - d))**2, 0)
eq3 = Eq((x-r3['x'])**2 + (y-r3['y'])**2 + (z-r3['z'])**2 - (c*(r3['d'] - d))**2, 0)
eq4 = Eq((x-r4['x'])**2 + (y-r4['y'])**2 + (z-r4['z'])**2 - (c*(r4['d'] - d))**2, 0)

sol = solve((eq1,eq2,eq3,eq4),(x,y,z,d))

for s in sol:
	print(s)