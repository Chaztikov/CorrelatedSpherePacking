# Sphere Pack Analysis 

import numpy as np

x,y,z,r = np.loadtxt("pack.out",unpack=True)

# Check Porosity
volume = 0.
for i in range(0,len(x)):
	volume = volume + 4./3.*np.pi*r[i]*r[i]*r[i]
print(volume)


# Check Distribution Parameters
print(np.mean(r))
print(np.std(r))

#Check Overlap
count = 0 
for i in range(0,len(x)):
	for j in range(i+1,len(x)):
		distance = np.sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) )
		if (distance < (r[i]+r[j])):
			print("OVERLAP!!!!"); print(i,j,x[i],y[i],z[i],x[j],y[j],z[j],distance-(r[i]+r[j]))
			count  = count + 1; print(count)