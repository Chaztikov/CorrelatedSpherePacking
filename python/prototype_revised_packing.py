# import sys, os
# cwd = os.getcwd()
# addedpath=cwd+"/periodic_kdtree/"
# print("Adding path: ", addedpath)
# sys.path.append(addedpath)
# print("New path including required packages: ", sys.path)


# import numpy as np
# import periodic_kdtree
# from periodic_kdtree import PeriodicCKDTree
# import scipy.spatial
# import time
# from mpl_toolkits.mplot3d import Axes3D
# import pickle
# import gzip
# import networkx as netx
# import getopt

# import matplotlib.pyplot as plt
# import matplotlib as mpl
# plt.style.use('ggplot')
# # %matplotlib inline
# mpl.rcParams['image.cmap'] = 'autumn'

import os,sys
cwd=os.getcwd()
sys.path.append(cwd)
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import run_test_packing
from run_test_packing import *
seed=0
xmax,xmin=1,0
ymax,ymin=2,0
zmax,zmin=4,0
ndimensions=3
radii_dist='lognormal'
radius_mu=1
radius_sig2=0.25
nsamples=1000
percentilemin,percentilemax=5,95

target_porosity=0.4



def compute_correction(x,r,ineighbors,pts,radii_scaled):
    print("COMPUTE CORRECTION")
    xn = pts
    rn = radii_scaled

    dxn = (xn - x[None,:]) 
    dn = np.linalg.norm(x[None,:]-xn,2,axis=1)
    correct = (r+ rn)
    correct = correct[:,None] - dxn / dn[:,None]
    correct = correct * dn[:,None]
    for c in correct:
        x += c
    return x

def overlap_correction(i, x, r, pts, radii_scaled, registered, pnorm=2):
#     pts = allpts[registered]
#     radii_scaled=allradii_scaled[registered]
#     search_radius = search_radius_factor_of_max_diameter * dmax    

    print(pts.shape)
    dist_collide = np.linalg.norm(x[None,:] - pts, pnorm, axis=1)
    my_collision_length = radii_scaled + r
    icollide = (dist_collide < my_collision_length)
    icollidesum = (icollide).sum()
        
    if( icollidesum < 1 ):
        reinit_flag = 0

    else:
        icollide = np.where(icollide)
        x = compute_correction(x,r,icollide,pts,radii_scaled)
        reinit_flag=0

    return x, reinit_flag

'''
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
'''

coord_min,coord_max,coord_minmax=np.min((xmin,ymin,zmin)), np.max((xmax,ymax,zmax)),np.min((xmax,ymax,zmax))
RandomState = np.random.RandomState(0)
print(radius_mu)
lognormal_sig2 = np.log( radius_sig2 / radius_mu**2 + 1 )
lognormal_mu = np.log( radius_mu ) - lognormal_sig2 / 2
Z = RandomState.lognormal( lognormal_mu, lognormal_sig2 ,nsamples)
Z.mean(),np.std(Z)
# RandomState.lognormal?
# nsamples=int(1e5)

# kdt = scipy.spatial.KDTree(pts) #,leafsize=set_leafsize_factor * num_neighbors )
# x,reinit_flag = overlap_correction(i, x, r, pts, radii_scaled, kdt, registered, unregistered, dmax, search_radius_factor_of_max_diameter, pnorm=2, eps=kdt_eps)

    
approx_samples_per_dim = int((nsamples/10)**(1/3))
bbox_rescale_factor = 2*Z.max() * approx_samples_per_dim
xmax *= bbox_rescale_factor
ymax *= bbox_rescale_factor
zmax *= bbox_rescale_factor

''' Detect Collisions and Translate Spheres '''
registered = []
unregistered = [i for i in range(nsamples)]
boundary = []
t_list = []
tlast = time.time()
for i,r in enumerate(Z):
#     pt = RandomState.uniform(0,1)
    pt = RandomState.uniform(0,1,[ndimensions])
    
    pt[0] = (pt[0] ) * (xmax-xmin) + xmin
    pt[1] = (pt[1] ) * (ymax-ymin) + ymin
    pt[2] = (pt[2] ) * (zmax-zmin) + zmin
    x = pt.copy()
    if(i%1000==0):
        print(i,x,r, time.time())
    if(i==0):
        registered.append(i)
        pts = np.array([x.tolist()])
        radii_scaled = np.array([r])
    else:
        x,reinit_flag = overlap_correction(i, x, r, pts, radii_scaled,registered)
        if(reinit_flag==1):
            if(periodic_geometry==1):
                break;
        elif(reinit_flag==0):                    
            registered.append(i)
#             pts[i] = x
            pts = np.vstack((pts,pt))
            radii_scaled = np.hstack((radii_scaled,r))
            Z[i] = r
    unregistered.remove(i)

    t_list.append(time.time() - tlast)
    tlast = time.time()
    
registered=np.array(registered)
domain_volume=xmax*ymax*zmax
save_filename = 'packing'
#save packing output
idx_points = np.arange(0,len(registered))

stacked_data = np.vstack((idx_points.astype(int), pts[:,0], pts[:,1], pts[:,2], radii_scaled[:])).T
np.savetxt(save_filename + ".txt", stacked_data, header="ID x y z r", fmt='%i,%E,%E,%E,%E')




pvolumes = radii_scaled**3 * 4 * np.pi / 3 if ndimensions==3 else radii_scaled**2 * np.pi 
porosity = 1 - (domain_volume - pvolumes.sum()) / domain_volume

print("\n domain size ", " [xmin,xmax] ", xmin,xmax, " [ymin,ymax] ", ymin,ymax, " [zmin,zmax] ", zmin,zmax, " volume " , domain_volume)
print('particle volumes: sum, mean, median, min, max', pvolumes.sum(), pvolumes.mean(), np.median(pvolumes) , pvolumes.min(), pvolumes.max())
print("\n \n \n porosity ", porosity)
print("\n number of spheres ", registered.shape)
print("\n number of registered spheres ", registered.shape)
print("\n number of unregistered spheres ", registered.shape )
print("\n sphere distribution parameters ", radius_mu, radius_sig2)
print("\n mean coordination number ", )
print("\n \n \n ")

'''
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
'''

import numpy as np

idx,x,y,z,r = np.loadtxt("packing.txt",delimiter=",",skiprows=1,unpack=True)

# Check Porosity
volume = 0.
for i in range(0,len(idx)):
	volume = volume + 4./3.*np.pi*r[i]*r[i]*r[i]
print(1.-volume)


# Check Distribution Parameters
print(np.mean(r))
print(np.std(r))

#Check Overlap
count = 0 
for i in range(0,len(idx)):
    add_to_count=0;
    for j in range(i+1,len(idx)):
        distance = np.sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) )
        if (distance < (r[i]+r[j])):
            add_to_count=1
    count  = count + add_to_count; 
    
print("Total Number of Overlaps ", count-1)= np.loadtxt("packing.txt",delimiter=",",skiprows=1,unpack=True)

# Check Porosity
volume = 0.
for i in range(0,len(idx)):
	volume = volume + 4./3.*np.pi*r[i]*r[i]*r[i]
print(1.-volume)


# Check Distribution Parameters
print(np.mean(r))
print(np.std(r))

#Check Overlap
count = 0 
for i in range(0,len(idx)):
    add_to_count=0;
    for j in range(i+1,len(idx)):
        distance = np.sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) )
        if (distance < (r[i]+r[j])):
            add_to_count=1
    count  = count + add_to_count; 
    
print("Total Number of Overlaps ", count-1)