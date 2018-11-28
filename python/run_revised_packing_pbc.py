#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 13:17:25 2018

@author: chaztikov
"""



"""
Created on Wed Nov 28 10:19:22 2018

@author: chaztikov
"""

import os,sys
cwd=os.getcwd()
sys.path.append(cwd)

import periodic_kdtree
from periodic_kdtree.periodic_kdtree import PeriodicKDTree, PeriodicCKDTree
import numpy as np
import scipy
import scipy.stats
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
# import run_test_packing
# from run_test_packing import *



def porosity_from_radii(radii_scaled, domain_volume,ndimensions = 3,verbose=1):
    pvolumes = radii_scaled**3 * 4 * np.pi / 3 if ndimensions==3 else radii_scaled**2 * np.pi
    porosity = (domain_volume - pvolumes.sum()) / domain_volume
    return porosity,pvolumes

def compute_current_domain_bounding_box(radii_scaled,pts,bbox_vertices,verbose=1):
    radii_scaled.max()
    domain_volume= (xmax-xmin) * (ymax-ymin) * (zmax-zmin)
    porosity,pvolumes = porosity_from_radii(radii_scaled, domain_volume)
    return bbox_vertices, domain_volume, porosity,pvolumes

def compute_correction(x,r,xn,rn,verbose=1):
    if(verbose ==2): print("compute_correction")
    if(verbose ==2): print(xn.shape, x.shape)
    if(verbose ==2): print(x.shape,rn.shape)
    if( xn.shape[0]>1):
        dxn = (xn - x[None,:] ).copy()
        dn = np.linalg.norm( dxn , 2 , axis=1)
        if(verbose ==2): print(dxn.shape,dn.shape)
        if(verbose ==2): print(r)
    else:
        dxn = (xn - x[:] ).copy()
        dn = np.linalg.norm( dxn , 2 , axis=1)
        if(verbose ==2): print(r.shape)
        if(verbose ==2): print(dxn.shape,dn.shape, (r+ rn).shape)

    correct = (rn).copy()

    if(verbose ==2): print(correct.shape)

    xnew = dxn.copy()
    xnew = xnew / dn[:,None]
    if(verbose ==2): print(xnew.shape)
    for i in range(3):
        xnew[:,i] = correct -  dxn[:,i] / dn
    if(verbose ==2): print(dxn.shape)

    xnew = np.mod(xnew,1)

    if(verbose ==2): print("xnew.shape" , xnew.shape)

    return xnew

def overlap_correction(x, r, pts, radii_scaled, RandomState, bbox_vertices,pnorm=2,periodic=1,deltar_min=1.001,verbose=1):
    reinit_flag=0

    if(verbose ==2): print(x.shape)
    if(verbose ==2): print("\n pts.shape \n", pts.shape, radii_scaled.shape)
    # assert(pts.shape[0], radii_scaled.shape[0])

    if(pts.shape[0]<2):
        return x, 0

    if(periodic):
        kdt=PeriodicCKDTree(np.array(periodic_bounds), pts)
        icollide = kdt.query_ball_point(x, deltar_min * np.max(r+radii_scaled))
    else:
        if(verbose ==2): print('aperiodic')
        icollide = np.where(((x-pts[:])**2).sum(axis=1) < (radii_scaled[:] + r) * deltar_min  )

    assert(pts.shape[0], radii_scaled.shape[0])


    if(verbose ==2): print(len(icollide))
    if(verbose ==2): print(icollide)



    if(verbose ==2): print("\n icollide \n", icollide , "\n")
    if( len(icollide) > 1 and type(icollide[0])==int ):
        if(verbose ==2): print(pts.shape, radii_scaled.shape)

        reinit_flag = 1
        periodic = 1

        if(periodic):
            if(verbose ==2): print( '\n periodic \n')

            if(verbose ==2): print(x.shape, r, pts.shape, radii_scaled.shape)
            x = compute_correction(x,r,pts[:],radii_scaled[:])
        else:
            x = compute_correction(x,r,pts[:],radii_scaled[:])

        xnew,reinit_flag = overlap_correction(x, r, pts, radii_scaled, RandomState, bbox_vertices, pnorm,periodic)



    return x, reinit_flag




seed=0
xmax,xmin=1,0
ymax,ymin=1,0
zmax,zmin=1,0
ndimensions=3
verbose=1
radii_dist='lognormal'
radii_dist='uniform'
radius_mu=0.001
radius_sig2=1e-10
tol_self_collision=1e-2

nsamples=int(1e3 * 8)
increment_print = int(nsamples/10)

percentilemin,percentilemax=10,90

target_porosity=0.40

do_rescale=1
periodic_geometry=True
kdt=None

'''
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
'''




domain_volume= (xmax-xmin) * (ymax-ymin) * (zmax-zmin)
coord_min,coord_max,coord_minmax=np.min((xmin,ymin,zmin)), np.max((xmax,ymax,zmax)),np.min((xmax,ymax,zmax))
RandomState = np.random.RandomState(0)
if(verbose ==2): print(radius_mu)

if(radii_dist=='lognormal'):
    lognormal_sig2 = np.log( radius_sig2 / radius_mu**2 + 1 )
    lognormal_mu = np.log( radius_mu ) - lognormal_sig2 / 2
    radii = RandomState.lognormal( lognormal_mu, lognormal_sig2 ,nsamples)
    radii.mean(),np.std(radii)

    '''truncate sampled distribution'''
    pmin = np.percentile(radii,percentilemin),
    pmax = np.percentile(radii,percentilemax)
    radii = radii[radii<pmax]
    radii = radii[radii>pmin]
    nsamplesclip = nsamples - radii.shape[0]
    nsamples = radii.shape[0]


else:
    radii = np.ones([nsamples]) * radius_mu



Z=None
pts=None
rmax=1
dmax=1
rscale=1

'''rescale radii to obtain desired porosity'''
rscale = ( (  np.sum( np.pi * (4./3.) * radii**3) ) / ( domain_volume * (1 - 0.9 * target_porosity) ) )**(1/3)
Z = radii/rscale
radius_mu_scaled = radius_mu / rscale
radius_sig2_scaled = radius_sig2 / rscale**2
rmax = ( Z ).max()
dmax = 2*rmax


# '''plot original and transformed distributions'''
# plt.hist(Z);
# plt.title( "Transformed Lognormal Parameters: \n" + str( np.round( scipy.stats.lognorm.fit( np.atleast_2d(Z) ), 3) ))
# plt.show()


# plt.hist(radii);
# plt.title( "Original Lognormal Parameters: \n" + str( np.round( scipy.stats.lognorm.fit( np.atleast_2d(radii) ), 3) ))
# plt.show()



bbox_vertices = [(xmin,ymin,zmin),(xmax,ymax,zmax)]
periodic_bounds=bbox_vertices[1]



''' Detect Collisions and Translate Spheres '''
registered = []
unregistered = [i for i in range(nsamples)]
boundary = []
t_list = []
tlast = time.time()
radii_scaled = np.array([])
pts=np.array([])
rands=[]
tlast = time.time()
for i,radius in enumerate(Z):

    r = radius

    pt=np.zeros([3])
    pt[0] = RandomState.uniform(0,xmax)
    pt[1] = RandomState.uniform(0,ymax)
    pt[2] = RandomState.uniform(0,zmax)

    x = pt.copy()
    x = np.atleast_2d(x)


    if(i>0 and i%increment_print==0):
        bbox_vertices, domain_volume, current_porosity, pvolumes = compute_current_domain_bounding_box(radii_scaled, pts, bbox_vertices)
        if(verbose ==1): print("\n Iteration",i,"\n Curr. Sph. Posn. ",x, "\n Sph. Rad.",r,
              "\n iter time", time.time()-tlast, "\n current porosity" , current_porosity,
              "\n bbox_vertices", bbox_vertices, "\n domain_volume" , domain_volume)


    else:
        assert( np.all( np.isfinite(pts)) )
        
        x,reinit_flag = overlap_correction(x, r, pts, radii_scaled, RandomState, bbox_vertices,2,1)

        if(reinit_flag==1):
            if(periodic_geometry==1):
                break;

        elif(reinit_flag==0):


            registered.append(i)

            if(i>0):
                pts = np.vstack((pts,x))
                radii_scaled = np.hstack((radii_scaled,r))
            else:
                pts = np.array([pt.tolist()])
                radii_scaled = np.array([radius])

            if(verbose ==2): print(pts.shape, radii_scaled.shape)

            try:
                unregistered.remove(i)
            except Exception as e:
                print("\n ",e)
                break;
    t_list.append(time.time() - tlast)
    tlast = time.time()





registered=np.array(registered)
domain_volume=xmax*ymax*zmax
bbox_vertices, domain_volume, porosity, pvolumes = compute_current_domain_bounding_box(radii_scaled, pts, bbox_vertices)
save_filename = 'packing'
#save packing output
idx_points = np.arange(0,len(registered)+0)






print("\n domain size ", " [xmin,xmax] ", xmin,xmax, " [ymin,ymax] ", ymin,ymax, " [zmin,zmax] ", zmin,zmax, " volume " , domain_volume)
print('particle volumes: sum, mean, median, min, max', pvolumes.sum(), pvolumes.mean(), np.median(pvolumes) , pvolumes.min(), pvolumes.max())
print("\n \n \n porosity ", porosity)
print("\n number of spheres ", registered.shape)
print("\n number of registered spheres ", registered.shape)
print("\n number of unregistered spheres ", registered.shape )
print("\n sphere distribution parameters ", radius_mu, radius_sig2)
print("\n mean coordination number ", )
print("\n \n \n ")

stacked_data = np.vstack((idx_points.astype(int), pts[:,0], pts[:,1], pts[:,2], radii_scaled[:])).T
np.savetxt(save_filename + ".txt", stacked_data, header="ID x y z r", fmt='%i,%E,%E,%E,%E')
'''
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
'''





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 13:13:16 2018

@author: chaztikov
"""
#import revised_packing_pbc.py
import os,sys
cwd=os.getcwd()
sys.path.append(cwd)
import numpy as np
import scipy
import scipy.stats

import matplotlib as mpl
import matplotlib.pyplot as plt
# import run_test_packing
# from run_test_packing import *

import numpy as np

idx,x,y,z,r = np.loadtxt("packing.txt",delimiter=",",skiprows=1,unpack=True)

# Check Porosity
volume = 0.
for i in range(0,len(idx)):
	volume = volume + 4./3.*np.pi*r[i]*r[i]*r[i]
print("check porosity " , 1.-volume)


# # Check Distribution Parameters
# print(np.mean(r))
# print(np.std(r))

#Check Overlap
import sklearn.metrics
from sklearn.metrics import pairwise_distances as pwdist
pts = np.vstack((x,y,z)).T
dist=pwdist(pts)
row,col=np.where(dist<r[:,None]+r[None,:])
[dist[ri,ci] for ri,ci in zip(row,col) if ri<ci]


# print("Total Number of Overlaps ", count-1)



from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


psize = (pvolumes-pvolumes.min() )
# psize = psize/psize.max()


import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import numpy as np

try:
    trace = go.Scatter3d(x,y,z,
        mode = 'markers',
        marker = dict(
            color = '#FFBAD2',
            line = dict(width = 1)
        )
    )
    data = [trace]
    py.iplot(data, filename='packing_revised.txt',fileopt='new')
except Exception as e:
    print("\n \n",e)
