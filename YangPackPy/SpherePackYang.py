# import gpflow
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
%matplotlib inline
import matplotlib as mpl
mpl.rcParams['image.cmap'] = 'autumn'
import periodic_kdtree
from periodic_kdtree import PeriodicCKDTree

import scipy.spatial
import time
ndimensions=3
xmin,xmax = 0,1
num_neighbors = 27
set_leafsize = num_neighbors * 2 
periodic_geometry = True
kdt_eps = 1e-2
pnorm = 2

radii_dist='lognormal'
radius_mu = 1
radius_sig2 = .05
nbins=100
nsamples = int(1e4)
show_hist=False
target_porosity = 0.44
percentilemin = 5
percentilemax = 95
zmax = ymax = xmax
find_all_neighbors = False

#forcing params
set_dt = 0.01
blfac = 0.01
# set_fac = pvol[:,None,None]
damping=1
tstep=0
tprint=100
tstepmax = 1000
force_abstol = 1e-7
fmin = 1e6
# fint = vel.copy() + 2*force_abstol
fabsmax = 2 * force_abstol
tstep=0
tprint=100
tstepmax = 10000
force_abstol = 1e-7
force_reltol = 1e-8
fabsmax = 2 * force_abstol
fmin = 1e6

Cmu = np.array([1,2,4])
Csig = np.diag([1,2,3])

max_overlap_factor = 1


import networkx as netx
def from_neighbors_to_edges(neighbors,return_type='nparray'):
#     return frozenset((irow,col)  for irow,row in enumerate(neighbors) for col in row if col!=irow)
    edges= ( (irow,col)  for irow,row in enumerate(neighbors) for col in row if col!=irow )
    if return_type=='nparray':
        edges=np.array(list(edges))
        return edges,len(edges)
    if return_type=='iterator':
        return edges,len(list(edges))
    
def get_adjacency_graph(edges, get_matrix=False, get_list=False, compute_laplacian=False):
    G = netx.Graph()
    G.add_edges_from(edges)
    adj_data = netx.json_graph.adjacency_data(G)
    Adj = netx.adjacency_graph(adj_data)
    if(get_matrix==True):
        return netx.to_scipy_sparse_matrix(Adj)
    if(get_list==True):
        return adj_data
    # Deg = netx.de
    if(compute_laplacian==True):
        Deg = np.diag(np.array(G.degree))
        GL = netx.laplacian_matrix(G)
        GLS = netx.laplacian_spectrum(G)
        return G,Adj,Deg,GL,GLS
    else:
        return G,Adj


'''form uniform hexahedral mesh'''
# ncells = nsamples
# np.m
ncells = nsamples
xx=np.linspace(xmin,xmax,ncells+1)
# nodes = np.meshgrid(xx,xx)
# cells = 
dcell = np.diff(xx)[0]
domain_volume = (xmax-xmin)**3



periodic_bounds = np.array([xmax,ymax,zmax])[:ndimensions]



#PART I
'''sample radii'''
# #sample radii 
# while target_porosity < v0:
# v0 = (radii**3 * 4*np.pi/(3 * (xmax-xmin)**3))
# v0 = v0/v0.sum()
if(radii_dist=='lognormal'):
    mu = np.log(radius_mu)
    sig2 = radius_sig2
    Z = np.random.lognormal(mu,sig2,nsamples)
    radii = np.exp(Z)

if(show_hist==True):
    plt.hist(radii,nbins)
    plt.show()

pmin = np.percentile(radii,percentilemin),
pmax = np.percentile(radii,percentilemax)
radii = radii[radii<pmax]
radii = radii[radii>pmin]
nsamplesclip = nsamples - radii.shape[0]
nsamples = radii.shape[0]# - nsamplesclip
radii = np.sort(radii)

print('1,99 percentiles ', pmin,pmax)
print('max, min, mean', radii.max(),radii.min(),radii.mean())



'''rescale radii'''
# rscale = ( 4*np.pi / (3 * target_porosity * domain_volume) * np.sum(radii**3))
rscale = ( (  np.sum(4*np.pi*(1/3.) * radii**3) ) / ( domain_volume * target_porosity ) )**(1/3)
radii_scaled = radii/rscale
radius_mu_scaled = radius_mu / rscale
radius_sig2_scaled = radius_sig2 / rscale**2 #?
rmax = ( radii_scaled ).max()
dmax = 2*rmax
delta = dcell - 2*rmax
search_radius = 2*dmax#+delta


    
'''sample isotropic gaussian correlation'''
Cnorm = np.random.multivariate_normal(Cmu,Csig,nsamples)#,[nsamples,ndimensions])

'''sample points'''
pts = np.random.uniform(xmin,xmax,[nsamples,ndimensions])

'''get nearest neighbor distances'''
t0 = time.time()

if(periodic_geometry==True):
    kdt = PeriodicCKDTree(periodic_bounds, pts)
else:
    kdt = scipy.spatial.KDTree(pts,leafsize=set_leafsize)
if(find_all_neighbors==True):
    dist,neighbors = kdt.query(pts, k=num_neighbors, eps=kdt_eps, p = pnorm)
    
    print("NNE time " , time.time()-t0)

    '''sort by func(distances) eg sum'''
    distsum = (dist[:,1:].sum(axis=1))
    distmean = (dist[:,1:].mean(axis=1))
    distmedian = np.median(dist[:,1:],axis=1)
    distmin = (dist[:,1:].min(axis=1))
    distmax = (dist[:,1:].max(axis=1))

    isort_pts = np.argsort(distmedian)
    isort_radii = np.argsort(isort_pts)

    sorted_radii = radii[isort_radii].copy()



    edges = from_neighbors_to_edges(neighbors)[0]

'''overlap'''
max_overlap = radius_mu_scaled * max_overlap_factor
#find neighbors I interesect (in 3/9/27 cells), use centers and radii to move away by dx if ||dx||<overlap_max
def compute_correction(x,r,ineighbors):
    xn = pts[ineighbors]
    rn = radii_scaled[ineighbors]
    dn = np.linalg.norm(x[:,None]-xn,2,axis=0)
    correct = (r[:,None] + rn) - dn * (xn - x[:,None]) / dn
    correct = np.sum(correct, axis=0)
    return correct
    
def overlap_correction(i, x, r, kdt, registered, unregistered, search_radius=2*dmax, pnorm=2, eps=kdt_eps):       
    ineighbors = kdt.query_ball_point(x, search_radius, pnorm, eps)
    ineighbors = np.array(ineighbors) 
    
    if(ineighbors.shape[0] < 2 ):
        return x

    if(ineighbors.shape[0]>2):
        dist_collide = np.linalg.norm(x[None,:] - pts[ineighbors], pnorm)#, axis=2)
    else:
        dist_collide = np.linalg.norm(x[:]-pts[ineighbors], pnorm)

    my_collision_length = radii_scaled[ineighbors] + r


    icollide = (dist_collide < my_collision_length)
    try:
        if( (icollide).sum() < 1 ):
            return x
        else:
            icollide = np.where(icollide)
            icollneigh = icollide[ineighbors]
            correct = compute_correction(x, r, icollneigh)
            x = x + correct
            x = overlap_correction(i, x, r, pts, periodic_bounds, set_leafsize, registered, unregistered, search_radius=2*dmax, pnorm=2, eps=kdt_eps)
            return x #hopefully...

    except Exception:
        print('Exception? ')
        return x
    
    return correct

def compute_all_overlap(pts, neighbors, nsamples, collision_length_factor):
    dist_collide = np.linalg.norm(pts[neighbors[:,0:1]]-pts[neighbors[:,1:]], axis=2)
    icollide = dist_collide < collision_length_factor
    icollide = np.where(icollide)
    pcoll, ncoll = icollide
    pcoll = np.unique(pcoll)
    ncoll = np.unique(ncoll)

    num_pcoll = np.unique(icollide[0]).shape[0] 
    num_ncoll = np.unique(icollide[1]).shape[0]
    (num_pcoll / nsamples), (num_ncoll / nsamples)
    return icollide, pcoll, ncoll, num_pcoll, num_ncoll

'''BC, EQ separation, pore throat size, collision distances, PD'''
boundary_layer = [blfac, xmax *(1- blfac)]
iboundary = np.any(((pts<boundary_layer[0]).astype(int)+(pts>boundary_layer[1]).astype(int)) , axis=1)
iinterior = np.all((pts>boundary_layer[0]).astype(int)*(pts<boundary_layer[1]).astype(int),axis=1)
print(iboundary.shape,iinterior.shape)
iinterior.sum(),iboundary.sum()

#formerly, scaled_radii was defined another way
scaled_radii = radii_scaled.copy()
scaled_radius_mu = radius_mu_scaled.copy()

num_inclusions_per_dim = int((nsamples)**(1/ndimensions))
if(find_all_neighbors==True):
    eq_length_factor = (scaled_radii[neighbors[:,0:1]] + scaled_radii[neighbors[:,1:]])
else:
    eq_length_factor = (4 * np.pi * (1/3) * (scaled_radii**3).mean())**(1/3) * np.array([[1]]) #???
pore_space_per_particle = (xmax - eq_length_factor.mean()*num_inclusions_per_dim)/num_inclusions_per_dim
medimean_eq_length = np.median( eq_length_factor.mean(axis=1))
porespace_per_dim = num_inclusions_per_dim * medimean_eq_length
porespace_per_particle  = (porespace_per_dim / (num_inclusions_per_dim - 1))/2
scaled_radius_diam = scaled_radius_mu*2
#set spacing
collision_length_factor = eq_length_factor.copy()# - pore_space_per_particle /2
eq_length_factor = eq_length_factor + pore_space_per_particle /2
horizon_factor = eq_length_factor * 1
# nneighbors = neighbors.shape[1]
tsteps = np.arange(0,tstepmax)

zmin = ymin = xmin
zmax = ymax = xmax

(pts.T[-1]-np.mod(pts.T[-1],xmax)).max()
# pts.T[-1][ (pts.T[-1] - radii.T)>xmax]
# (p - radii.T)>xmax
(pts.T[-1]+radii_scaled > xmax).sum()
((pts**2).sum(axis=1)+radii_scaled**2 > xmax**2).sum()

def where_boundary_intersect(pts,radii_scaled):
    cond = (pts**2).sum(axis=1)-radii_scaled**2 > xmax**2
    cond = cond + ((pts**2).sum(axis=1)-radii_scaled**2 < xmin**2)
    cond = cond + ((pts**2).sum(axis=1)+radii_scaled**2 > ymax**2)
    cond = cond +((pts**2).sum(axis=1)-radii_scaled**2 < ymin**2)
    cond = cond + ((pts**2).sum(axis=1)+radii_scaled**2 > zmax**2)
    cond = cond +((pts**2).sum(axis=1)-radii_scaled**2 < zmin**2)
    # cond = cond + (pts**2).sum(axis=1)+radii_scaled**2 > zmax**2
    # cond = cond + (pts**2).sum(axis=1)-radii_scaled**2 < -zmax**2
    return cond, cond.sum()
cond, conds = where_boundary_intersect(pts,radii_scaled)
conds
# (np.linalg.norm(pts,2,axis=1)+radii_scaled > xmax).sum()
# pts.min()

def get_reflected_pts(pts,radii_scaled,xmin,xmax):
    sz0 = pts.shape[0]
    ptsr = pts**2 - radii_scaled[:,None]**2 
    ip1 = np.any(ptsr < xmin,axis=1)
    pts1 = np.mod(ptsr[ip1],xmax-xmin)
    pts = np.vstack((pts,pts1))
    radii_scaled = np.concatenate((radii_scaled,radii_scaled[ip1]))

    ptsr = pts**2 + radii_scaled[:,None]**2
    ip1 = np.any(ptsr > xmax,axis=1)
    pts1 = np.mod(ptsr[ip1],xmax-xmin)
    pts = np.vstack((pts,pts1))
    radii_scaled = np.concatenate((radii_scaled,radii_scaled[ip1]))
    sz = np.unique( pts ,axis=0).shape[0]
    print(sz,sz0)
    done = 1 if( sz==sz0) else 0
        
    return pts, radii_scaled,flag

flag=0
# while flag==0:
pts, radii_scaled,flag = get_reflected_pts(pts,radii_scaled,xmin,xmax)
print(radii_scaled.shape)
nsamples = radii_scaled.shape[0]
assert(radii_scaled.shape[0] == pts.shape[0])

# ip2 = np.any(ptsr > xmax,axis=1)
# pts[ip1]
# pts[ip2]
# pts[ip2]+radii_scaled[ip2]
# ip2.sum()
# ip1.sum()

''' REVERSE, PROBABLY DO THIS EARLIER '''
radii_scaled = radii_scaled[::-1]
pts = pts[::-1,:]

registered = []
unregistered = [i for i in range(nsamples)]
boundary = []
t_list = []
tlast = time.time()
for i,(x,r) in enumerate(zip( pts , radii_scaled)):
    if(i%1000==0):
        print(i,x,r, time.time())
    if(i==0):
        registered.append(i)
        unregistered.remove(i)
        pts[i] = x
        radii_scaled[i] = r
    else:
        x = overlap_correction(i, x, r, kdt, registered, unregistered, search_radius=2*dmax, pnorm=2, eps=kdt_eps)
#         x = overlap_correction(i, x, r, pts, periodic_bounds, set_leafsize, registered, unregistered, 
#                                search_radius=2*dmax, pnorm=2, eps=kdt_eps)
        registered.append(i)
        unregistered.remove(i)
        pts[i] = x
        radii_scaled[i] = r

    t_list.append(time.time() - tlast)
    tlast = time.time()
t_list = np.array(t_list)

import pickle
import gzip
savedict = {'points':pts,
           'radii_scaled':radii_scaled,
           'target_porosity':target_porosity,
            'boundary':boundary
           }
filename = 'uncorr_pack_dim_'+str(ndimensions)+'_nparticles_'+str(nsamples)
with gzip.open(filename + '_.pkl','wb') as f:
    pickle.dump(savedict,f)
    


radii_scaled
registered = np.array(registered)
pvolumes = radii_scaled**3 * 4 * np.pi / 3 if ndimensions==3 else radii_scaled**2 * np.pi 
psize = ((pvolumes-pvolumes.min()))
psize = psize/psize.max()
psize = psize * nsamples**(1/ndimensions) #/16



from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(16,16))
ax = fig.add_subplot(111, projection='3d')
X,Y,Z = pts[:,0],pts[:,1],pts[:,2]
cube = ax.scatter(X, Y, Z, zdir='z', s=psize, c=psize,alpha=0.4)
cbar = fig.colorbar(cube, shrink=0.6, aspect=5)
plt.savefig(filename+'_3d.png',dpi=200)
plt.show()



# len(),nsamples
plt.figure(figsize=(16,16))
plt.scatter(pts[:,0],pts[:,1],s=psize,c=psize,alpha=0.6)
plt.colorbar()
# plt.xlim(0,0.1)
# plt.ylim(0,0.1)
plt.savefig(filename+'_.png',dpi=200)
plt.show()


# plt.figure(figsize=(16,8))
# plt.plot(t_list,'.',alpha=0.4,ms=4);
# plt.ylabel('Time Elapsed Per Iteration ')
# plt.xlabel('Iteration (Sequential Addition)')
# plt.savefig(filename+'_tiame.png',dpi=200)
# plt.show()

