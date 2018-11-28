import sys, os
cwd = os.getcwd()
addedpath=cwd+"/periodic_kdtree/"
print("Adding path: ", addedpath)
sys.path.append(addedpath)
print("New path including required packages: ", sys.path)


import numpy as np
import periodic_kdtree
from periodic_kdtree import PeriodicCKDTree
import scipy.spatial
import time
from mpl_toolkits.mplot3d import Axes3D
import pickle
import gzip
import networkx as netx
import getopt

import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('ggplot')
# %matplotlib inline
mpl.rcParams['image.cmap'] = 'autumn'


def read_input_parameters(param_filename='Parameters.in'):
    input_parameters = np.genfromtxt(param_filename,delimiter='=',comments='#',dtype=str)
    input_parameters[:,0]
    input_parameters[:,1]

    # input_parameter_names, input_parameter_values = np.split( input_parameters, [1], axis=1)
    input_parameter_names = input_parameters[:,0]
    input_parameter_values = input_parameters[:,1]
    parameters = dict({})
    for name, value in zip(input_parameter_names, input_parameter_values):
    #     print(name)
        #assume parameter is numeric
        try:
            parameters.setdefault(name, (eval(value)))
        except Exception as e:
            #assume paraemter is string
            try:
                parameters.setdefault(name, value)
            except Exception as e:
                print(e,'\n',name)

    parameters.setdefault('set_leafsize', int(parameters['set_leafsize_factor'] * parameters['num_neighbors']) )
    input_parameter_names=  parameters.keys()
    return parameters

def parameter_values(parameters, *args):
#     print(args)
    return args


def compute_correction(x,r,ineighbors,pts,radii_scaled):
    print("COMPUTE CORRECTION")

    # if(len(ineighbors)>1):
    #     1;
    # else:
    #     print(len(ineighbors))

    xn = pts[ineighbors]
    # print("xn")
    # print(xn)
    rn = radii_scaled[ineighbors]
    # print("rn")
    # print(rn)
    # print(xn.shape,x.shape)
    dxn = (xn - x[None,:]) 
    dn = np.linalg.norm(x[None,:]-xn,2,axis=1)
    # print("dn, dxn")
    # print(dn, dxn.shape)
    # print(rn.shape,r.shape)

    correct = (r+ rn)
    # print("correct")
    # print(correct)
    # correct = correct[:,None] - dxn * dn[:,None]
    # correct = correct / dn[:,None]

    correct = correct[:,None] - dxn / dn[:,None]
    correct = correct * dn[:,None]
    for c in correct:
        x += c
    # correct = np.sum(correct, axis=0)
    print("correct")
    print(x)

    return x
# xn = pts[ineighbors]
# rn = radii_scaled[ineighbors]
# dn = np.linalg.norm(x[:,None]-xn,2,axis=0)
# correct = (r[:,None] + rn) - (xn - x[:,None]) * 2 * dn
# correct = np.sum(correct, axis=0)
# return correct

def overlap_correction(i, x, r, pts,radii_scaled, kdt, registered, unregistered, dmax, search_radius_factor_of_max_diameter=2, pnorm=2, eps=1e-2):       
    search_radius = search_radius_factor_of_max_diameter * dmax
    ineighbors = kdt.query_ball_point(x, search_radius, pnorm, eps)
    ineighbors = np.array(ineighbors) 
    
    if(ineighbors.shape[0] < 2 ):
        reinit_flag=0
        # return x,0

    print(ineighbors.shape[0])
    if(ineighbors.shape[0]>2):
        print(pts.shape)
        dist_collide = np.linalg.norm(x[None,:] - pts[ineighbors], pnorm, axis=1)
        my_collision_length = radii_scaled[ineighbors] + r
        icollide = (dist_collide < my_collision_length)
        icollidesum = (icollide).sum()

    else:
        dist_collide = np.linalg.norm(x[:]-pts[ineighbors], pnorm)
        icollidesum=0
    try:
        if( icollidesum < 1 ):
            reinit_flag = 0
            # return x, reinit_flag
        else:
            print("CORRECT")
            icollide = np.where(icollide)
            print(icollide)
            print(ineighbors)

            icollneigh = icollide
            # icollneigh = np.intersect_1d(icollide,ineighbors)
            print(icollneigh)

            # correct = compute_correction(x,r,icollide,pts,radii_scaled)
            x = compute_correction(x,r,icollide,pts,radii_scaled)
            print("CORRECTION COMPUTED")
            reinit_flag=0

            # x = x + correct
            # x,reinit_flag = overlap_correction(i, x, r, pts,radii_scaled, kdt, registered, unregistered, dmax, search_radius_factor_of_max_diameter, pnorm, eps)
        
    except Exception:
        print(' \n \n Exception During Collision Check: Re-initialize with new random seed \n \n ')
    return x, reinit_flag
# print('x ', x)
# print('icollidesum ', icollidesum)
# print('r ', r)
# print('registered ', len(registered))
# print('unregistered ', len(unregistered))
# reinit_flag = 1
# return x, reinit_flag

# return correct,0


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


def where_boundary_intersect(parameters, pts,radii_scaled):
    xmin,xmax=parameters['xmin'],parameters['xmax']
    ymin,ymax=parameters['ymin'],parameters['ymax']
    zmin,zmax=parameters['zmin'],parameters['zmax']
    
    cond = (pts**2).sum(axis=1)-radii_scaled**2 > xmax**2
    cond = cond + ((pts**2).sum(axis=1)-radii_scaled**2 < xmin**2)
    cond = cond + ((pts**2).sum(axis=1)+radii_scaled**2 > ymax**2)
    cond = cond +((pts**2).sum(axis=1)-radii_scaled**2 < ymin**2)
    cond = cond + ((pts**2).sum(axis=1)+radii_scaled**2 > zmax**2)
    cond = cond +((pts**2).sum(axis=1)-radii_scaled**2 < zmin**2)
    return cond, cond.sum()


def get_reflected_pts(pts,idx_points, boundary, boundary_indices, boundary_radii, radii_scaled,xmin,xmax):
    sz0 = pts.shape[0]
    ptsr = pts**2 - radii_scaled[:,None]**2 
    #indices of spheres crossing lower boundary
    ip1 = np.any(ptsr < xmin,axis=1)
    #image spheres
    pts1 = np.mod(ptsr[ip1],xmax-xmin)

    # boundary = np.vstack((boundary, pts1.copy()))
    # boundary_indices = np.vstack((boundary_indices, ip1.copy()));
    # boundary_radii = np.vstack((boundary_radii, radii_scaled[ip1]))

    boundary.extend(pts1.tolist())
    boundary_indices.extend(ip1.tolist())
    boundary_radii.extend(radii_scaled[ip1].tolist())

    #tack on image spheres
    pts = np.vstack((pts,pts1.copy()))
    #tack on image indices
    if(len(ip1)>0):
        print(idx_points.shape)
        idx_points = np.hstack((idx_points,ip1[np.where(ip1)[0]]))
    #append image sphere radii
    radii_scaled = np.concatenate((radii_scaled,radii_scaled[ip1]))

    #DO IT AGAIN FOR THE UPPER BOUNDARIES
    ptsr = pts**2 + radii_scaled[:,None]**2
    ip1 = np.any(ptsr > xmax,axis=1)
    pts1 = np.mod(ptsr[ip1],xmax-xmin)

    # boundary = np.vstack((boundary, pts1.copy()))
    # boundary_indices = np.vstack((boundary_indices, ip1.copy()));
    # boundary_radii = np.vstack((boundary_radii, radii_scaled[ip1]))

    boundary.extend(pts1.tolist())
    boundary_indices.extend(ip1.tolist())
    boundary_radii.extend(radii_scaled[ip1].tolist())


    pts = np.vstack((pts,pts1.copy()))
    if(len(ip1)>0):
        print(idx_points.shape)
        idx_points = np.hstack((idx_points,ip1[np.where(ip1)[0]]))
    radii_scaled = np.concatenate((radii_scaled,radii_scaled[ip1]))
    sz = np.unique( pts ,axis=0).shape[0]
    print(sz,sz0)
    flag = 1 if( sz==sz0) else 0
        
    return pts, radii_scaled, idx_points, boundary, boundary_indices, boundary_radii, flag


def Run_Correlated_Sphere_Packing(input_parameters_filename="Parameters.in", seed_increment = 0, seed=None,periodic_geometry=None,nsamples=None,target_porosity=None):
    parameters = read_input_parameters(input_parameters_filename)
    reinit_flag=0

    '''debugging, ignore this'''
    print(','.join(np.array([name for name in parameters.keys()])))

    try:
        if(seed is not None):
            parameters['seed'] = seed
    except Exception as e:
        print('no seed specified')

    try:
        if(periodic_geometry is not None):
            parameters['periodic_geometry'] = periodic_geometry
    except Exception as e:
        print('no periodicv geometry specified')

    try:
        if(nsamples is not None):
            parameters['nsamples'] = nsamples
    except Exception as e:
        print('no nsampels specified')

    try:
        if(target_porosity is not None):
            parameters['target_porosity'] = target_porosity
    except Exception as e:
        print('no porosity specified')



    '''increment in case of reinitialization'''
    print(" seed_increment ", seed_increment, " type ", type(seed_increment))
    parameters['seed']+=seed_increment

    '''use input parameters'''
    periodic_geometry,ndimensions,xmin,xmax,ymin,ymax,zmin,zmax,seed,radius_mu,radius_sig2,Cmu,Csig,filename,radii_dist,nsamples,target_porosity,nbins,num_neighbors,set_leafsize_factor,kdt_eps,pnorm,search_radius_factor_of_max_diameter,find_all_neighbors,percentilemin,percentilemax,show_hist,set_dt,blfac,damping,tstep,tprint,tstepmax,force_abstol,fmin,force_absmax_factor,force_reltol,max_overlap_factor,set_leafsize_factor \
    =parameter_values(parameters, *parameters.values())
#     print(xmin,xmax)
    RandomState = np.random.RandomState(seed)

    '''form uniform hexahedral mesh'''
    ncells = nsamples
    xx=np.linspace(xmin,xmax,ncells+1)
    dcell = np.diff(xx)[0]
    domain_volume = (xmax-xmin) * (ymax-ymin) * (zmax-zmin)

    periodic_bounds = np.array([xmax,ymax,zmax])[:ndimensions]



#lognormal
# probability distr. fcn.
# cumulative distr. fcn.
# inverse cdf 
# q05 = np.exp( )
#get rmin, rmin+rrange=:rmax given quantiles
#with quantiles, get truncated CDF
# given mean and variance, have median

    #PART I
    '''sample radii'''
    # #sample radii 
    # while target_porosity < v0:
    # v0 = (radii**3 * 4*np.pi/(3 * (xmax-xmin)**3))
    # v0 = v0/v0.sum()
    if(radii_dist=='lognormal'):
        # V prior error, fix: mean of lognormal is mean of associated normal distribution

        #lognormal:
        # mean: exp(mu+sig2/2) 
        #median exp(mu) 
        #mode exp(mu - sig2) 
        #variance = mean^2 (-1+E^\[Sigma]^2)}
        #quantile q = E^(\[Mu] - Sqrt[2] \[Sigma]^2 InverseErfc[2 q])
        #Quantiles 
        # q05 = E^(\[Mu] - Sqrt[2]*E^(\[Mu] - Sqrt[2]*\[Sigma]^2*InverseErfc[1/10])*\[Sigma])
        # q95 = E^(\[Mu] - Sqrt[2]*E^(\[Mu] + Sqrt[2]*\[Sigma]^2*InverseErfc[1/10])*\[Sigma])} - -E^(\[Mu] - Sqrt[2]*E^(\[Mu] - Sqrt[2]*\[Sigma]^2*InverseErfc[1/10])*\[Sigma])

        # lognormal_sig2 = np.log( radius_sig2 / radius_mu**2 + 1 )
        # lognormal_mu = np.log( radius_mu ) - lognormal_sig2 / 2

        # Z = RandomState.lognormal( lognormal_mu , lognormal_sig2 ,nsamples)
        Z = RandomState.lognormal( np.log(radius_mu) , radius_sig2 ,nsamples)
        # lognormal_sig2 = np.log( radius_sig2 / radius_mu**2 + 1 )
        # lognormal_mu = np.log( radius_mu ) - lognormal_sig2 / 2
        # Z = RandomState.lognormal( lognormal_mu, lognormal_sig2 ,nsamples)
        radii = np.exp(Z)

    if(show_hist==True):
        plt.hist(radii,nbins)
        plt.savefig('hist.png')

    pmin = np.percentile(radii,percentilemin),
    pmax = np.percentile(radii,percentilemax)
    radii = radii[radii<pmax]
    radii = radii[radii>pmin]
    nsamplesclip = nsamples - radii.shape[0]
    nsamples = radii.shape[0]# - nsamplesclip
    radii = np.sort(radii)


    
    
    print('1,99 percentiles ', pmin,pmax)
    print('max, min, mean', radii.max(),radii.min(),radii.mean())



    '''rescale radii to obtain desired porosity'''
    rscale = ( (  np.sum(4*np.pi*(1/3.) * radii**3) ) / ( domain_volume * (1-target_porosity) ) )**(1/3)
    radii_scaled = radii/rscale
    radius_mu_scaled = radius_mu / rscale
    radius_sig2_scaled = radius_sig2 / rscale**2
    rmax = ( radii_scaled ).max()
    dmax = 2*rmax
    delta = dcell - 2*rmax
    search_radius = 2*dmax#+delta

    # print("IINTERIOR " , iinterior)


    #TBD
    '''sample isotropic gaussian correlation'''
    Cnorm = RandomState.multivariate_normal(Cmu,Csig,nsamples)#,[nsamples,ndimensions])

    '''sample points uniformly in space '''
    pts = RandomState.uniform(0,1,[nsamples,ndimensions])
    pts[:,0] = (pts[:,0] ) * (xmax-xmin) + xmin
    pts[:,1] = (pts[:,1] ) * (ymax-ymin) + ymin
    pts[:,2] = (pts[:,2] ) * (zmax-zmin) + zmin


    '''sort radii and points '''
    radii_scaled = radii_scaled[::-1]
    pts = pts[::-1,:]
    
    '''get nearest neighbor distances'''
    t0 = time.time()

    if(periodic_geometry==True):
        kdt = PeriodicCKDTree(periodic_bounds, pts)
    else:
        kdt = scipy.spatial.KDTree(pts) #,leafsize=set_leafsize_factor * num_neighbors )
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


    '''BC, EQ separation, pore throat size, collision distances, PD'''
    boundary_layer = [blfac, xmax *(1- blfac)]
    # IINT
    print("BOUINDARY LAYER" , boundary_layer)
    iboundary = np.any(((pts<boundary_layer[0]).astype(int)+(pts>boundary_layer[1]).astype(int)) , axis=1)
    iinterior = np.all((pts>boundary_layer[0]).astype(int) * (pts<boundary_layer[1]).astype(int),axis=1)
    # print("IBOUNDARY " , iboundary)
    # print("IINTERIOR " , iinterior)
    print('\n \n \n ' , " NUMBER OF BOUNDARY SPHERES ", iboundary.sum(), " TOTAL CONSIDERED ", iboundary.shape, '\n \n \n ' , "NUMBER OF INTERIOR SPHERES", iinterior.sum(), " TOTAL CONSIDERED ",  iinterior.shape, '\n \n \n ')
    # ,iboundary.sum()




    ''' Detect Collisions and Translate Spheres '''
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
            pts[i] = x
            radii_scaled[i] = r
        else:
            x,reinit_flag = overlap_correction(i, x, r, pts, radii_scaled, kdt, registered, unregistered, dmax, search_radius_factor_of_max_diameter, pnorm=2, eps=kdt_eps)
            if(reinit_flag==1):
                if(periodic_geometry==1):
                    break;
            elif(reinit_flag==0):                    
                registered.append(i)
                pts[i] = x
                radii_scaled[i] = r
        unregistered.remove(i)

        t_list.append(time.time() - tlast)
        tlast = time.time()

    # if(reinit_flag==1):
    #     print(" \n Reinitializing Simulation \n")
    #     return Run_Correlated_Sphere_Packing(seed_increment+1)
    # else:
    #     print("\n No collisions found, continuing.. \n")
    t_list = np.array(t_list)

    registered = np.array(registered)




    print("STORE BOUNDARY SPHERE DATA")
    indices_boundary = np.where(iboundary)
    print("STORE INTERIOR SPHERE DATA")
    indices_interior = np.where(iinterior)

    interior_points = pts[indices_interior]
    boundary_points = pts[indices_boundary]

    print("SET POINT IDs (to keep track of boundary image spheres)")
    idx_points = np.arange(0,len(registered))



    boundary, boundary_indices, boundary_radii = [], [], [] #np.array([]), np.array([]), np.array([])
    if(periodic_geometry==1):
        print("COPY BOUNDARY POINTS TO IMAGE SPHERES ACROSS PERIODIC BOUNDARIES")
        print("NUM POINTS BEFORE BOUNDARY IMAGE COPY" , pts.shape)
        flag=0
        # while flag==0:
        # pts, radii_scaled,idx_points, flag = get_reflected_pts(pts,idx_points, radii_scaled,xmin,xmax)


        pts, radii_scaled, idx_points, boundary, boundary_indices, boundary_radii, flag = get_reflected_pts(pts,idx_points, boundary, boundary_indices, boundary_radii, radii_scaled,xmin,xmax)
        print(radii_scaled.shape)
        nsamples = radii_scaled.shape[0]
        assert(radii_scaled.shape[0] == pts.shape[0])
        print("NUM POINTS AFTER " , pts.shape)

    # pts = pts[registered]
    # radii_scaled = radii_scaled[registered]
    
    # print(np.all((pts[0]<xmax , pts[1]<ymax , pts[2]<zmax , pts[0]>xmin , pts[1]>ymin , pts[2]>zmin) , axis=2))
    # idbounds = np.where( np.all(pts[0]<xmax , pts[1]<ymax , pts[2]<zmax , pts[0]>xmin , pts[1]>ymin , pts[2]>zmin))[0]
    # print(idbounds)
    # print(idx_points.shape, pts.shape, radii_scaled.shape)

    # seed=parameters['seed']
    # nsamples = parameters['nsamples']
    # ndimensions = parameters['ndimensions']
    # print('seed', seed)
    save_filename = 'packing'
    #save packing output
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
    return parameters, radii_scaled, registered, unregistered, pts, pvolumes, idx_points, boundary, boundary_indices, boundary_radii



def main():
    argumentList=sys.argv[1:]
    unixOptions = "f:s:g:n:p:"  
    gnuOptions = ["seed", "periodic_geometry", "nsamples","target_porosity"]
    print(" \n \n \n ")

    # seed

    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions)
    except getopt.error as e:  
        print(str(e))
        sys.exit(2)

    # print("arguments", arguments)
    # print("values", values)

    seed,periodic_geometry,nsamples,target_porosity = None,None,None,None
    seed_increment=0


    for currentArgument, currentValue in arguments:  
        if currentArgument in ("-f", "--file"):
            print ("input_parameters_filename ", currentValue)
            input_parameters_filename = currentValue
        else:
            input_parameters_filename = "Parameters.in"

        if currentArgument in ("-s", "--seed","--randomstate"):
            print ("seed ", currentValue)
            seed = int(currentValue)
        if currentArgument in ("-g", "--periodic_geometry"):
            print ("periodic_geometry ", currentValue)
            periodic_geometry = bool(currentValue)
        if currentArgument in ("-n", "--nsamples","--nspheres"):
            print ("nspheres ", currentValue)
            nsamples = int(currentValue)
        if currentArgument in ("-p", "--porosity","--target_porosity"):
            print ("target_porosity ", currentValue)
            target_porosity = float(currentValue)

    print(" \n \n \n ")

    # try:
    # 	parameters, radii_scaled, registered, unregistered, pts, pvolumes = Run_correlated_Sphere_Packing(input_parameters_filename, seed_increment, seed,periodic_geometry,nsamples,target_porosity)
    # except Exception as e:
    parameters, radii_scaled, registered, unregistered, pts, pvolumes, idx_points, boundary, boundary_indices, boundary_radii = Run_Correlated_Sphere_Packing('Parameters.in', seed_increment, seed,periodic_geometry,nsamples,target_porosity)

    print(idx_points.shape, pts.shape, radii_scaled.shape)
    seed=parameters['seed']
    nsamples = parameters['nsamples']
    ndimensions = parameters['ndimensions']
    print('seed', seed)
    save_filename = 'test_packing'
    #save packing output
    stacked_data = np.vstack((idx_points.astype(int), pts[:,0], pts[:,1], pts[:,2], radii_scaled[:])).T
    np.savetxt(save_filename + ".txt", stacked_data, header="ID x y z r", fmt='%i,%E,%E,%E,%E')


    savedict = {'points':pts,
               'radii_scaled':radii_scaled,
               'target_porosity':target_porosity,
                'boundary':boundary,
                'boundary_indices':boundary_indices, 
                'boundary_radii':boundary_radii
               }
    filename = 'uncorr_pack_dim_'+str(ndimensions)+'_nparticles_'+str(nsamples)
    with gzip.open(filename + '_.pkl','wb') as f:
        pickle.dump(savedict,f)



    psize = ((pvolumes-pvolumes.min()))
    psize = psize/psize.max()
    psize = psize * nsamples



    fig = plt.figure(figsize=(16,16))
    ax = fig.add_subplot(111, projection='3d')
    X,Y,Z = pts[:,0],pts[:,1],pts[:,2]
    cube = ax.scatter(X, Y, Z, zdir='z', s=psize, c=psize,alpha=0.4)
    cbar = fig.colorbar(cube, shrink=0.6, aspect=5)
    plt.savefig(save_filename+'_3d.png',dpi=200)
    # plt.show()



    # len(),nsamples
    plt.figure(figsize=(16,16))
    plt.scatter(pts[:,0],pts[:,1],s=psize,c=psize,alpha=0.6)
    plt.colorbar()
    plt.savefig(save_filename+'_.png',dpi=200)
    # plt.show()


    # plt.figure(figsize=(16,8))
    # plt.plot(t_list,'.',alpha=0.4,ms=4);
    # plt.ylabel('Time Elapsed Per Iteration ')
    # plt.xlabel('Iteration (Sequential Addition)')
    # plt.savefig(filename+'_tiame.png',dpi=200)
    # plt.show()






if __name__ == "__main__":
    main()



