#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 13:13:16 2018

@author: chaztikov
"""
import revised_packing_pbc.py
import os,sys
cwd=os.getcwd()
sys.path.append(cwd)
import numpy as np
import scipy
import scipy.stats

import matplotlib as mpl
import matplotlib.pyplot as plt
import run_test_packing
from run_test_packing import *

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

N = 1000
trace = go.Scattergl(
    x = np.random.randn(N),
    y = np.random.randn(N),
    mode = 'markers',
    marker = dict(
        color = '#FFBAD2',
        line = dict(width = 1)
    )
)
data = [trace]
py.iplot(data, filename='packing_revised.txt',fileopt='new')