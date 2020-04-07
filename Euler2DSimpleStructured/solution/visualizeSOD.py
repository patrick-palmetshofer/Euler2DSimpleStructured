import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

gridpath = "../mesh/sodX.grd"
solpath = "raw/sodX600.res"

with open(gridpath) as f:
    first_line = f.readline()
	
print(first_line)

with open(gridpath) as f:
    first_line = f.readline()
	
print(first_line)

linevec = first_line.split()

nxi = int(linevec[0])
neta = int(linevec[1])
	


npoints = [int(vn) for vn in first_line.split()]

grid = pd.read_csv(gridpath,skiprows=1,header=None,delimiter="\s")
solution = pd.read_csv(solpath,skiprows=1,delimiter=",\t")


X =  grid.loc[:,0].values.reshape((nxi,neta))
Y =  grid.loc[:,1].values.reshape((nxi,neta))
data=[0]*4
for i in range(4):
    data[i] = solution.iloc[:,i].values.reshape((nxi+1,neta+1))

fig1,axes = plt.subplots(2,2)
axes = axes.reshape((4,1))
fig1.set_figheight(15)
fig1.set_figwidth(20)
for i,ax in enumerate(axes):
    d = data[i][0:-1,1]
    axes[i,0].plot(X[:,0],d)
    
# gridpath = "../mesh/sodY.grd"
# solpath = "raw/sodY2.res"

# with open(gridpath) as f:
#     first_line = f.readline()
# 	
# print(first_line)

# with open(gridpath) as f:
#     first_line = f.readline()
# 	
# print(first_line)

# linevec = first_line.split()

# nxi = int(linevec[0])
# neta = int(linevec[1])
# 	


# npoints = [int(vn) for vn in first_line.split()]

# grid = pd.read_csv(gridpath,skiprows=1,header=None,delimiter="\s")
# solution = pd.read_csv(solpath,skiprows=1,delimiter=",\t")


# X =  grid.loc[:,0].values.reshape((nxi,neta))
# Y =  grid.loc[:,1].values.reshape((nxi,neta))
# data=[0]*4
# for i in range(4):
#     data[i] = solution.iloc[:,i].values.reshape((nxi+1,neta+1))

# fig1,axes = plt.subplots(2,2)
# axes = axes.reshape((4,1))
# fig1.set_figheight(15)
# fig1.set_figwidth(20)
# for i,ax in enumerate(axes):
#     d = data[i][1,0:-1]
#     axes[i,0].plot(Y[0,:],d)


# fig1, axes = plt.subplots(2,2)
# axes = axes.reshape((4,1))
# fig1.set_figheight(15)
# fig1.set_figwidth(20)

# plasma = cm.get_cmap('plasma', 256)
# #newcolors = plasma(np.linspace(np.amin(data[2]), np.amax(data[2]), 256))
# #cmap = ListedColormap(newcolors)

# for i,ax in enumerate(axes):
#     d = data[i]
    # axes[i,0].pcolormesh(X,Y,d,cmap="plasma")
