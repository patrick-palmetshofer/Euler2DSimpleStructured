import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

gridpath = "../mesh/Grid20deg.grd"
solpath = "raw/Grid20deg10000.res"

with open(gridpath) as f:
    first_line = f.readline()
	
print(first_line)

with open(gridpath) as f:
    first_line = f.readline()
	
print(first_line)

linevec = first_line.split()

nxi = int(linevec[0])
neta = int(linevec[1])
	


ncells = [int(vn) for vn in first_line.split()]

grid = pd.read_csv(gridpath,skiprows=1,header=None,delimiter="\s")
solution = pd.read_csv(solpath,skiprows=1,delimiter=",\t")


X =  grid.loc[:,0].values.reshape((nxi,neta))
Y =  grid.loc[:,1].values.reshape((nxi,neta))
data=[0]*4
for i in range(4):
    data[i] = solution.iloc[:,i].values.reshape((nxi-1,neta-1))


fig1, axes = plt.subplots(2,2)
axes = axes.reshape((4,1))
fig1.set_figheight(15)
fig1.set_figwidth(20)

plasma = cm.get_cmap('plasma', 256)
#newcolors = plasma(np.linspace(np.amin(data[2]), np.amax(data[2]), 256))
#cmap = ListedColormap(newcolors)

for i,ax in enumerate(axes):
    d = data[i]
    axes[i,0].pcolormesh(X,Y,d,cmap="plasma")
