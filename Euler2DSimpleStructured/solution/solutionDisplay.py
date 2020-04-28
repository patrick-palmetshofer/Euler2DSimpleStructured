import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

deg = "10"
it = ""
#1041.566000 555.500000 867.970000
vel = "555.500000"

gridpath = "../mesh/Grid"+deg+"deg.grd"
solpath = "raw/withlimiter/Grid"+deg+"deg"+vel+".res"
residualpath = "raw/withlimiter/Grid"+deg+"deg"+vel+"_residualsL2.csv"

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
residuals = pd.read_csv(residualpath,delimiter=",\t")


X =  grid.loc[:,0].values.reshape((nxi,neta))
Y =  grid.loc[:,1].values.reshape((nxi,neta))
data=[0]*4
for i in range(4):
    data[i] = solution.iloc[:,i].values.reshape((nxi+1,neta+1))


fig1, axes = plt.subplots(2,2)
axes = axes.reshape((4,1))
fig1.set_figheight(15)
fig1.set_figwidth(22)

cmap = "plasma"
titles = ["Density [kg/m^3]", "u [m/s]", "v [m/s]", "T [K]"]

for i,ax in enumerate(axes):
    d = data[i][1:-1,1:-1]
    ax = axes[i,0]
    plot = ax.pcolormesh(X,Y,d,cmap=cmap)
    ax.axis('equal')
    fig1.colorbar(plot, ax = ax)
    ax.set_title(titles[i])

# ADD CORRECT CELL MIDPOINTS
fig2, axes2 = plt.subplots(2,2)
axes2 = axes2.reshape((4,1))
fig2.set_figheight(15)
fig2.set_figwidth(22)
for i,ax in enumerate(axes2):
    d = data[i][1:,1:]
    ax = axes2[i,0]
    plot = ax.contourf(X,Y,d,cmap=cmap)
    ax.axis('equal')
    fig2.colorbar(plot, ax = ax)
    ax.set_title(titles[i])

R = 1005*(1-1/1.4)
gamma = 1.4
Ma = np.sqrt((data[1]**2+data[2]**2)/(gamma*R*data[3]))
figMa, axMa = plt.subplots()
figMa.set_figheight(15)
figMa.set_figwidth(22)
axMa.axis('equal')
plot = axMa.contourf(X,Y,Ma[1:,1:],cmap=cmap)
figMa.colorbar(plot, ax = axMa)
axMa.set_title("Mach number [-]")

figRes, axRes = plt.subplots()
figRes.set_figheight(15)
figRes.set_figwidth(22)
for col in residuals.columns:
    axRes.plot(residuals.index,residuals[col],label=col)
axRes.set_yscale("log")
axRes.set_ylabel("Normalized Residual")
axRes.set_ylabel("Iteration")
axRes.legend()
#axRes.axis([0,1000,1e-6,1])
axRes.grid()