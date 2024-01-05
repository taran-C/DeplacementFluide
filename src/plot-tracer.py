import numpy as np
from matplotlib import pyplot as plt, animation
from netCDF4 import Dataset
import os

dirname = os.path.dirname(__file__)
input_file = os.path.join(dirname, '../output/a.nc')

nc_file = Dataset(input_file, "r")

k = 4
q = nc_file.variables["q"][k,:,:,:]
i = range(len(q[:,0,0]))
j = range(len(q[0,:,0]))
t = 10
fig,ax = plt.subplots()

cax = ax.pcolormesh(i,j,q[:,:,0], vmin=-1, vmax=4)
fig.colorbar(cax)

def animate(t):
    cax.set_array(q[:,:,t])
    #print(q[:,:,t].max())
    #print(t)

animator = animation.FuncAnimation(fig, animate, interval=5, frames = len(q[0,0,:]))
animator.save('output/tracer.gif', dpi=200, writer='magick')
plt.show()