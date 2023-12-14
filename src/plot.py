import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os

dirname = os.path.dirname(__file__)
input_file = os.path.join(dirname, '../output/a.nc')

nc_file = Dataset(input_file, "r")

x = nc_file.variables["x(t)"][0,0,0,:]
y = nc_file.variables["y(t)"][0,0,0,:]
z = nc_file.variables["z(t)"][0,0,0,:]

ax = plt.axes(projection='3d')

ax.plot(x,y,z, lw = 0.5)

plt.show()

