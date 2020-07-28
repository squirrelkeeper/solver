import os
import numpy as np
import matplotlib.pyplot as pl

from scipy.interpolate import griddata
from matplotlib import colors



path = "./"
all_files = os.listdir(path)
dat_files = list(filter(lambda x: x[-4:] == '.dat', all_files))


file_name = "out_sweep.swp.dat"

data = np.loadtxt(file_name)
data = np.asarray(list(zip(*data)))


x = data[0]
y = data[1]
z = data[2]

grid_x, grid_y = np.meshgrid(np.linspace(min(x),max(x),100),np.linspace(min(y),max(y),100))

points = np.array(list(zip(x,y)))

gridded = griddata(points, z, (grid_x, grid_y), method='nearest')

gridded = np.flip(gridded, axis=0)







cmap = colors.ListedColormap(['white', 'red', 'blue', 'black'])
bounds=[0.5,1.5,2.5,3.5,1500]
norm = colors.BoundaryNorm(bounds, cmap.N)







fig = pl.imshow(gridded, extent=(min(x),max(x),min(y),max(y)), cmap=cmap, norm=norm, interpolation="none", aspect="auto")

pl.colorbar(fig, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0, 1, 2])

pl.xlabel(r"$J_g$")
pl.ylabel(r"$q_0$")


pl.show()



'''
plt.clim(0,20.0)

cb1 = plt.colorbar()



cb1.set_label(r"$\sigma_{lt}$ [fs]")

plt.title(r"$\alpha_g=1$, $\alpha_q=1$")

plt.savefig("lttj_Jg_q0_ag1 _aq1.png",dpi=500)
'''
