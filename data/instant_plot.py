import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import griddata

file_name = "test_hom_D0.2.ts.dat"


data = np.loadtxt(file_name)


data = np.asarray(list(zip(*data)))


t  = data[0]
ER = data[1]
EI = data[2]
G = data[3]
Q = data[4]

I = data[6]

fig, (ax1, ax2, ax3) = pl.subplots(3)

ax1.plot(t, ER)
ax1.plot(t, EI)

ax2.plot(t, I)
ax2.plot(t, np.full(len(t), 1.68488))

t_scat = []
I_scat = []


index=[
3043,4104,
13226,14140,
23396,24249,
33448,34459,
43518,44595,
53660,54672,
63801,64766,
73920,74876,
84041,84938,
94144,95002,
104299,105171,
114403,115357,
124621,125464,
134703,135589,
144765,145795,
155036,156053,
165147,166168,
175223,176341,
185318,186422,
195492,196460,
205624,206520]

for i in range(0, len(index)):
	t_scat.append(t[index[i]])
	I_scat.append(I[index[i]])



ax2.scatter(t_scat, I_scat, color='r', zorder=3)

ax3.plot(t, G)
ax3.plot(t, Q)


#pl.plot(data)

pl.show()




#pl.savefig("bil_cc_reverse.ts.dat.png", dpi=350)
