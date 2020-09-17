import numpy as np
import os
import matplotlib.pyplot as pl
from scipy.interpolate import griddata


def get_files(suffix,path='./',m=True):
	all_files = os.listdir(path)
	dat_files = list(filter(lambda x: x[-len(suffix):] == suffix, all_files))
	dat_files.sort()
	if m:
		print(str(len(dat_files))+ " files found")
	return dat_files




dat_files = get_files(".ts.dat")


file_name = dat_files[0]

print(file_name)
start = 0
end = 95000



data = np.loadtxt(
	file_name,
	max_rows = end
)


data = np.asarray(list(zip(*data)))



t  = data[0]
ER = data[1]
EI = data[2]
G  = data[3]
Q  = data[4]
J  = data[5]

I  = data[6]

fig, (ax1, ax2) = pl.subplots(2)

ax1.plot(
	t,
	I,
	linewidth = 0.5,
	color = 'k',
)

ax1.plot(
	t,
	G,
	linewidth = 0.5,
	color = 'b',
)

ax1.plot(
	t,
	Q,
	linewidth = 0.5,
	color = 'r',
)

ax1.plot(
	t,
	J,
	linewidth = 0.5,
	color = 'g',
)


ax2.plot(t, ER)
ax2.plot(t, EI)

#pl.plot(data)

pl.show()



#pl.savefig("bil_cc_reverse.ts.dat.png", dpi=350)
