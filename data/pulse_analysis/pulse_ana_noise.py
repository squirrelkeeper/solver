import numpy as np
import matplotlib.pyplot as pl

def binning(t, I, n):
	dt = (max(t)-min(t))/len(t)
	
	i_max = int(1/dt/n)
	j_max = int(np.floor(len(I)/i_max))
	
	bin_t = np.zeros(j_max-1)
	bin_I = np.zeros(j_max-1)
	
	for j in range(0, j_max-1):
		bin_t[j] = np.mean(t[j*i_max:(j+1)*i_max])
		bin_I[j] = np.mean(I[j*i_max:(j+1)*i_max])

	return bin_t, bin_I

def sliding_win(t, I):
	print(8)


file_name = "test_hom_D0.2.ts.dat"


data = np.loadtxt(file_name)


data = np.asarray(list(zip(*data)))


t  = data[0]
ER = data[1]
EI = data[2]
G = data[3]
Q = data[4]
J = data[5]
I = data[6]


bin_t, bin_I = binning(t, I, 4)


trig_line = np.full(len(bin_t), max(I)*0.5)

fig, (ax1, ax2) = pl.subplots(2)

ax1.plot(t, I, linewidth=0.5)

ax1.plot(bin_t, trig_line, linewidth=0.5)


trig_line = np.full(len(bin_t), max(bin_I)*0.5)


ax2.plot(bin_t, bin_I, linewidth=0.5)
ax2.plot(bin_t, trig_line, linewidth=0.5)


#ax2.plot(t, G, linewidth=0.5)
#ax2.plot(t, Q, linewidth=0.5)


#pl.plot(data)

pl.show()



#pl.savefig("bil_cc_reverse.ts.dat.png", dpi=350)

