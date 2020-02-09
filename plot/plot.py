import numpy as np
import matplotlib.pyplot as pl

data = np.loadtxt("ts_data.dat")

data = list(zip(*data))

t = np.asarray(data[0])

ER = np.asarray(data[1])
EI = np.asarray(data[2])

G = np.asarray(data[3])
Q = np.asarray(data[4])
J = np.asarray(data[5])

N = np.asarray(data[5])

I = (ER*ER+EI*EI)/(N*N)

pl.plot(t, I, c='k', linewidth=0.5)
#pl.plot(t, EI, c='magenta', linewidth=0.5)



#pl.plot(t, G, c='b', linewidth=0.5)
#pl.plot(t, Q, c='cyan', linewidth=0.5)



pl.show()




#pl.savefig("test_out.ts.dat.pdf", dpi=500)
