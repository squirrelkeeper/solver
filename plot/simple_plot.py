import numpy as np
import matplotlib.pyplot as pl

data = np.loadtxt("../data/template_ts_ba9257.dat")
data = list(zip(*data))

t = np.asarray(data[0])

ER = np.asarray(data[1])
EI = np.asarray(data[2])

G = np.asarray(data[3])
Q = np.asarray(data[4])
J = np.asarray(data[5])

I = np.asarray(data[6])

t = t[len(I)-int(len(I)/10)::]
G = G[len(I)-int(len(I)/10)::]
Q = Q[len(I)-int(len(I)/10)::]
I = I[len(I)-int(len(I)/10)::]

pl.plot(t, I, c='k', linewidth=0.5)
pl.plot(t, G, c='r', linewidth=0.5)
pl.plot(t, Q, c='b', linewidth=0.5)




#pl.plot(t, G, c='b', linewidth=0.5)
#pl.plot(t, Q, c='cyan', linewidth=0.5)



#pl.show()




#pl.savefig("lina_out.ts.dat.png", dpi=350)
pl.savefig("lasse_out.ts.dat.png", dpi=350)
