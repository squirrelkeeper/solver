import numpy as np
import matplotlib.pyplot as pl

file_name = "test_hom_D0.ts.dat"


data = np.loadtxt(file_name)


data = np.asarray(list(zip(*data)))


t  = data[0]
ER = data[1]
EI = data[2]
G = data[3]
Q = data[4]

fig, (ax1, ax2, ax3) = pl.subplots(3)

ax1.plot(t, ER)
ax1.plot(t, EI)

ax2.plot(t, ER*ER+EI*EI)

ax3.plot(t, G)
ax3.plot(t, Q)


print(np.sqrt(sum(ER*ER+EI*EI)*0.0001)*0.2*25000)

#pl.plot(data)

pl.show()



#pl.savefig("bil_cc_reverse.ts.dat.png", dpi=350)
