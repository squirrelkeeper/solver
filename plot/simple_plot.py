import os
import numpy as np
import matplotlib.pyplot as pl

data_path = "../data/"
data_name = "template_ts_ba"
data_suffix = ".dat"

files_all = os.listdir(data_path)
files_dat = list(filter(lambda x: x[-4:] == ".dat", files_all))

file_name = data_path + files_dat[1]

data = np.loadtxt(file_name)
data = list(zip(*data))

t = np.asarray(data[0])

ER = np.asarray(data[1])
EI = np.asarray(data[2])

G = np.asarray(data[3])
Q = np.asarray(data[4])
J = np.asarray(data[5])

I = np.asarray(data[6])

'''
t = t[len(I)-int(len(I)/10)::]
G = G[len(I)-int(len(I)/10)::]
Q = Q[len(I)-int(len(I)/10)::]
I = I[len(I)-int(len(I)/10)::]
'''

#I_max  = np.full(len(t), 8.15138)
#I_mid  = np.full(len(t), 4.07569)
#I_mean = np.full(len(t), 0.310664)
'''
pl.plot(t, I_max, c='grey', linewidth=0.5)
pl.plot(t, I_mid, c='grey', linewidth=0.5)
pl.plot(t, I_mean, c='grey', linewidth=0.5)
'''

#thres_front  = np.full(len(t), (4.07569 - 0.310664)*0.5)
#thres_tail = np.full(len(t), (4.07569 - 0.310664)*0.25)


#pl.plot(t, thres_front, c='grey', linewidth=0.5)
#pl.plot(t, thres_tail, c='grey', linewidth=0.5)


pl.plot(t, ER*ER+EI*EI, c='k', linewidth=0.5)
#l.plot(t, G, c='r', linewidth=0.5)
#pl.plot(t, Q, c='b', linewidth=0.5)



pl.show()




#pl.savefig("lina_out.ts.dat.png", dpi=350)
#pl.savefig("lasse_out.ts.dat.png", dpi=350)
