import os
import numpy as np
import matplotlib.pyplot as pl

def get_files(suffix,path='./',m=True):
	all_files = os.listdir(path)
	dat_files = list(filter(lambda x: x[-len(suffix):] == suffix, all_files))
	dat_files.sort()
	if m:
		print(str(len(dat_files))+ " files found")
	return dat_files


file_path = "./"
files_dat = get_files(".ts.dat")


for i in range(len(files_dat)):
	file_name = file_path + files_dat[0]

	data = np.loadtxt(file_name)
	data = list(zip(*data))

	t = np.asarray(data[0])

	I = np.asarray(data[1])


	pl.plot(t, I, c='k', linewidth=0.5)



	pl.show()


