import numpy as np

pts = 100

par1_name = "K"
par1_start = 0.0
par1_stop = 8.0
par1_interval = par1_stop - par1_start 
par1_incr = par1_interval / pts


par2_name = "tau"
par2_start = 0.0
par2_stop = 5.0

par2 = np.linspace(par2_start, par2_stop, pts, endpoint=False)




pxl_per_job = 10
jobs_per_arfile = 1000


base_str1 = "-m:lscan[" + par1_name + ","
base_str2 = ","
base_str3 = "," + str(pxl_per_job) + "," + par2_name + "] -" + par2_name + " "
base_str4 = " -rea 100 -D 0.2 -int_time 1500 -out_time 100 -wLP " +str(np.pi*0.75) 



start_list = []
stop_list = []

start = par1_start
stop  = par1_start


for i in range(pts):
	stop += par1_incr
	if i%pxl_per_job == pxl_per_job-1:
		start_list.append(start)
		stop_list.append(stop)
		start = stop
	elif i == pts -1:
		start_list.append(start)
		stop_list.append(par1_stop)



str_list = []

for j in range(pts):
	for i in range(len(start_list)):
		s = base_str1 + str(start_list[i])
		s+= base_str2 + str(stop_list[i])
		s+= base_str3 + str(par2[j])
		s+= base_str4
		str_list.append(s)

argfile_list = []

arg_str = ""

for i in range(len(str_list)):
	arg_str += str_list[i]
	arg_str += '\n'
	if i%jobs_per_arfile == jobs_per_arfile-1:
		argfile_list.append(arg_str)
		arg_str = ""
	elif i == len(str_list) - 1:
		argfile_list.append(arg_str)


for i in range(len(argfile_list)):
	f = open("argfile"+str(i), "w")
	f.write(argfile_list[i])
	f.close()
