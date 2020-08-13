import numpy as np





other_pars = ""

pts = 100

par1_name = "K"
par1_start = 0.0
par1_stop = 5.0
par1_interval = par1_stop - par1_start 
par1_incr = par1_interval / pts


par2_name = "tau"
par2_start = 0.0
par2_stop = 5.0



par2 = np.linspace(par2_start, par2_stop, pts)



pxl_per_job = 10
jobs_per_line = int(np.ceil(pts/pxl_per_job))
jobs_per_arfile = 1000

pxl = pts*pts
jobs = int(np.ceil(pts/pxl_per_job)*pts)

argfile_num = int(np.ceil(jobs/1000))




args_list = []


for j in range(pts):
	for i in range(jobs_per_line):
		this_start = pxl_per_job * i * par1_incr
		this_stop = pxl_per_job * (i+1) * par1_incr
		this_pts = pxl_per_job
		if i == jobs_per_line-1:
			this_start = pxl_per_job * i * par1_incr
			this_stop = par1_stop
			this_pts = pts%pxl_per_job
		if this_pts != 0:
			args = "-m:lscan["
			args += par1_name
			args += ","
			args += str(this_start)
			args += ","
			args += str(this_stop)
			args += ","
			args += str(this_pts)
			args += ","
			args += par2_name
			args += "] -"
			args += par2_name
			args += " "
			args += str(j*(par2_stop-par2_start)/pts + par2_start)
			args += " -rea "
			args += str(100)
			args += " -int_time "
			args += str(1500)
			args += " -out_time "
			args += str(1000)
			args += '\n'
			args_list.append(args)


args = ""





argfile_list = []


for j in range(argfile_num):
	for i in range(jobs_per_arfile):
		if i + jobs_per_arfile * j == len(args_list)-1:
			break
		args += args_list[i + jobs_per_arfile * j]
	argfile_list.append(args)
	args = ""



for i in range(len(argfile_list)):
	f = open("argfile"+str(i), "w")
	f.write(argfile_list[i])
	f.close()



'''
for j in range(argfile_num):
	args = ""
	for i in range(jobs_per_arfile):
		if i+j*jobs_per_arfile == len(args_list) - 1:
			break
		args += args_list[i+j*jobs_per_arfile]
		if i == len(args_list)-1:
			break
	

'''	



	



