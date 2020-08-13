import numpy as np





other_pars = "-K 0.0 -D 0.2 -rea 100 "

pts = 101

par1_name = "g"
par1_start = 30.0
par1_stop = 100.0

par1 = np.linspace(par1_start, par1_stop, pts)


par2_name = "Jg"
par2_start = 2.5
par2_stop = 5.0

par2 = np.linspace(par2_start, par2_stop, pts)

args_list = []

total_counter = 0
file_counter = 1

args = ""

for i in range(len(par1)):
	for j in range(len(par2)):
		total_counter += 1
		
		args += other_pars
		args += "-" + par1_name + " "
		args += str(par1[i]) + " "
		
		args += "-" + par2_name + " "
		args += str(par2[j]) + '\n'
		
		if(total_counter == 999 or (i == len(par1)-1 and j == len(par2)-1)):
			file_counter += 1
			total_counter = 0
			
			args_list.append(args)
			args = ""


for i in range(len(args_list)):
	f1 = open("argfile" + str(i+1), "w")
	f1.write(args_list[i])
	f1.close()


'''		
print(argfile3)

f1 = open("argfile1", "w")
f1.write(argfile1)
f1.close()

f2 = open("argfile2", "w")
f2.write(argfile2)
f2.close()

f3 = open("argfile3", "w")
f3.write(argfile3)
f3.close()
'''
