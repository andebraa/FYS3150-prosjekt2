import numpy as np
import matplotlib.pyplot as plt 

def readfile(filename):

	infile = open(filename, 'r')
	a = infile.readline().split()
	
	n = int(a[1])
	rhomax = int(a[3])
	
	res = np.zeros(n)
	for i, line in enumerate(infile):
		res[i] = float(line)
		
	mean = np.mean(res)
	return mean, n, rhomax
	

mean_list_main = []
rhomax_list_main= []
n_list_main= []

n_ = 5
rhomax = 3
while n_ < 405:
	print(n_)
	mean_ = []
	n_list = []
	rhomax_list = []
	rhomax = 3
	while rhomax < 50:
		print(rhomax)
		mean, n, rhomax = readfile("res"+str(n_)+str(rhomax)+".txt")
		mean_.append(mean)
		n_list.append(n)
		rhomax_list.append(rhomax)
	
		rhomax += 5
	
	mean_list_main.append(mean_)	#hver liste i mean_list er en egen rekke med rho for fast n. 
	rhomax_list_main.append(rhomax_list)
	n_list_main.append(n_list)

	plt.plot(rhomax_list, mean_, label=("n = %3.1f" %n_))
	plt.title("Error in the eigenvalue solver algorithm for increasing Rho_max")
	plt.xlabel("rho_max [dimentionless]")
	plt.ylabel("Mean error for the n eigenvalues")
	plt.legend()
	n_ += 100
plt.show()

print(mean_list_main)
print(n_list_main)
print(rhomax_list_main)


	
