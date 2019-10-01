import numpy as np
import matplotlib.pyplot as plt 

def readfile(filename):

	infile = open(filename, 'r')
	a = infile.readline().split()
	
	n = int(a[1])
	rhomax = int(a[3])
	omega = float(a[5])
	
	res = np.zeros(n)
	for i, line in enumerate(infile):
		res[i] = float(line)
		

	return res, n, rhomax, omega
	

mean_list_main = []
rhomax_list_main= []
n_list_main= []

omega_list = ['0.010000', '0.500000', '1.000000', '5.000000']

n_ = 10
rhomax = "15"
omega = "0.000000"

mean_ = []
n_list = []
r = np.linspace(0,float(rhomax), n_)


res, n, rhomax, omega = readfile("eigvec"+str(n_)+(rhomax)+(omega)+".txt")


	
plt.plot(r, res, label=("$\omega_r$= %3.2f" %( omega)))
plt.title("Displacement of a buckling beam under a applied force")
plt.ylabel("displacement [dimentionless]")
plt.xlabel("")
plt.legend()

plt.show()


	
