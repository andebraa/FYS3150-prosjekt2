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
rhomax = 15


mean_ = []
n_list = []
r = np.linspace(0,rhomax, n_)


for omega in omega_list:
	print(rhomax)
	res, n, rhomax, omega = readfile("eigvec"+str(n_)+str(rhomax)+(omega)+".txt")


	
	plt.plot(r, res, label=("$\omega_r$= %3.2f" %( omega)))
	plt.title("The wavefunction for the ground state for warying frequencies, and n = 100")
	plt.ylabel("$ \psi $")
	plt.xlabel("rho [dimentionless]")
	plt.legend()

plt.show()

print(mean_list_main)
print(n_list_main)
print(rhomax_list_main)


	
