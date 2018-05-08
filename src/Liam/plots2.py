#Plots of the wavefunctions
import matplotlib.pyplot as plt
import numpy as np

A = np.loadtxt("QD/A.csv",delimiter=',')
V = np.loadtxt("QD/V.csv",delimiter=',')

dim = A.shape[0]
A_diag = np.zeros(dim)
for i in range(dim):
	A_diag[i] = A[i][i]

print(A_diag)


min = A_diag[0]
for i in range(dim):
	if A_diag[i] < min:
		min = A_diag[i]
		j = i

print("j: ",j,"min: ",min)

min2 = A_diag[0]
for i in range(dim):
	if min < A_diag[i] < min2:
		min2 = A_diag[i]
		j2 = i
print("j2: ",j2,"min2: ",min2)

min3 = A_diag[0]
for i in range(dim):
	if min2 < A_diag[i] < min3:
		min3 = A_diag[i]
		j3 = i
print("j3: ",j3,"min3: ",min3)

end = 500
rvals = np.linspace(0,15,end)
plt.plot(rvals,V[:end,j],label="$E_0 = $"+str(min))
plt.plot(rvals,V[:end,j2],label="$E_1 = $"+str(min2))
plt.plot(rvals,V[:end,j3],label="$E_2 = $"+str(min3))
plt.xlabel("r")
plt.ylabel("Wavefunction Amplitude")
plt.legend()
plt.savefig("plots/QDAmplitude.pdf")
plt.show()

pi = 3.1415926535897932384626433832975
for i in range(dim):
	V[i,j] *= 4*pi*V[i,j]
	V[i,j2] *= 4*pi*V[i,j2]
	V[i,j3] *= 4*pi*V[i,j3]

plt.plot(rvals,V[:end,j],label="$E_0 = $"+str(min))
plt.plot(rvals,V[:end,j2],label="$E_1 = $"+str(min2))
plt.plot(rvals,V[:end,j3],label="$E_2 = $"+str(min3))
plt.legend()
plt.xlabel("r")
plt.ylabel("Probability Density")
plt.savefig("plots/QDDensity.pdf")
plt.show()
