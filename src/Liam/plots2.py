#Plots of the wavefunctions
import matplotlib.pyplot as plt
import numpy as np

A = np.loadtxt("A.csv",delimiter=',')
V = np.loadtxt("V.csv",delimiter=',')

dim = A.shape[0]
A_diag = np.zeros(dim)
for i in range(dim):
	A_diag[i] = A[i][i]

A_diag *= dim/0.16
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

plt.plot(V[:,j])
plt.plot(V[:,j2])
plt.plot(V[:,j3])
plt.show()