import matplotlib.pyplot as plt
import numpy as np

jacobi_convergence = [
[2,1],
[3,8],
[4,6],
[5,25],
[6,31],
[7,55],
[8,71],
[9,92],
[10,108],
[20,494],
[50, 3030],
[100, 11474],
[200, 43704],
[500, 249460],
[1000, 926120],
[2000, 3392385]]

jacobi_convergence = np.array(jacobi_convergence)


plt.rcParams["font.family"] = "serif"
plt.loglog(jacobi_convergence[:,0],jacobi_convergence[:,1],label="Jacobi")
plt.minorticks_on()
plt.grid(which='minor', linestyle='--', linewidth=0.5)
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.ylabel("Matrix Operations")
plt.xlabel("Matrix Size")

plt.legend()
plt.tight_layout()
plt.savefig("convergence.pdf")
plt.show()
