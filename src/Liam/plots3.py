#plots3.py

import matplotlib.pyplot as plt

dim = [10, 20, 50, 100, 200, 500]
classical = [20, 146, 1829, 8401, 36040, 246585]
cyclic = [100, 800, 7500, 40000, 180000, 1500000]

plt.loglog(dim, classical)
plt.loglog(dim, cyclic)
plt.ylabel("Rotations")
plt.xlabel("Dimension")
plt.savefig("plots/Operations.pdf")
plt.show()

classical = [0.000247771, 0.001454993, 0.012223932, 0.111000046, 1.481069326, 62.638318191]
cyclic = [0.000314408, 0.001396732, 0.011375362, 0.107470026, 1.229996526, 72.663446155]

plt.loglog(dim, classical)
plt.loglog(dim, cyclic)
plt.ylabel("Seconds")
plt.xlabel("Dimension")
plt.savefig("plots/Timing.pdf")
plt.show()

dim = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
house = [91531,
479381,
10830105,
52730033,
145743628,
920983899,
4319594902,
9018110786,
137441120730]

plt.loglog(dim, house)
plt.ylabel("Time (ns)")
plt.xlabel("Dimension")
plt.savefig("plots/House.pdf")
plt.show()
