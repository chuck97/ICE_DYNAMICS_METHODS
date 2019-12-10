import numpy as np
import matplotlib.pyplot as plt
  
resudal1 = np.genfromtxt('./resudal/resud00001.txt')
resudal2 = np.genfromtxt('./resudal/resud00002.txt')
resudal3 = np.genfromtxt('./resudal/resud00003.txt')
resudal4 = np.genfromtxt('./resudal/resud00004.txt')
resudal5 = np.genfromtxt('./resudal/resud00005.txt')
resudal6 = np.genfromtxt('./resudal/resud00006.txt')
resudal7 = np.genfromtxt('./resudal/resud00007.txt')
resudal8 = np.genfromtxt('./resudal/resud00008.txt')
resudal9 = np.genfromtxt('./resudal/resud00009.txt')
resudal10 = np.genfromtxt('./resudal/resud00010.txt')

plt.plot(resudal1)
plt.plot(resudal2)
plt.plot(resudal3)
plt.plot(resudal4)
plt.plot(resudal5)
plt.plot(resudal6)
plt.plot(resudal7)
plt.plot(resudal8)
plt.plot(resudal9)
plt.plot(resudal10)

plt.show()

