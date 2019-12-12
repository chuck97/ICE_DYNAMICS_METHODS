import numpy as np
import matplotlib.pyplot as plt

a = int(input("number of time step: "))

if ((a // 10) == 0):
	
	b = str('0000') + str(a)

else:
  
   if((a // 100) == 0):
	   
	   b = str('000') + str(a)	   
    
   else:
	   
	   if((a // 1000) == 0):
		   
		   b = str('00') + str(a)
		   
	   else:
		   
		   b = str('0') + str(a) 	   
        
      
adr_VP = str('./resudal/resud') + b + str('.txt')
adr_mEVP = ('/home/ibrae.ac.ru/s.s.petrov/Desktop/SEA_ICE/mEVP_method/resudal/resud') + b + str('.txt')   

resudal_VP = np.genfromtxt(adr_VP)
resudal_mEVP = np.genfromtxt(adr_mEVP)

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('k, local iteration number')
ax1.set_ylabel(r'$\frac{||res_k||_2}{||res_n||_2}, VP-Picard$', size=20, color=color)
ax1.plot(resudal_VP, label='VP-Picard(30) + Anderson', color=color)
ax1.tick_params(axis='y', labelcolor=color)


ax2 = ax1.twinx()

color = 'tab:blue'

ax2.set_ylabel(r'$\frac{||res_k||_2}{||res_n||_2}, mEVP(500)$', size=20, color=color) 
ax2.plot(resudal_mEVP, label='mEVP(500)', color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()


plt.show()

