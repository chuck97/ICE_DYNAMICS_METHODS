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
        
      
adr = str('./resudal/resud') + b + str('.txt')
   

resudal = np.genfromtxt(adr)

plt.plot(resudal)
plt.show()

