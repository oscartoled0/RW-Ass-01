# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import BEMmodel as BEM

""" ---- User defined variables ---- """
# mode = 'cosinus' # 2 possible modes: constant and cosinus distribution of annulus. 
                  # Write 'cosinus' or 'constant' if you want to change the distribution.
plotter = False # If you don't want to plot all the default plots please enter False
                  
N = np.array([2, 4, 8, 16, 32, 64, 128])
TSR = 8

# Constant mode
CT = np.zeros(N.size)
CP = np.zeros(N.size)
for i in range(N.size):
    [results, CT[i], CP[i], Uinf, Radius] = BEM.ExecuteBEM(N[i], plotter, 'constant', TSR)
    
# Cosinus mode    
CT_cos = np.zeros(N.size)
CP_cos = np.zeros(N.size)
for i in range(N.size):
    [results, CT_cos[i], CP_cos[i], Uinf, Radius] = BEM.ExecuteBEM(N[i], plotter, 'cosinus', TSR)
    
plt.figure(figsize=(12, 6))
plt.plot(N, CT_cos, 'r--', label='Cosine distribution')
plt.plot(N, CT, 'b--', label = 'Uniform distribution')
plt.grid()
plt.xlabel('N (number of annulus)', fontsize=20)
plt.ylabel(r'$C_T$',fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Convergence.pdf', format='pdf', dpi=1000)
