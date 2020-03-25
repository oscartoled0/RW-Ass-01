# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import BEMmodel as BEM

""" ---- User defined variables ---- """
PR = True
WT = False
# mode = 'cosinus' # 2 possible modes: constant and cosinus distribution of annulus. 
                  # Write 'cosinus' or 'constant' if you want to change the distribution.
plotter = False # If you don't want to plot all the default plots please enter False
                  
N = np.array([10, 15, 100, 200, 300, 400])

# Constant mode
CT = np.zeros(N.size)
for i in range(N.size):
    [results, CT[i], CP, Uinf, Radius] = BEM.ExecuteBEM(PR, WT, N[i], plotter, mode = 'constant')
    
# Cosinus mode    
CT_cos = np.zeros(N.size)
for i in range(N.size):
    [results, CT_cos[i], CP, Uinf, Radius] = BEM.ExecuteBEM(PR, WT, N[i], plotter, mode = 'cosinus')
    
plt.figure(figsize=(12, 6))
plt.plot(N, CT_cos, 'r--', label='Cosine distribution')
plt.plot(N, CT, 'b--', label = 'Uniform distribution')
plt.grid()
plt.xlabel('N (number of annulus)', fontsize=20)
plt.ylabel(r'$C_T$',fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Convergence.pdf', format='pdf', dpi=1000)
