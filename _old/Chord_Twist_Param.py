import matplotlib.pyplot as plt
import numpy as np

# This function calculates the baseline coeffciients for the optimization of 
# chord and twist

c_R = np.loadtxt('c_R.txt')
twi = np.loadtxt('twist.txt')

N = 2;
# Curve fit: coefficients from higher to lower

c_R_fit = np.polyfit(c_R[:,0], c_R[:,1], N)
twi_fit = np.polyfit(twi[:,0], twi[:,1], N)

# Comparison between real and fitted data

#r_R = np.linspace(0,1,100)
#c_R_fitted = c_R_fit[0]*r_R**5 + c_R_fit[1]*r_R**4 + c_R_fit[2]*r_R**3 + \
#             c_R_fit[3]*r_R**2 + c_R_fit[4]*r_R + c_R_fit[5]
#             
#twi_fitted = twi_fit[0]*r_R**5 + twi_fit[1]*r_R**4 + twi_fit[2]*r_R**3 + \
#             twi_fit[3]*r_R**2 + twi_fit[4]*r_R + twi_fit[5]

r_R = np.linspace(0,1,100)
c_R_fitted = c_R_fit[0]*r_R**2 + c_R_fit[1]*r_R + c_R_fit[2]

twi_fitted = twi_fit[0]*r_R**2 + twi_fit[1]*r_R + twi_fit[2]
# c_R
fig1 = plt.figure(figsize=(12, 6))
plt.plot(c_R[:,0], c_R[:,1], 'k-', label=r'c_R Real')
plt.plot(twi[:,0], twi[:,1], 'b-', label=r'$\theta$ Real')
plt.plot(r_R, c_R_fitted, 'g--', label=r'c_R')
plt.plot(r_R, twi_fitted, 'g--', label=r'$\theta$')
plt.xlabel(r'r_R', fontsize=20)
plt.ylabel(r'c_R and $\theta$', fontsize=20)
plt.grid()
plt.legend(fontsize=20)