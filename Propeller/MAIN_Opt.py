# import necessary libraries
#from platypus import NSGAII, Problem, Real
import matplotlib.pyplot as plt
import numpy as np
import BEM_Opt as BEM_Opt

""" ---- Code variables ---- """
N = 100 # annulus spacing

""" ---- Optimization ---- """

# BACKGROUND: This optimization is based on the concept of optimizing the performance
# of the rotor for maximum Cp with two approaches:
# Propeller case: Ct = Free choice and Maximum Efficiency at all radial positions
    
a_opt = 0.3
NB = 6
Radius = 0.7
Uinf = 60
Omega = 1200/60*2*np.pi
TSR = Omega*Radius/Uinf

Variables = [NB, Radius, Uinf, TSR, Omega]
    
finish = False  
CT0 = 4*a_opt*(1-a_opt) 
count=1
while (finish==False):
    
    [CP, CT, fnorm, ftan, gamma, chord, twist, r_R, phi, a, aline, vmag2] = BEM_Opt.ExecuteBEM_opt(N, Variables)
    count=count+1
    if (abs(CT0-CT)<0.002):
        
        finish = True
        
    else:
        
        Variables[0]=Variables[0] + 1.7*(CT0-CT)
        
        if (count/100-np.int(count/100) == 0):
            display(CT-CT0)
        
""" Plotter """

# Chord distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Chord distribution in Propeller', fontsize=20)
plt.plot(r_R, chord, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$c(r/R)$ [m]', fontsize=20)
plt.savefig('Chord_distrib_WT.pdf', format='pdf', dpi=1000)

# Twist distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Twist distribution in Propeller', fontsize=20)
plt.plot(r_R, twist, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\theta(r/R)$ [deg]', fontsize=20)
plt.savefig('Twist_distrib_WT.pdf', format='pdf', dpi=1000)

# Fnorm distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Normal force distribution in Propeller', fontsize=20)
plt.plot(r_R, fnorm, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$f_n(r/R)$ [N]', fontsize=20)
plt.savefig('fnorm_distrib_WT.pdf', format='pdf', dpi=1000)

# Circulation distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Circulation distribution in Propeller', fontsize=20)
plt.plot(r_R, gamma, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\Gamma(r/R)$ [1/s]', fontsize=20)
plt.savefig('Gamma_WT.pdf', format='pdf', dpi=1000)

# Inflow angle distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Inflow distribution in Propeller', fontsize=20)
plt.plot(r_R, a, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$a(r/R)$ [1/s]', fontsize=20)
plt.savefig('a_WT.pdf', format='pdf', dpi=1000)

    

